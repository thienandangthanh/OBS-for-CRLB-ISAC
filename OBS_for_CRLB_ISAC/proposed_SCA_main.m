clc
clear

%% Introaduction
% This script provides a basic implementation of the simulation code for our paper. 
% 
% - By adjusting the weight coefficients, we obtain:
%   - Figure 2: Convergence behavior
%   - Figure 3: Tradeoff region
%
% - By modifying the number of sensing streams, we obtain Figure 4.
% - By changing the number of transmit and receive antennas, we obtain Figure 5 and Table 1.
% - By varying the number of users, we obtain Figure 6.
% - By adjusting the transmit power, we obtain Figure 7.



%% system parameters
%number of transmit antennas
Nth=4;
Ntv=4;
Nt=Nth*Ntv;

%number of receive antennas
Nrh=5;
Nrv=4;
Nr=Nrh*Nrv;

% number of communication user and sensing target
K=4;
M=2;



% number of radar signals
num_sensing_streams=Nt;

% tolerance
tolerance=1e-5;

% matrix to hold the convergence behavior
Con=[];

% time slot
L=64;

% weights for balance communications and sensing
delta_s=1;
delta_c=0.25;

% transmit power and noise
Pt=db2pow(10-30); %dBm
noise_c=db2pow(0-30); %dBm
noise_s=db2pow(0-30); %dBm

%% generate channels and angles to estimate, and initial point
% radar cross-section
rng(4)
alpha=0.1*(1+0.2*randn(M,1)).*exp(1j*2*pi*rand(M,1));

% angle to estimate
theta=-pi/3+2*pi/3*rand(M,1);
phi=-pi/3+2*pi/3*rand(M,1);


% communications channels
H=1/sqrt(2)*(randn(Nt,K)+1j*randn(Nt,K));

% construct steer vectors and the corresponding derivativer vectors
tic
[A,dAtheta,dAphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nth,Ntv);
[B,dBtheta,dBphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nrh,Nrv);

U=diag(alpha);

% generate initial feasiable points

% we use MRT directions for communication beamforming matrix
Wc=H;

% we use a random matrix for sensing beamforming matrix
Ws=randn(Nt,num_sensing_streams)+1j*randn(Nt,num_sensing_streams);

% alternatively, we can solve the CRLB minimization (i.e., let delta_c=0)
% to obtain an initial point.
% Ws=initial_Ws(L,noise_s,Nt,num_sensing_streams,A,dAtheta,dAphi,B,dBtheta,dBphi,U);

% normalize W
W=[Wc,Ws];
W=W*sqrt(Pt/trace(W*W'));


%% optimization begin
% calculate fisher information matrix according to eq(8)
FIM=calculateFIM(L,noise_s,W,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
W_last=W;

%main loop
for count=1:2000

% update auxilary variable according to (16), step 4 in Algorithm 1
T_k=sum(square_abs(H'*W(:,1:K)),2)+noise_c*ones(K,1);
alpha_k=T_k./(T_k-square_abs(diag(H'*W(:,1:K))))-1;
beta_k=sqrt(1+alpha_k).*diag(H'*W(:,1:K))./T_k;
Sigma1=diag(sqrt(1+alpha_k).*beta_k);
Sigma2=diag(square_abs(beta_k));

% calculate the inverse of the FIM, step 5 in Algorithm 1
CRBM=inv(FIM);

% PHI=CRBM*CRBM;and construct matrix Q according to eq(24)
Q=construct_matrixQ(L,noise_s,CRBM*CRBM,A,dAtheta,dAphi,B,dBtheta,dBphi,U);

% construct matrix C1 and C2
C1=[delta_c*H*Sigma1,zeros(Nt,num_sensing_streams)];
C2=-0.5*delta_s*(Q+Q')+delta_c*H*Sigma2*H';

%dom eigvalue
mu=abs(eigs(C2,1,'LM'));
C2=mu*eye(Nt)-C2;

% update W, step 6 in Algorithm 1
W=C1+C2*W;
W=W*sqrt(Pt/trace(W*W'));


% calculate the objective values
FIM=calculateFIM(L,noise_s,W,A,dAtheta,dAphi,B,dBtheta,dBphi,U);

 T_k=sum(square_abs(H'*W),2)+noise_c*ones(K,1);

 rate_p=log(T_k./(T_k-square_abs(diag(H'*W))));
 % CRBM=delta_s*trace(inv(FIM));

 obj=delta_c*sum(log(T_k./(T_k-square_abs(diag(H'*W)))))-delta_s*trace(inv(FIM));
 Con=[Con;[sum(log(T_k./(T_k-square_abs(diag(H'*W))))),-trace(inv(FIM)),obj]];

 if norm(W-W_last)<tolerance
     break
 else
     W_last=W;
 end

end


%% result check
% calculate final objective values
Rx=W*W';
CRB=trace(inv(FIM));
SR=sum(log(T_k./(T_k-square_abs(diag(H'*W)))));


% plot figures to check convergence 
toc
figure(1)
hold on
plot(Con(:,1))
xlabel('iteration')
ylabel('Sum Rate (nat/s/Hz)')
figure(2)
hold on
plot(Con(:,2))
xlabel('iteration')
ylabel('trace of inverse of FIM')
figure(3)
hold on
plot(Con(:,3))
xlabel('iteration')
ylabel('total objective value')










function Ws=initial_Ws(L,noise_s,Nt,num_sensing_streams,A,dAtheta,dAphi,B,dBtheta,dBphi,U)

Ws=randn(Nt,num_sensing_streams)+1j*randn(Nt,num_sensing_streams);


for iter=1:5
    FIM=calculateFIM(L,noise_s,Ws,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
    CRBM=inv(FIM);
    Q=construct_matrixQ(L,noise_s,CRBM*CRBM,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
Ws=(Q+Q')*Ws;
Ws=Ws*sqrt(1/trace(Ws*Ws'));
end


end

function Q=construct_matrixQ(L,noise_s,Phi,A,dAtheta,dAphi,B,dBtheta,dBphi,U)
M=size(U,1);

phi11=Phi(1:M,1:M);
phi12=Phi(1:M,M+1:2*M);
phi13=Phi(1:M,2*M+1:3*M);
phi14=Phi(1:M,3*M+1:4*M);
phi22=Phi(M+1:2*M,M+1:2*M);
phi23=Phi(M+1:2*M,2*M+1:3*M);
phi24=Phi(M+1:2*M,3*M+1:4*M);
phi33=Phi(2*M+1:3*M,2*M+1:3*M);
phi34=Phi(2*M+1:3*M,3*M+1:4*M);
phi44=Phi(3*M+1:4*M,3*M+1:4*M);
% Compute Q11

Q11 = (A * U' * (phi11 .* (dBtheta' * dBtheta)) * U * A')  + ...
      (dAtheta * U' * (phi11.* (B' * dBtheta)) * U * A')  + ...
      (A * U' * (phi11.* (dBtheta' * B)) * U * dAtheta' )  + ...
      (dAtheta * U' * (phi11.* (B' * B)) * U * dAtheta') ;

% Compute F12
Q12 = (A * U' * (2*phi12.* (dBtheta' * dBphi)) * U * A' )  + ...
    ( dAtheta * U' * (2*phi12.* (B' * dBphi)) * U * A')  + ...
      (A * U' * (2*phi12.* (dBtheta' * B)) * U * dAphi') + ...
      (dAtheta * U' * (2*phi12.* (B' * B)) * U * dAphi') ;

% Compute F13
Q13 = (A * U' * ((2*phi13+2j*phi14).* (dBtheta' * B)) * A' )  + ...
      (dAtheta * U' * ((2*phi13+2j*phi14).* (B' * B)) * A') ;

% Compute F22
Q22 = (A * U' * (phi22.* (dBphi' * dBphi)) * U * A' )  + ...
      (dAphi * U' * (phi22.* (B' * dBphi)) * U * A')  + ...
      (A * U' * (phi22.* (dBphi' * B)) * U * dAphi' )  + ...
      (dAphi * U' * (phi22.* (B' * B)) * U * dAphi' ) ;

% Compute F23
Q23 = (A * U' * ((2*phi23+2j*phi24).* (dBphi' * B)) *A')  + ...
      (dAphi * U' * ((2*phi23+2j*phi24).* (B' * B)) *A' ) ;

% Compute F33
Q33 = ( A * ((phi33+phi44+2j*phi34).* (B' * B))*A') ;

Q=2*L/noise_s*(Q11+Q12+Q13+Q22+Q23+Q33);
end
