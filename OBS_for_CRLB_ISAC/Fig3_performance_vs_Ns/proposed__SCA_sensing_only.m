clc
clear


Nth=4;
Ntv=4;
Nt=Nth*Ntv;

Nrh=5;
Nrv=4;
Nr=Nrh*Nrv;
I_out=1;
L=128;
K_max=4;
M_max=3;
m_target=3;
stream_test=zeros(6,1);

for m_stream=1:8
num_sensing_streams=m_stream;

tolerance=1e-5;

delta_all=[-7:0.05:4.8];

Per_all=zeros(2,length(delta_all));
tic
for channel=1
for de=1
    % :length(delta_all)

delta_s=1;
delta_c=0;
% 10^(delta_all(de));
K=2;

Pt=db2pow(20-30); %dBm
noise_c=db2pow(0-30); %dBm
noise_s=db2pow(0-30); %dBm



rng(1,'twister')
alpha_all=0.1*(1+0.2*randn(M_max,1)).*exp(1j*2*pi*rand(M_max,1));
theta_all=-pi/3+2*pi/3*rand(M_max,1);
phi_all=-pi/3+2*pi/3*rand(M_max,1);

H_all=1/sqrt(2)*(randn(Nt,K_max, I_out)+1j*randn(Nt,K_max,I_out));

H=squeeze(H_all(:,1:K,channel));
alpha=alpha_all(1:m_target);
theta=theta_all(1:m_target);
phi=phi_all(1:m_target);


kappa=2*L/noise_s;


% Q=zeros(Nt,Nt,K);
% for k=1:K
%     Q(:,:,k)=H(:,k)*H(:,k)';
% end

[A,dAtheta,dAphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nth,Ntv);
[B,dBtheta,dBphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nrh,Nrv);

U=diag(alpha);
equivalent_channel=B*U*A';
[eigvector,~]=eigs(equivalent_channel'*equivalent_channel,1,'largestabs');

Wc=delta_c*H./vecnorm(H);

Ws=randn(Nt,num_sensing_streams)+1j*randn(Nt,num_sensing_streams);
W=[Wc,Ws];
W=W*sqrt(Pt/trace(W*W'));
FIM=calculateFIM(L,noise_s,W,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
W_last=W;
Con=[];

for count=1:2000

T_k=sum(square_abs(H'*W(:,1:K)),2)+noise_c*ones(K,1);
alpha_k=T_k./(T_k-square_abs(diag(H'*W(:,1:K))))-1;
beta_k=sqrt(1+alpha_k).*diag(H'*W(:,1:K))./T_k;
Sigma1=diag(sqrt(1+alpha_k).*beta_k);
Sigma2=diag(square_abs(beta_k));

CRBM=inv(FIM);

% PHI=CRBM*CRBM;
Q=construct_matrixQ(L,noise_s,CRBM*CRBM,A,dAtheta,dAphi,B,dBtheta,dBphi,U);

C2=-0.5*delta_s*(Q+Q')+delta_c*H*Sigma2*H';

mu=abs(eigs(C2,1,'LM'));

C1=[delta_c*H*Sigma1,zeros(Nt,num_sensing_streams)];
C2=mu*eye(Nt)-C2;


W=C1+C2*W;
W=W*sqrt(Pt/trace(W*W'));


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


Rx=W*W';
CRB=trace(inv(FIM))
SR=sum(log(T_k./(T_k-square_abs(diag(H'*W)))));

Per(:,de)=[CRB;SR];

stream_test(m_stream)=obj;

end

% Per_all=Per+Per_all;
end

end
% plot(Per(1,:),Per(2,:))

save data_sensing.mat

% toc
% figure(1)
% plot(Con(:,1))
% xlabel('iteration')
% ylabel('Sum Rate (nat/s/Hz)')
% figure(2)
% plot(Con(:,2))
% xlabel('iteration')
% ylabel('trace of inverse of FIM')
% figure(3)
% plot(Con(:,3))
% xlabel('iteration')
% ylabel('total objective value')


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


