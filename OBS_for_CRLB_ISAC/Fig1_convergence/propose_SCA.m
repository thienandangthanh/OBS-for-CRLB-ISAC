clc
clear

Nth=4;
Ntv=4;
Nt=Nth*Ntv;

Nrh=5;
Nrv=4;
Nr=Nrh*Nrv;

channel_number=50;
K=4;
M=2;

rng(0)
theta=-pi/3+2*pi/3*rand(M,1);
phi=-pi/3+2*pi/3*rand(M,1);

num_sensing_streams=Nt;
tolerance=1e-5;

delta_all=-7:0.05:4.8;

Per_all=zeros(2,length(delta_all));
Con=cell(3,1);

for channel=2
for de=1:3
lin=[0.05,0.1,0.15];
disp([channel,de])
delta_s=1;
delta_c=lin(de);
%10^(delta_all(de));


Pt=db2pow(10-30); %dBm
noise_c=db2pow(0-30); %dBm
noise_s=db2pow(0-30); %dBm



rng(channel)

alpha=0.1*(1+0.2*randn(M,1)).*exp(1j*2*pi*rand(M,1));


L=30;

% theta=[-pi/6, pi/6, pi/10];
% phi=[-pi/3, pi/6, pi/10];
 

kappa=2*L/noise_s;

H=1/sqrt(2)*(randn(Nt,K)+1j*randn(Nt,K));


tic
[A,dAtheta,dAphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nth,Ntv);
[B,dBtheta,dBphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nrh,Nrv);

U=diag(alpha);

Wc=randn(Nt,K)+1j*randn(Nt,K);
% Ws=delta_s*repmat(eigvector,1,num_sensing_streams);
 Ws=randn(Nt,num_sensing_streams)+1j*randn(Nt,num_sensing_streams);

% Ws=initial_Ws(L,noise_s,Nt,num_sensing_streams,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
W=[Wc,Ws];


W=W*sqrt(Pt/trace(W*W'));

FIM=calculateFIM(L,noise_s,W,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
W_last=W;


for count=1:2000

T_k=sum(square_abs(H'*W(:,1:K)),2)+noise_c*ones(K,1);
alpha_k=T_k./(T_k-square_abs(diag(H'*W(:,1:K))))-1;
beta_k=sqrt(1+alpha_k).*diag(H'*W(:,1:K))./T_k;
Sigma1=diag(sqrt(1+alpha_k).*beta_k);
Sigma2=diag(square_abs(beta_k));

CRBM=inv(FIM);

% PHI=CRBM*CRBM;
Q=construct_matrixQ(L,noise_s,CRBM*CRBM,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
C1=[delta_c*H*Sigma1,zeros(Nt,num_sensing_streams)];
C2=-0.5*delta_s*(Q+Q')+delta_c*H*Sigma2*H';

%dom eigvalue
mu=abs(eigs(C2,1,'LM'));
C2=mu*eye(Nt)-C2;
for iter=1
W=C1+C2*W;

% gradient
% mu=5e-5;
% W= W + mu *(C1-C2*W );

W=W*sqrt(Pt/trace(W*W'));
end

FIM=calculateFIM(L,noise_s,W,A,dAtheta,dAphi,B,dBtheta,dBphi,U);

 T_k=sum(square_abs(H'*W),2)+noise_c*ones(K,1);

 rate_p=log(T_k./(T_k-square_abs(diag(H'*W))));
 % CRBM=delta_s*trace(inv(FIM));
 obj=delta_c*sum(log(T_k./(T_k-square_abs(diag(H'*W)))))-delta_s*trace(inv(FIM));
 Con{de}=[Con{de};[sum(log(T_k./(T_k-square_abs(diag(H'*W))))),-trace(inv(FIM)),obj]];

 if norm(W-W_last)<tolerance
     break
 else
     W_last=W;
 end

end


Rx=W*W';
CRB=trace(inv(FIM));
SR=sum(log(T_k./(T_k-square_abs(diag(H'*W)))));

Per(:,de)=[CRB;SR];
toc
figure(1)
hold on
plot(Con{de}(:,1))
xlabel('iteration')
ylabel('Sum Rate (nat/s/Hz)')
figure(2)
hold on
plot(Con{de}(:,2))
xlabel('iteration')
ylabel('trace of inverse of FIM')
figure(3)
hold on
plot(Con{de}(:,3))
xlabel('iteration')
ylabel('total objective value')

end

% Per_all=Per+Per_all;
end

% save data
% plot(Per_all(1,:)/channel_number,Per_all(2,:)/channel_number)
% plot(Per(1,:),Per(2,:))
% 

save data_convergence Con


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
