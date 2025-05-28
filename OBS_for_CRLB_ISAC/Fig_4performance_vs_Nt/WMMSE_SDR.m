clc
clear





K=4;
M=2;



% num_sensing_streams=max([0,  M-K]);

tolerance=1e-4;

I_in=3;
I_out=100;


SR_all=zeros(I_out,I_in);
CRB_all=zeros(I_out,I_in);
Time_all=zeros(I_out,I_in);


parfor k_par=1:I_out*I_in
    

channel=mod((k_par-1),I_out)+1;
weight=floor((k_par-1)/I_out)+1;

% progress_bar=((channel-1)*I_out+weight)/(I_out*length(delta_all))*100;
disp([channel,weight])


delta_s=1;
delta_c=0.25;

Nth=4;
Ntv=2^(weight+1);
Nt=Nth*Ntv;

Nrh=5;
Nrv=2^(weight+1);
Nr=Nrh*Nrv;

num_sensing_streams=Nt;

Pt=db2pow(10-30); %dBm
noise_c=db2pow(0-30); %dBm
noise_s=db2pow(0-30); %dBm



 % alpha=0.1*exp(1j*2*pi/3)*ones(M,1);


L=30;
% generate channels
rng(1,'twister')
alpha=0.1*(1+0.2*randn(M,1)).*exp(1j*2*pi*rand(M,1));
theta=-pi/3+2*pi/3*rand(M,1);
phi=-pi/3+2*pi/3*rand(M,1);

H_all=1/sqrt(2)*(randn(I_out,Nt,K)+1j*randn(I_out,Nt,K));

% theta=[-pi/6, pi/6, pi/10];
% phi=[-pi/3, pi/6, pi/10];
 

kappa=2*L/noise_s;

H=squeeze(H_all(channel,:,:));

tic

[A,dAtheta,dAphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nth,Ntv);
[B,dBtheta,dBphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nrh,Nrv);

U=diag(alpha);

Wc=delta_c*H./vecnorm(H);
Ws=delta_s*(randn(Nt,num_sensing_streams)+1j * randn(Nt,num_sensing_streams));

W=[Wc,Ws];
W=W*sqrt(Pt/trace(W*W'));

Rs=Ws*Ws';

W_last=W;


E=eye(4*M);

for count=1:2000

T_k=sum(square_abs(H'*Wc),2)+noise_c*ones(K,1)+diag(H'*Rs*H);
alpha_k=T_k./(T_k-square_abs(diag(H'*Wc)))-1;
beta_k=sqrt(1+alpha_k).*diag(H'*Wc)./T_k;
Sigma1=diag(sqrt(1+alpha_k).*beta_k);
Sigma2=diag(square_abs(beta_k));


C1=H*Sigma1;
C2=H*Sigma2*H';

[Wc,Rs,Rx]=update_W(Pt,delta_s,delta_c,Nt,K,M,num_sensing_streams,E,H,L,noise_s,A,dAtheta,dAphi,B,dBtheta,dBphi,U,C1,C2);

 Ws = chol(Rs, 'lower');
 W=[Wc,Ws]; 
 W=W*sqrt(Pt/trace(W*W'));
 T_k=sum(square_abs(H'*W),2)+noise_c*ones(K,1);
 rate_p=log(T_k./(T_k-square_abs(diag(H'*W))));

 FIM=calculateFIM(L,noise_s,Rx,A,dAtheta,dAphi,B,dBtheta,dBphi,U);


 obj=delta_c*sum(log(T_k./(T_k-square_abs(diag(H'*W)))))-delta_s*trace(inv(FIM));
%  Con=[Con;[sum(log(T_k./(T_k-square_abs(diag(H'*W))))),-trace(inv(FIM)),obj]];

 

 if norm(W-W_last)<tolerance
     break
 else
     W_last=W;
 end

end


Rx=W*W';
FIM=calculateFIM(L,noise_s,Rx,A,dAtheta,dAphi,B,dBtheta,dBphi,U);


T_k=sum(square_abs(H'*W),2)+noise_c*ones(K,1);
SR=sum(log(T_k./(T_k-square_abs(diag(H'*W)))));

CRB_all(k_par)=trace(inv(FIM));
SR_all(k_par)=SR;
Time_all(k_par)=toc;
end

save data_WMMSESDR_Nt
% Per_all=Per+Per_all;
% end

figure(1)
plot(mean(CRB_all),mean(SR_all))


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

function [Wc,Rs,Rx]=update_W(Pt,delta_s,delta_c,Nt,K,M,num_sensing_streams,E,H,L,noise_s,A,dAtheta,dAphi,B,dBtheta,dBphi,U,C1,C2)
cvx_clear
cvx_begin sdp quiet
variable Wc(Nt,K) complex 
variable Ws(Nt,num_sensing_streams) complex
variable t(4*M, 1)

variable Rc(Nt,Nt) complex hermitian
variable Rs(Nt,Nt) complex hermitian

Rx=Rc+Rs;
obj=cvx(0);
for k=1:K
    obj=obj+real(H(:,k)'*Rs*H(:,k));
end
obj=delta_s*sum(t)-delta_c*(2*real(trace(Wc'*C1))-quad_form(Wc(:),kron(eye(K),C2) )-obj  );
minimize(obj)

FIM= calculateFIM(L,noise_s,Rx,A,dAtheta,dAphi,B,dBtheta,dBphi,U);

subject to

trace(Rx)<=Pt;

[Rc, Wc
    Wc', eye(K)]==hermitian_semidefinite(Nt+K);

[Rs, Ws
    Ws', eye(num_sensing_streams)]==hermitian_semidefinite(Nt+num_sensing_streams);

for m=1:4*M
    [FIM, E(:,m)
        E(:,m)', t(m)]==hermitian_semidefinite(4*M+1);
end

cvx_end
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


