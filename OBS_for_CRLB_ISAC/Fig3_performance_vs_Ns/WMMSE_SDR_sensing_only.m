clc
clear

K_max=4;

Nth=4;
Ntv=4;
Nt=Nth*Ntv;

Nrh=5;
Nrv=4;
Nr=Nrh*Nrv;

M_max=3;
m_target=2;
M=m_target;
num_sensing_streams=0;

tolerance=1e-4;

I_in=1;
I_out=100;



SR_all=zeros(length(Ns_all),I_out*I_in);
CRB_all=zeros(length(Ns_all),I_out*I_in);
Time_all=zeros(length(Ns_all),I_out*I_in);


parfor k_par=1:I_out*I_in
    

channel=mod((k_par-1),I_out)+1;
weight=floor((k_par-1)/I_out)+1;


disp([channel,weight])


delta_s=1;
delta_c=0;

K=weight;
Pt=db2pow(20-30); %dBm
noise_c=db2pow(0-30); %dBm
noise_s=db2pow(0-30); %dBm


L=128; kappa=2*L/noise_s;

% generate channels
 rng(1,'twister')
 alpha_all=0.1*(1+0.2*randn(M_max,1)).*exp(1j*2*pi*rand(M_max,1));
 theta_all=-pi/3+2*pi/3*rand(M_max,1);
 phi_all=-pi/3+2*pi/3*rand(M_max,1);

 H_all=1/sqrt(2)*(randn(Nt,K_max, I_out)+1j*randn(Nt,K_max,I_out));

 H=squeeze(H_all(:,1:K,channel));

 alpha=alpha_all(1:m_target);
 theta=theta_all(1:m_target);
 phi=phi_all(1:m_target);



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
Con=[];
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
 % T_k=sum(square_abs(H'*W),2)+noise_c*ones(K,1);
 % rate_p=log(T_k./(T_k-square_abs(diag(H'*W(:,1:K)))));

 FIM=calculateFIM(L,noise_s,Rx,A,dAtheta,dAphi,B,dBtheta,dBphi,U);


 obj=delta_c*sum(log(T_k./(T_k-square_abs(diag(H'*W(:,1:K))))))-delta_s*trace(inv(FIM));
 Con=[Con;[sum(log(T_k./(T_k-square_abs(diag(H'*W(:,1:K)))))),-trace(inv(FIM)),obj]];

 
 if norm(W-W_last)<tolerance
     break
 else
     W_last=W;
 end

end


Rx=Rx-Wc*Wc';
[Eigenvector,eigenvalues]=eig(Rx);
[eigenvalues, index]=sort(diag(eigenvalues),'descend');
Eigenvector=Eigenvector(:,index);
CRB_temp=zeros(9,1);
SR_temp=zeros(9,1);
for m_stream=0:8
    Ws=Eigenvector(:,1:m_stream)*sqrt(diag(eigenvalues(1:m_stream)));
    W=[Wc,Ws];
    W=W*sqrt(Pt/trace(W*W'));
    T_k=sum(square_abs(H'*W),2)+noise_c*ones(K,1);
    rate_p=log(T_k./(T_k-square_abs(diag(H'*W(:,1:K)))));

    FIM=calculateFIM(L,noise_s,W*W',A,dAtheta,dAphi,B,dBtheta,dBphi,U);
    CRB_temp(m_stream+1)=trace(inv(FIM));
    SR_temp(m_stream+1)=sum(log(T_k./(T_k-square_abs(diag(H'*W(:,1:K))))));
    
end

CRB_all(:, k_par)=CRB_temp;
SR_all(:,k_par)=SR_temp;
end

save data_WMMSE_SDR_sensing


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


