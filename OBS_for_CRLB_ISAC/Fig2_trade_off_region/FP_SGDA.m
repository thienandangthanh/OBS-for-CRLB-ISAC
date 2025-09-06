clc
clear



gamma_W=10e-6;
gamma_Y=0.02;
gamma_W_hat=0.02;
p_W=1;

%transmit power and noise
Pt=db2pow(10-30); %dBm
noise_c=db2pow(0-30); %dBm
noise_s=db2pow(0-30); %dBm
L=30;
kappa=2*L/noise_s;

K=4;
M=2;

%antenna for transmit
Nth=4;
Ntv=4;
Nt=Nth*Ntv;

%antenna for receiver
Nrh=5;
Nrv=4;
Nr=Nrh*Nrv;

tolerance=1e-5;
delta_all=-7:0.2:4.8;

I_in=length(delta_all);
I_out=50;

SR_all=zeros(I_out,I_in);
CRB_all=zeros(I_out,I_in);
Time_all=zeros(I_out,I_in);
Obj_all=zeros(I_out,I_in);


for k_par=1 : I_out*I_in

    channel=mod((k_par-1),I_out)+1;
    weight=floor((k_par-1)/I_out)+1;

    disp([channel,weight])

    delta_s=1;
    delta_c= 10^(delta_all(weight));

    if delta_c>1
        delta_s=1/delta_c;
        delta_c=1;
    end



    num_sensing_streams=Nt;
    rng(1,'twister')
    alpha=0.1*(1+0.2*randn(M,1)).*exp(1j*2*pi*rand(M,1));
    
    theta=-pi/3+2*pi/3*rand(M,1);
    phi=-pi/3+2*pi/3*rand(M,1);

    H_all=1/sqrt(2)*(randn(I_out,Nt,K)+1j*randn(I_out,Nt,K));

    H=squeeze(H_all(channel,:,1:K));


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
Y=inv(FIM);
W_last=W;
W_hat=W;
Con=[];

for count=1:2000

% update alpha, beta with closed form
T_k=sum(square_abs(H'*W(:,1:K)),2)+noise_c*ones(K,1);
alpha_k=T_k./(T_k-square_abs(diag(H'*W(:,1:K))))-1;
beta_k=sqrt(1+alpha_k).*diag(H'*W(:,1:K))./T_k;

% calculate the gradient of W
Sigma1=diag(sqrt(1+alpha_k).*beta_k);
Sigma2=diag(square_abs(beta_k));

CRBM=inv(FIM);

Q=construct_matrixQ(L,noise_s,CRBM*CRBM,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
C1=[delta_c*H*Sigma1,zeros(Nt,num_sensing_streams)];
C2=0.5*delta_s*(Q+Q')-delta_c*H*Sigma2*H';

% update W
gradient_W=2*(C1+C2*W);
W=W+gamma_W*gradient_W+gamma_W*p_W*(W-W_hat);
W=W*sqrt(Pt/trace(W*W'));

% update the gradient of Y and Y
gradient_Y=2*delta_s*(FIM*Y-Y);
Y=Y-gamma_Y*gradient_Y;

% update W_hat
W_hat=W_hat+gamma_W_hat*(W-W_hat);


% calculate the objective functions
FIM=calculateFIM(L,noise_s,W,A,dAtheta,dAphi,B,dBtheta,dBphi,U);

T_k=sum(square_abs(H'*W),2)+noise_c*ones(K,1);

rate_p=log(T_k./(T_k-square_abs(diag(H'*W))));

 obj=delta_c*sum(log(T_k./(T_k-square_abs(diag(H'*W)))))-delta_s*trace(inv(FIM));
 Con=[Con;[sum(log(T_k./(T_k-square_abs(diag(H'*W))))),-trace(inv(FIM)),obj]];

 if norm(W-W_last)<tolerance
     break
 else
     W_last=W;
 end

end

 Rx=W*W';
    FIM=calculateFIM(L,noise_s,W,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
    CRB=inv(FIM);
    SR=sum(log(T_k./(T_k-square_abs(diag(H'*W(:,1:K))))));

    CRB_all(k_par)=trace(inv(FIM));
    SR_all(k_par)=SR;
    Time_all(k_par)=toc;
    Obj_all(k_par)=delta_c*SR-delta_s*trace(inv(FIM));

end

save data_FP_SGDA

% 
% figure(1)
% hold on
% plot(mean(Obj_all))
% plot(mean(SR_all))
% plot(mean(CRB_all))
% legend('OBJ', 'SR', 'CRB')
% xlabel('Number of Communications Users')
% ylabel('Total Objective Value')
% 
% figure(2)
% plot(mean(Time_all))
% xlabel('Number of Communications Users')
% ylabel('Total Objective Value')



% 
toc
figure(1)
plot(Con(:,1))
xlabel('iteration')
ylabel('Sum Rate (nat/s/Hz)')
figure(2)
plot(Con(:,2))
xlabel('iteration')
ylabel('trace of inverse of FIM')
figure(3)
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


