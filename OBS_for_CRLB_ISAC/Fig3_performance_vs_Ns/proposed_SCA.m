clc
clear

%antenna for transmit
Nth=4;
Ntv=4;
Nt=Nth*Ntv;

%antenna for receiver
Nrh=5;
Nrv=4;
Nr=Nrh*Nrv;

% number of users and targets
K_max=4;
M_max=3;
m_target=3;

%transmit power and noise
Pt=db2pow(10-30); %dBm
noise_c=db2pow(0-30); %dBm
noise_s=db2pow(0-30); %dBm
L=128;
kappa=2*L/noise_s;

tolerance=1e-4;

Ns_all=[0,1,2,3,4,5,6,7,8];

I_in=length(Ns_all);
I_out=50;

SR_all=zeros(I_out,I_in);
CRB_all=zeros(I_out,I_in);
Time_all=zeros(I_out,I_in);
Obj_all=zeros(3,I_out,I_in);

for user=2:4

    for k_par=1 :I_out*I_in

        channel=mod((k_par-1),I_out)+1;
        weight=floor((k_par-1)/I_out)+1;

        disp([channel,weight])
        delta_s=1;
        delta_c=0.25;
        % 10^(delta_all(weight));

        K=user-1;
        num_sensing_streams=weight-1;
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
        % equivalent_channel=B*U*A';
        % [eigvector,~]=eigs(equivalent_channel'*equivalent_channel,1,'largestabs');

        % Wc=delta_c*H./vecnorm(H);
        Wc=delta_c*H./vecnorm(H);

        %Wc=randn(Nt,K)+1j*randn(Nt,K);
        % Ws=delta_s*repmat(eigvector,1,num_sensing_streams);
        Ws=randn(Nt,num_sensing_streams)+1j*randn(Nt,num_sensing_streams);

        % Ws=initial_Ws(L,noise_s,Nt,num_sensing_streams,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
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
            % rank(Q+Q')
            C1=[delta_c*H*Sigma1,zeros(Nt,num_sensing_streams)];
            C2=0.5*delta_s*(Q+Q')-delta_c*H*Sigma2*H';

            %dom eigvalue

            mu=abs(eigs(H*Sigma2*H',1,'LM'));
            C2=delta_c*mu*eye(Nt)+C2;

            for iter=1
                W=C1+C2*W;

                W=W*sqrt(Pt/trace(W*W'));
            end

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
        FIM=calculateFIM(L,noise_s,W,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
        CRB=inv(FIM);
        SR=sum(log(T_k./(T_k-square_abs(diag(H'*W(:,1:K))))));

        CRB_all(k_par)=trace(inv(FIM));
        SR_all(k_par)=SR;
        Time_all(k_par)=toc;
        Obj_all(user,channel,weight)=delta_c*SR-delta_s*trace(inv(FIM));

    end

    save data Obj_all

    % figure(1)
    % hold on
    % plot(Ns_all,mean(Obj_all))
    % save data_PSLA_twolayer
    % plot(Per_all(1,:)/channel_number,Per_all(2,:)/channel_number)
    % figure(m_target)
    % hold off
    % plot(CRB_all(1,:),SR_all(1,:))

end

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

