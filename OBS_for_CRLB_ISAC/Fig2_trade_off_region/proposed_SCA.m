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
K=4;
M_max=2;

%transmit power and noise
Pt=db2pow(10-30); %dBm
noise_c=db2pow(0-30); %dBm
noise_s=db2pow(0-30); %dBm
L=30;
kappa=2*L/noise_s;

num_sensing_streams=Nt;
tolerance=1e-4;
delta_all=-7:0.2:4.8;

I_in=length(delta_all);
I_out=50;

% Initialize result storage
SR_all=zeros(I_out,I_in);
CRB_all=zeros(I_out,I_in);
Time_all=zeros(I_out,I_in);

% Total number of parameter combinations
total_combinations = I_out * I_in;

fprintf('Starting parameter sweep over %d combinations (%d x %d)\n', total_combinations, I_out, I_in);
fprintf('Delta range: [%.1f, %.1f] dB with %.1f dB steps\n', delta_all(1), delta_all(end), delta_all(2)-delta_all(1));

% Fixed number of targets for this analysis
m_target = 2;
M = m_target;

% Generate all random parameters once at the beginning for consistency
rng(1,'twister')
alpha_all=0.1*(1+0.2*randn(M_max,1)).*exp(1j*2*pi*rand(M_max,1));
theta_all=-pi/3+2*pi/3*rand(M_max,1);
phi_all=-pi/3+2*pi/3*rand(M_max,1);
H_all=1/sqrt(2)*(randn(I_out,Nt,K)+1j*randn(I_out,Nt,K));

% Extract target parameters
alpha=alpha_all(1:m_target);
theta=theta_all(1:m_target);
phi=phi_all(1:m_target);

% Construct steering matrices (same for all scenarios)
[A,dAtheta,dAphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nth,Ntv);
[B,dBtheta,dBphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Nrh,Nrv);
U=diag(alpha);

% Store convergence data for the last scenario (for plotting)
Con_last = [];

% Main parameter sweep loop
for k_par = 1:total_combinations

    % Display progress
    if mod(k_par-1, 100) == 0 || k_par == 1
        fprintf('Processing combination %d/%d (%.1f%%)\n', k_par, total_combinations, 100*k_par/total_combinations);
    end

    % Map linear index to 2D indices
    channel = mod((k_par-1), I_out) + 1;
    weight = floor((k_par-1) / I_out) + 1;

    % Set trade-off parameters
    delta_s = 1;
    delta_c = 10^(delta_all(weight));

    % Extract channel realization
    H = squeeze(H_all(channel,:,:));

    % Start timing for this scenario
    tic;

    % Initialize beamforming matrices with fresh random values for each scenario
    rng(k_par + 1000, 'twister'); % Different seed for each scenario
    Wc = randn(Nt,K) + 1j*randn(Nt,K);
    Ws = randn(Nt,num_sensing_streams) + 1j*randn(Nt,num_sensing_streams);

    % Combine and normalize beamforming matrices
    W = [Wc, Ws];
    W = W * sqrt(Pt/trace(W*W'));

    % Initial FIM calculation
    FIM = calculateFIM(L,noise_s,W,A,dAtheta,dAphi,B,dBtheta,dBphi,U);
    W_last = W;
    Con = [];

    % SCA iteration loop
    for count = 1:2000

        % Update auxiliary variables for communication
        T_k = sum(square_abs(H'*W(:,1:K)),2) + noise_c*ones(K,1);
        alpha_k = T_k./(T_k-square_abs(diag(H'*W(:,1:K)))) - 1;
        beta_k = sqrt(1+alpha_k).*diag(H'*W(:,1:K))./T_k;
        Sigma1 = diag(sqrt(1+alpha_k).*beta_k);
        Sigma2 = diag(square_abs(beta_k));

        % Update FIM and construct matrix Q for sensing
        CRBM = inv(FIM);
        Q = construct_matrixQ(L,noise_s,CRBM*CRBM,A,dAtheta,dAphi,B,dBtheta,dBphi,U);

        % Construct matrices C1 and C2 for SCA update
        C1 = [delta_c*H*Sigma1, zeros(Nt,num_sensing_streams)];
        C2 = 0.5*delta_s*(Q+Q') - delta_c*H*Sigma2*H';

        % Dominant eigenvalue calculation (Type B construction)
        mu = abs(eigs(H*Sigma2*H',1,'LM'));
        C2 = delta_c*mu*eye(Nt) + C2;

        % SCA update step
        W = C1 + C2*W;
        W = W * sqrt(Pt/trace(W*W'));

        % Update FIM
        FIM = calculateFIM(L,noise_s,W,A,dAtheta,dAphi,B,dBtheta,dBphi,U);

        % Calculate performance metrics for convergence tracking
        T_k = sum(square_abs(H'*W),2) + noise_c*ones(K,1);
        obj = delta_c*sum(log(T_k./(T_k-square_abs(diag(H'*W))))) - delta_s*trace(inv(FIM));

        % Store convergence data
        sum_rate = sum(log(T_k./(T_k-square_abs(diag(H'*W)))));
        trace_inv_fim = -trace(inv(FIM));
        Con = [Con; [sum_rate, trace_inv_fim, obj]];

        % Check convergence
        if norm(W-W_last) < tolerance
            break
        else
            W_last = W;
        end
    end

    % Final calculations for this scenario
    T_k = sum(square_abs(H'*W),2) + noise_c*ones(K,1);
    SR = sum(log(T_k./(T_k-square_abs(diag(H'*W)))));
    CRB_trace = trace(inv(FIM));
    computation_time = toc;

    % Store results in proper matrix locations
    SR_all(channel, weight) = SR;
    CRB_all(channel, weight) = CRB_trace;
    Time_all(channel, weight) = computation_time;

    % Store convergence data for the last scenario (for plotting)
    if k_par == total_combinations
        Con_last = Con;
    end
end

fprintf('Parameter sweep completed!\n');
fprintf('Average computation time per scenario: %.3f seconds\n', mean(Time_all(:)));
fprintf('Total computation time: %.1f minutes\n', sum(Time_all(:))/60);

% Save to the specified file
save('data_proposed_SCA.mat');
fprintf('Results saved to data_proposed_SCA.mat\n');

% Display summary statistics
fprintf('\nSummary Statistics:\n');
fprintf('Sum Rate - Min: %.4f, Max: %.4f, Mean: %.4f\n', min(SR_all(:)), max(SR_all(:)), mean(SR_all(:)));
fprintf('CRB Trace - Min: %.4e, Max: %.4e, Mean: %.4e\n', min(CRB_all(:)), max(CRB_all(:)), mean(CRB_all(:)));

% Plot convergence for the last scenario
if ~isempty(Con_last)
    figure(1)
    plot(Con_last(:,1), 'b-', 'LineWidth', 1.5)
    grid on
    xlabel('Iteration')
    ylabel('Sum Rate (nat/s/Hz)')
    title('Sum Rate Convergence (Last Scenario)')

    figure(2)
    plot(Con_last(:,2), 'r-', 'LineWidth', 1.5)
    grid on
    xlabel('Iteration')
    ylabel('Trace of inverse FIM')
    title('Sensing Performance Convergence (Last Scenario)')

    figure(3)
    plot(Con_last(:,3), 'g-', 'LineWidth', 1.5)
    grid on
    xlabel('Iteration')
    ylabel('Total Objective Value')
    title('Objective Function Convergence (Last Scenario)')
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

