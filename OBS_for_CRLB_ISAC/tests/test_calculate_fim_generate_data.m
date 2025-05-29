function test_calculate_fim_generate_data()
% Generate test data using MATLAB calculateFIM implementation
%
% This script generates test parameters, computes FIM using MATLAB's
% calculateFIM function, and saves the data for Python validation.
% The output file is saved in the same directory as this script.

clc
clear

fprintf('=== MATLAB Test Data Generation ===\n');

% Get the directory where this script is located
script_path = mfilename('fullpath');
[script_dir, ~, ~] = fileparts(script_path);

% Initialize test parameters (same as Python tests for consistency)
Nth = 4;
Ntv = 4;
Nt = Nth * Ntv;

Nrh = 5;
Nrv = 4;
Nr = Nrh * Nrv;

M = 2;  % Number of targets
K = 4;  % Number of communication streams

% Use fixed seed for reproducibility (same as Python)
rng(42, 'twister');

% Generate test angles
theta = -pi/3 + 2*pi/3 * rand(M,1);
phi = -pi/3 + 2*pi/3 * rand(M,1);

% Generate test parameters
L = 30;
noise_s = 1e-3;

% Generate complex reflection coefficients
alpha = 0.1 * (1 + 0.2 * randn(M,1)) .* exp(1j * 2 * pi * rand(M,1));

fprintf('Test parameters:\n');
fprintf('  L = %g\n', L);
fprintf('  noise_s = %g\n', noise_s);
fprintf('  M (targets) = %d\n', M);
fprintf('  Array config: %dx%d transmit, %dx%d receive\n', Nth, Ntv, Nrh, Nrv);

% Construct steering matrices and derivatives
tic;
[A, dAtheta, dAphi] = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nth, Ntv);
[B, dBtheta, dBphi] = construct_steer_matrix_and_derivative_steer_matrix(theta, phi, Nrh, Nrv);
fprintf('Steering matrix computation time: %.6f seconds\n', toc);

% Create diagonal matrix U
U = diag(alpha);

% Generate test beamforming matrix W
num_sensing_streams = Nt;
Wc = randn(Nt, K) + 1j * randn(Nt, K);
Ws = randn(Nt, num_sensing_streams) + 1j * randn(Nt, num_sensing_streams);

W = [Wc, Ws];

% Normalize power
Pt = 1.0;
W = W * sqrt(Pt / trace(W * W'));

fprintf('Generated beamforming matrix W: %dx%d\n', size(W,1), size(W,2));

% Calculate FIM using MATLAB implementation
fprintf('\nCalculating FIM using MATLAB reference implementation...\n');
tic;
FIM_matlab = calculateFIM(L, noise_s, W, A, dAtheta, dAphi, B, dBtheta, dBphi, U);
matlab_time = toc;
fprintf('MATLAB FIM computation time: %.6f seconds\n', matlab_time);

% Validate MATLAB results
fprintf('\nValidating MATLAB FIM properties:\n');
if ishermitian(FIM_matlab)
    fprintf('✓ FIM is Hermitian\n');
else
    fprintf('✗ FIM is NOT Hermitian\n');
end

eig_vals = eig(FIM_matlab);
if all(real(eig_vals) >= -1e-10)
    fprintf('✓ FIM is positive semi-definite (min eigenvalue: %.2e)\n', min(real(eig_vals)));
else
    fprintf('✗ FIM is NOT positive semi-definite (min eigenvalue: %.2e)\n', min(real(eig_vals)));
end

fprintf('FIM condition number: %.2e\n', cond(FIM_matlab));
fprintf('FIM dimensions: %dx%d\n', size(FIM_matlab,1), size(FIM_matlab,2));

% Test covariance matrix input equivalence
Rx = W * W';
FIM_matlab_Rx = calculateFIM(L, noise_s, Rx, A, dAtheta, dAphi, B, dBtheta, dBphi, U);
covar_diff = max(abs(FIM_matlab - FIM_matlab_Rx), [], 'all');
if covar_diff < 1e-12
    fprintf('✓ W and Rx inputs produce identical results (max diff: %.2e)\n', covar_diff);
else
    fprintf('✗ W and Rx inputs differ (max diff: %.2e)\n', covar_diff);
end

% Prepare reference data for Python validation
matlab_reference_data = struct();
matlab_reference_data.L = L;
matlab_reference_data.noise_s = noise_s;
matlab_reference_data.W = W;
matlab_reference_data.A = A;
matlab_reference_data.dAtheta = dAtheta;
matlab_reference_data.dAphi = dAphi;
matlab_reference_data.B = B;
matlab_reference_data.dBtheta = dBtheta;
matlab_reference_data.dBphi = dBphi;
matlab_reference_data.U = U;
matlab_reference_data.FIM_matlab_reference = FIM_matlab;
matlab_reference_data.theta = theta;
matlab_reference_data.phi = phi;
matlab_reference_data.alpha = alpha;
matlab_reference_data.Nth = Nth;
matlab_reference_data.Ntv = Ntv;
matlab_reference_data.Nrh = Nrh;
matlab_reference_data.Nrv = Nrv;
matlab_reference_data.M = M;
matlab_reference_data.K = K;
matlab_reference_data.matlab_computation_time = matlab_time;
matlab_reference_data.matlab_condition_number = cond(FIM_matlab);

% Save reference data in the same directory as this script
output_file = fullfile(script_dir, 'test_calculate_fim_data.mat');
save(output_file, 'matlab_reference_data');

fprintf('\n=== Test Data Generation Complete ===\n');
fprintf('Reference data saved to: %s\n', output_file);
fprintf('File size: %.1f KB\n', dir(output_file).bytes / 1024);
fprintf('\nUse this data in Python validation by loading:\n');
fprintf('  scipy.io.loadmat(''%s'')\n', output_file);

end 