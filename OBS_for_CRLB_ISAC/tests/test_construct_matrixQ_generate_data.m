function test_construct_matrixQ_matlab()
% TEST_CONSTRUCT_MATRIXQ_MATLAB Generate test data for construct_matrixQ function
%
% This function creates test inputs and computes the output using MATLAB's
% construct_matrixQ function. The results are saved to a .mat file that
% can be used by the Python test to verify equivalence.

fprintf('Generating test data for construct_matrixQ function...\n');

% Get the directory where this script is located
script_path = mfilename('fullpath');
[script_dir, ~, ~] = fileparts(script_path);

% Set random seed for reproducibility
rng(42);

% Test case 1: Small matrices
M1 = 2;
Nt1 = 8;
Nr1 = 10;

% Generate test inputs
L1 = 30;
noise_s1 = 1e-3;

% Create random matrices with proper dimensions
Phi1 = randn(4*M1, 4*M1) + 1j*randn(4*M1, 4*M1);
Phi1 = Phi1 + Phi1'; % Make Hermitian

A1 = randn(Nt1, M1) + 1j*randn(Nt1, M1);
dAtheta1 = randn(Nt1, M1) + 1j*randn(Nt1, M1);
dAphi1 = randn(Nt1, M1) + 1j*randn(Nt1, M1);

B1 = randn(Nr1, M1) + 1j*randn(Nr1, M1);
dBtheta1 = randn(Nr1, M1) + 1j*randn(Nr1, M1);
dBphi1 = randn(Nr1, M1) + 1j*randn(Nr1, M1);

alpha1 = randn(M1, 1) + 1j*randn(M1, 1);
U1 = diag(alpha1);

% Compute output
Q1 = construct_matrixQ(L1, noise_s1, Phi1, A1, dAtheta1, dAphi1, B1, dBtheta1, dBphi1, U1);

% Test case 2: Larger matrices with different parameters
M2 = 3;
Nt2 = 12;
Nr2 = 15;

L2 = 50;
noise_s2 = 2e-3;

Phi2 = randn(4*M2, 4*M2) + 1j*randn(4*M2, 4*M2);
Phi2 = Phi2 + Phi2'; % Make Hermitian

A2 = randn(Nt2, M2) + 1j*randn(Nt2, M2);
dAtheta2 = randn(Nt2, M2) + 1j*randn(Nt2, M2);
dAphi2 = randn(Nt2, M2) + 1j*randn(Nt2, M2);

B2 = randn(Nr2, M2) + 1j*randn(Nr2, M2);
dBtheta2 = randn(Nr2, M2) + 1j*randn(Nr2, M2);
dBphi2 = randn(Nr2, M2) + 1j*randn(Nr2, M2);

alpha2 = randn(M2, 1) + 1j*randn(M2, 1);
U2 = diag(alpha2);

Q2 = construct_matrixQ(L2, noise_s2, Phi2, A2, dAtheta2, dAphi2, B2, dBtheta2, dBphi2, U2);

% Test case 3: Edge case with very small values
M3 = 1;
Nt3 = 4;
Nr3 = 5;

L3 = 10;
noise_s3 = 1e-6;

Phi3 = 1e-8 * (randn(4*M3, 4*M3) + 1j*randn(4*M3, 4*M3));
Phi3 = Phi3 + Phi3'; % Make Hermitian

A3 = 1e-3 * (randn(Nt3, M3) + 1j*randn(Nt3, M3));
dAtheta3 = 1e-3 * (randn(Nt3, M3) + 1j*randn(Nt3, M3));
dAphi3 = 1e-3 * (randn(Nt3, M3) + 1j*randn(Nt3, M3));

B3 = 1e-3 * (randn(Nr3, M3) + 1j*randn(Nr3, M3));
dBtheta3 = 1e-3 * (randn(Nr3, M3) + 1j*randn(Nr3, M3));
dBphi3 = 1e-3 * (randn(Nr3, M3) + 1j*randn(Nr3, M3));

alpha3 = 1e-2 * (randn(M3, 1) + 1j*randn(M3, 1));
U3 = diag(alpha3);

Q3 = construct_matrixQ(L3, noise_s3, Phi3, A3, dAtheta3, dAphi3, B3, dBtheta3, dBphi3, U3);

% Save all test data
test_data = struct();

% Test case 1
test_data.test1.inputs.L = L1;
test_data.test1.inputs.noise_s = noise_s1;
test_data.test1.inputs.Phi = Phi1;
test_data.test1.inputs.A = A1;
test_data.test1.inputs.dAtheta = dAtheta1;
test_data.test1.inputs.dAphi = dAphi1;
test_data.test1.inputs.B = B1;
test_data.test1.inputs.dBtheta = dBtheta1;
test_data.test1.inputs.dBphi = dBphi1;
test_data.test1.inputs.U = U1;
test_data.test1.output = Q1;

% Test case 2
test_data.test2.inputs.L = L2;
test_data.test2.inputs.noise_s = noise_s2;
test_data.test2.inputs.Phi = Phi2;
test_data.test2.inputs.A = A2;
test_data.test2.inputs.dAtheta = dAtheta2;
test_data.test2.inputs.dAphi = dAphi2;
test_data.test2.inputs.B = B2;
test_data.test2.inputs.dBtheta = dBtheta2;
test_data.test2.inputs.dBphi = dBphi2;
test_data.test2.inputs.U = U2;
test_data.test2.output = Q2;

% Test case 3
test_data.test3.inputs.L = L3;
test_data.test3.inputs.noise_s = noise_s3;
test_data.test3.inputs.Phi = Phi3;
test_data.test3.inputs.A = A3;
test_data.test3.inputs.dAtheta = dAtheta3;
test_data.test3.inputs.dAphi = dAphi3;
test_data.test3.inputs.B = B3;
test_data.test3.inputs.dBtheta = dBtheta3;
test_data.test3.inputs.dBphi = dBphi3;
test_data.test3.inputs.U = U3;
test_data.test3.output = Q3;

% Save to file
output_file = fullfile(script_dir, 'test_construct_matrixQ_data.mat');
save(output_file, 'test_data');

fprintf('Test data saved to test_construct_matrixQ_data.mat\n');
fprintf('Generated %d test cases\n', 3);

% Display summary
fprintf('\nTest case summary:\n');
fprintf('Test 1: M=%d, Nt=%d, Nr=%d, L=%.1f, noise_s=%.1e\n', M1, Nt1, Nr1, L1, noise_s1);
fprintf('Test 2: M=%d, Nt=%d, Nr=%d, L=%.1f, noise_s=%.1e\n', M2, Nt2, Nr2, L2, noise_s2);
fprintf('Test 3: M=%d, Nt=%d, Nr=%d, L=%.1f, noise_s=%.1e\n', M3, Nt3, Nr3, L3, noise_s3);

end
