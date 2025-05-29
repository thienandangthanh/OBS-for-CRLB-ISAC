% test_steer_matrix_generate_data.m
% Generate test data for steering matrix function verification

clc; clear;

% Get the directory where this script is located
script_path = mfilename('fullpath');
[script_dir, ~, ~] = fileparts(script_path);

% Test cases with different parameters
test_cases = struct();

% Test Case 1: Small array
test_cases.case1.Mx = 3;
test_cases.case1.My = 3;
test_cases.case1.theta = [-pi/6, pi/6];
test_cases.case1.phi = [-pi/3, pi/4];

% Test Case 2: Medium array (same as main script)
test_cases.case2.Mx = 4;
test_cases.case2.My = 4;
test_cases.case2.theta = [-pi/3 + 2*pi/3*rand(2,1)];
test_cases.case2.phi = [-pi/3 + 2*pi/3*rand(2,1)];

% Test Case 3: Larger array
test_cases.case3.Mx = 5;
test_cases.case3.My = 4;
test_cases.case3.theta = [-pi/4, 0, pi/4];
test_cases.case3.phi = [-pi/2, 0, pi/2];

% Generate results for each test case
test_data = struct();

for case_name = fieldnames(test_cases)'
    case_name = case_name{1};
    params = test_cases.(case_name);
    
    fprintf('Generating data for %s...\n', case_name);
    
    % Call the MATLAB function
    [A, dAtheta, dAphi] = construct_steer_matrix_and_derivative_steer_matrix(...
        params.theta, params.phi, params.Mx, params.My);
    
    % Store inputs and outputs
    test_data.(case_name).inputs.theta = params.theta;
    test_data.(case_name).inputs.phi = params.phi;
    test_data.(case_name).inputs.Mx = params.Mx;
    test_data.(case_name).inputs.My = params.My;
    
    test_data.(case_name).outputs.A = A;
    test_data.(case_name).outputs.dAtheta = dAtheta;
    test_data.(case_name).outputs.dAphi = dAphi;
    
    % Store some metadata
    test_data.(case_name).metadata.array_size = [params.Mx, params.My];
    test_data.(case_name).metadata.num_angles = length(params.theta);
    test_data.(case_name).metadata.output_size = size(A);
end

% Save test data
output_file = fullfile(script_dir, 'test_steer_matrix_data.mat');
save(output_file, 'test_data');
fprintf('Test data saved to test_steer_matrix_data.mat\n');
