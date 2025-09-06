function FIM=calculateFIM(L,noise_s,W_or_Rx,A,dAtheta,dAphi,B,dBtheta,dBphi,U)
% calculateFIM - Calculate Fisher Information Matrix for ISAC system
%
% This function computes the Fisher Information Matrix (FIM) for 
% Integrated Sensing and Communication (ISAC) systems, used to derive
% the Cramer-Rao Lower Bound (CRLB) for parameter estimation.
%
% Inputs:
%   L        - Number of symbols/snapshots
%   noise_s  - Sensing noise power  
%   W_or_Rx  - Either beamforming matrix W (Nt x total_streams) or receive covariance matrix Rx (Nt x Nt)
%   A        - Transmit steering matrix (Nt x M)
%   dAtheta  - Partial derivative of A with respect to elevation angle theta (Nt x M)
%   dAphi    - Partial derivative of A with respect to azimuth angle phi (Nt x M)
%   B        - Receive steering matrix (Nr x M) 
%   dBtheta  - Partial derivative of B with respect to elevation angle theta (Nr x M)
%   dBphi    - Partial derivative of B with respect to azimuth angle phi (Nr x M)
%   U        - Diagonal matrix of reflection coefficients (M x M)
%
% Output:
%   FIM      - Fisher Information Matrix (4M x 4M)
%              Block structure: [theta_params; phi_params; real_alpha; imag_alpha]
%
% The FIM is constructed for parameter vector:
% [theta_1, ..., theta_M, phi_1, ..., phi_M, Re(alpha_1), ..., Re(alpha_M), Im(alpha_1), ..., Im(alpha_M)]
%
% Reference: 
% - ISAC optimization and CRLB analysis for target parameter estimation

% Determine if input is W (beamforming matrix) or Rx (covariance matrix)
[Nt, num_cols] = size(W_or_Rx);
if Nt == num_cols
    % Square matrix - assume it's Rx (covariance matrix)
    Rx = W_or_Rx;
else
    % Non-square matrix - assume it's W (beamforming matrix)
    Rx = W_or_Rx * W_or_Rx';
end

% Get number of targets
M=size(U,1);

% Initialize FIM based on input type (for CVX compatibility)
if isa(W_or_Rx,'double')
    FIM=zeros(4*M);
else
    FIM=cvx(zeros(4*M));
end

% Compute intermediate Fisher information terms
% F11: theta-theta block
F11 = (U * A' * Rx * A * U').' .* (dBtheta' * dBtheta) + ...
      (U * A' * Rx * dAtheta * U').' .* (B' * dBtheta) + ...
      (U * dAtheta' * Rx * A * U').' .* (dBtheta' * B) + ...
      (U * dAtheta' * Rx * dAtheta * U').' .* (B' * B);

% F12: theta-phi block  
F12 = (U * A' * Rx * A * U').' .* (dBtheta' * dBphi) + ...
    (U * A' * Rx * dAtheta * U').' .* (B' * dBphi) + ...
      (U * dAphi' * Rx * A * U').' .* (dBtheta' * B) + ...
      (U * dAphi' * Rx * dAtheta * U').' .* (B' * B);

% F13: theta-alpha block (complex coefficients)
F13 = (A' * Rx * A * U').' .* (dBtheta' * B) + ...
      (A' * Rx * dAtheta * U').' .* (B' * B);

% F22: phi-phi block
F22 = (U * A' * Rx * A * U').' .* (dBphi' * dBphi) + ...
      (U * A' * Rx * dAphi * U').' .* (B' * dBphi) + ...
      (U * dAphi' * Rx * A * U').' .* (dBphi' * B) + ...
      (U * dAphi' * Rx * dAphi * U').' .* (B' * B);

% F23: phi-alpha block (complex coefficients)  
F23 = (A' * Rx * A * U').' .* (dBphi' * B) + ...
      (A' * Rx * dAphi * U').' .* (B' * B);

% F33: alpha-alpha block (complex coefficients)
F33 = (A' * Rx * A).' .* (B' * B);

% Fill the FIM blocks
% Real parts for angle parameters and real parts of complex coefficients
FIM(1:M,1:M)=real(F11);                    % theta-theta  
FIM(1:M,M+1:2*M)=real(F12);                % theta-phi
FIM(1:M,2*M+1:3*M)=real(F13);              % theta-Re(alpha)
FIM(1:M,3*M+1:4*M)=-imag(F13);             % theta-Im(alpha)
FIM(M+1:2*M,M+1:2*M)=real(F22);            % phi-phi
FIM(M+1:2*M,2*M+1:3*M)=real(F23);          % phi-Re(alpha)  
FIM(M+1:2*M,3*M+1:4*M)=-imag(F23);         % phi-Im(alpha)
FIM(2*M+1:3*M,2*M+1:3*M)=real(F33);        % Re(alpha)-Re(alpha)
FIM(2*M+1:3*M,3*M+1:4*M)=-imag(F33);       % Re(alpha)-Im(alpha)
FIM(3*M+1:4*M,3*M+1:4*M)=real(F33);        % Im(alpha)-Im(alpha)

% Make the matrix Hermitian by symmetry
FIM=triu(FIM)+triu(FIM)'-FIM.*eye(4*M);

% Apply scaling factor
FIM=2*L/noise_s*FIM;

end 