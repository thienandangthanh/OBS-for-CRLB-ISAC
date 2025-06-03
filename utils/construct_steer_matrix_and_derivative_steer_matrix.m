function [A,dAtheta,dAphi]=construct_steer_matrix_and_derivative_steer_matrix(theta,phi,Mx,My)
    % Construct steering matrix and its derivatives for UPA (Uniform Planar Array)
    %
    % Inputs:
    %   theta - Elevation angle (in radians)
    %   phi   - Azimuth angle (in radians)
    %   Mx    - Number of elements along x-axis
    %   My    - Number of elements along y-axis
    %
    % Outputs:
    %   A       - Steering matrix
    %   dAtheta - Derivative of steering matrix with respect to theta
    %   dAphi   - Derivative of steering matrix with respect to phi

    % Generate indices for x-axis and y-axis
    ix = (0:Mx-1)'; % Column vector for x-axis indices
    iy = (0:My-1)'; % Column vector for y-axis indices

    % Number of steering vectors
    M = length(theta);

    % Preallocate outputs
    A = zeros(Mx*My,M);
    dAtheta = zeros(Mx*My,M);
    dAphi = zeros(Mx*My,M);

    % Compute steering vectors along x and y axes
    for m=1:M
        % Steering vectors for x and y dimensions
        ax = 1/sqrt(Mx)* exp(1j * pi * ix * sin(theta(m)) * sin(phi(m)));
        ay = 1/sqrt(My)* exp(1j * pi * iy * cos(phi(m)));

        % Derivatives of ax and ay
        daxtheta = 1j*pi*ix*cos(theta(m))*sin(phi(m)).*ax;
        daxphi = 1j*pi*ix*sin(theta(m))*cos(phi(m)).*ax;
        dayphi = -1j*pi*iy*sin(phi(m)).*ay;

        % Compute the full UPA steering vector using Kronecker product
        A(:,m) = kron(ay, ax); % Kronecker product combines ax and ay
        dAtheta(:,m) = kron(ay,daxtheta);
        dAphi(:,m) = kron(ay,daxphi)+kron(dayphi,ax);
    end
end
