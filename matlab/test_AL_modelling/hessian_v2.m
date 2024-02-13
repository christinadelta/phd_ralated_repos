function H = hessian_v2(f, x)

% created: @christinadelta 22/01/2024
% modified: February 2024

% f is the function handle
% x is the point at which to evaluate the Hessian (e.g., x_opt)

% Update: updated hessian matrix using central difference.

% Compute the function value at the original point x
f_0 = f(x);
n   = numel(x);   % how many parameters?
H   = zeros(n);   % initialize the Hessian matrix
h   = 1e-4;       % Step size for finite difference

%% compute the matrix 

% iterate over all combinations of the elements in x, 
% to compute each element of the Hessian matrix
for i = 1:n
    for j = 1:n
        x_ij_plus           = x; 
        x_ij_minus          = x;  % Initialize x copies for plus and minus adjustments
        x_ij_plus([i, j])   = x_ij_plus([i, j]) + h;  % Increment both i and j by h
        x_ij_minus([i, j])  = x_ij_minus([i, j]) - h; % Decrement both i and j by h

        x_i_plus            = x; 
        x_i_plus(i)         = x_i_plus(i) + h; % Increment i by h
        x_i_minus           = x; 
        x_i_minus(i)        = x_i_minus(i) - h; % Decrement i by h

        x_j_plus            = x; 
        x_j_plus(j)         = x_j_plus(j) + h; % Increment j by h
        x_j_minus           = x; 
        x_j_minus(j)        = x_j_minus(j) - h; % Decrement j by h

        % Compute function values for central difference approximation
        f_ij_plus   = f(x_ij_plus);
        f_ij_minus  = f(x_ij_minus);

        f_i_plus    = f(x_i_plus);
        f_i_minus   = f(x_i_minus);

        f_j_plus    = f(x_j_plus);
        f_j_minus   = f(x_j_minus);

        % Update the Hessian matrix using central difference
        H(i,j)      = (f_ij_plus - f_i_plus - f_j_plus + 2*f_0 - f_ij_minus + f_i_minus + f_j_minus) / (2*h^2);
    end
end

end % end of function