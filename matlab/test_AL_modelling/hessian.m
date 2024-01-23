function H = hessian(f, x)

% created: @christinadelta 22/01/2024

    % f is the function handle
    % x is the point at which to evaluate the Hessian (e.g., x_opt)
    
    n = numel(x);   % how many parameters?
    H = zeros(n);   % initialize the Hessian matrix
    h = 1e-4;       % Step size for finite difference
    
    %% compute the matrix 
    
    % iterate over all combinations of the elements in x, 
    % to compute each element of the Hessian matrix
    for i = 1:n

        % loop over parameters, compute eigenvalues 
        for j = 1:n

            % calculate the function value at a point where both the ith 
            % and jth parameters are incremented by h. It's part of the 
            % finite difference approximation.
            x_ij = x;
            x_ij(i) = x_ij(i) + h;
            x_ij(j) = x_ij(j) + h;
            f_ij = f(x_ij);
            
            % the function value is calculated at a point where only the 
            % ith parameter is incremented by h
            x_i = x;
            x_i(i) = x_i(i) + h;
            f_i = f(x_i);

            % calculate the function value at a point where only the jth 
            % parameter is incremented by h.
            x_j = x;
            x_j(j) = x_j(j) + h;
            f_j = f(x_j);

            % Compute the function value at the original point x
            f_0 = f(x);

            % This is the core of the finite difference approximation. 
            % It computes the second-order partial derivative with respect 
            % to the ith and jth parameters. This formula is derived from 
            % the second-order Taylor expansion and calculates an approximation 
            % of the second derivative by using the values of the function 
            % at x, x_i, x_j, and x_ij.
            H(i,j) = (f_ij - f_i - f_j + f_0) / h^2;

        end % end of inner loop
    end % end of outer loop

end % end of function
