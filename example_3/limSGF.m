function trZint = limSGF( ...
    A11, B11, C11, ...
    A22, B22, C22, ...
    B12, C21, ...
    tol, maxit, m)

% [C11 A11 B11              ]
% [    C11 A11 B12          ]
% [        C21 A22 B22      ]
% [            C22 A22 B22  ]
% [             ... ... ... ]   
% Inputs:
%
% A11, B11, C11 : left bulk blocks, nLeft-by-nLeft
% A22, B22, C22 : right bulk blocks, nRight-by-nRight
% B12 : left surface row, right surface column, size nLeft-by-nRight
% C21 : right surface row, left surface column, size nRight-by-nLeft
% tol   : reduced Riccati tolerance, default 1e-6
% maxit : maximum Riccati iterations, default 8000
% m     : number of Hutchinson vectors, default 200
%
% Output:
%
% trZint : Hutchinson estimate of
%
%          trace(G11) + trace(G22)
%
%          where [G11 G12; G21 G22] is the Green function
%          restricted to the two joined surface blocks.


    % Default parameters
    if nargin < 9 ,tol = 1e-6; end
    if nargin < 10 ,maxit = 8000; end
    if nargin < 11 ,m = 200;end
    
    % ---------------------------------------------------------
    % 1. QR decompositions of B11, C11, B22, C22, B12, C21
    % X = U_X*V_X'
    % ---------------------------------------------------------
    rankTol = 1e-12;
    blocks = {B11, C11, B22, C22, B12, C21};
    U_blocks = cell(1, 6);V_blocks = cell(1, 6);
    ranks = zeros(1, 6);
    
    for j = 1:6
        block = blocks{j};
        [Q, R, ~] = qr(block, 0);
        diagR = abs(diag(R));
        r = sum(diagR > rankTol * max(diagR));
        if r == 0
            error([num2str(block),' is numerically zero']);
        end
        U_blocks{j} = Q(:, 1:r);
        V_blocks{j} = block' * U_blocks{j};
        ranks(j) = r;
    end
    
    U_B11 = U_blocks{1};V_B11 = V_blocks{1};
    U_C11 = U_blocks{2};V_C11 = V_blocks{2};
    U_B22 = U_blocks{3};V_B22 = V_blocks{3};
    U_C22 = U_blocks{4};V_C22 = V_blocks{4};
    U_B12 = U_blocks{5};V_B12 = V_blocks{5};
    U_C21 = U_blocks{6};V_C21 = V_blocks{6};
    
    r_B11 = ranks(1);r_C11 = ranks(2);
    r_B22 = ranks(3);r_C22 = ranks(4);
    r_B12 = ranks(5);r_C21 = ranks(6);
    
    clear blocks U_blocks V_blocks Q R block
    
    % ---------------------------------------------------------
    % 2. LU factorizations and reduced Riccati quantities
    %
    % Left:
    % Lt_left has size r_C11-by-r_B11.
    % gLeft = inv(A11 + U_C11*Lt_left*V_B11')
    %
    % Right:
    % Lt_right has size r_B22-by-r_C22.
    % gRight = inv(A22 + U_B22*Lt_right*V_C22')
    % ---------------------------------------------------------
    nLeft  = size(A11, 1);nRight = size(A22, 1);
    nInt   = nLeft + nRight;
    
    [L11, U11, p11, q11] = lu(sparse(A11), 'vector');
    [L22, U22, p22, q22] = lu(sparse(A22), 'vector');
    
    A11_U_B11 = solve_A(U_B11, L11, U11, p11, q11);
    A11_U_C11 = solve_A(U_C11, L11, U11, p11, q11);
    
    A22_U_B22 = solve_A(U_B22, L22, U22, p22, q22);
    A22_U_C22 = solve_A(U_C22, L22, U22, p22, q22);
    
    % Left orientation: C11*gLeft*B11
    
    P_left = V_C11' * A11_U_B11;
    Q_left = V_C11' * A11_U_C11;
    R_left = V_B11' * A11_U_C11;
    T_left = V_B11' * A11_U_B11;
    
    % Right orientation: B22*gRight*C22
    
    P_right = V_B22' * A22_U_C22;
    Q_right = V_B22' * A22_U_B22;
    R_right = V_C22' * A22_U_B22;
    T_right = V_C22' * A22_U_C22;

    % ---------------------------------------------------------
    % 3. Reduced Riccati iterations
    % Lt_new = -P + Q*(I + Lt*R)^{-1}*Lt*T
    % ---------------------------------------------------------
    Lt_left = zeros(r_C11, r_B11, 'like', A11);
    I_left  = eye(r_C11, 'like', A11);

        for iter_left = 1:maxit
            K_left = I_left + Lt_left * R_left;
            Y_left = K_left \ (Lt_left * T_left);
            Lt_left_new = -P_left + Q_left * Y_left;
    
            diff_left = norm(Lt_left_new - Lt_left, 'fro') / (norm(Lt_left_new, 'fro')+eps);
            Lt_left = Lt_left_new;
            if diff_left < tol
                break;
            end
        end
    
    if diff_left >= tol
        warning('Riccati iteration on left did not converge within maxit steps');
    end

    Lt_right = zeros(r_B22, r_C22, 'like', A22);
    I_right  = eye(r_B22, 'like', A22);

        for iter_right = 1:maxit
            K_right = I_right + Lt_right * R_right;
            Y_right = K_right \ (Lt_right * T_right);
            Lt_right_new = -P_right + Q_right * Y_right;
    
            diff_right = norm(Lt_right_new - Lt_right, 'fro') / (norm(Lt_right_new, 'fro')+ eps);
            Lt_right = Lt_right_new;
            if diff_right < tol
                break;
            end
    
        end
    
        if diff_right >= tol
            warning('Riccati iteration on left did not converge within maxit steps');
        end
    
    % ---------------------------------------------------------
    % 4. Interface Woodbury preparation
    %
    % D = blkdiag(A11, A22)
    %
    % Zint = [A11 + U_C11*Lt_left*V_B11',           B12; 
    %                  C21,               A22 + U_B22*Lt_right*V_C22']
    %
    % Write:Zint = D - U_bd*H_bd*V_bd'
    % ---------------------------------------------------------
    
    U_bd = [U_C11, zeros(nLeft, r_B22, 'like', A11), U_B12, zeros(nLeft, r_C21, 'like', A11); ...
           zeros(nRight, r_C11, 'like', A22), U_B22, zeros(nRight, r_B12, 'like', A22), U_C21];
    
    V_bd = [V_B11, zeros(nLeft, r_C22, 'like', A11), zeros(nLeft, r_B12, 'like', A11), V_C21; ...
            zeros(nRight, r_B11, 'like', A22), V_C22, V_B12, zeros(nRight, r_C21, 'like', A22)];
    
    lT_bd = blkdiag(-Lt_left, -Lt_right, -eye(r_B12, 'like', A11), -eye(r_C21, 'like', A22));
    
    A_int = [ solve_A(U_bd(1:nLeft, :), L11, U11, p11, q11); ...
               solve_A(U_bd(nLeft+1:end, :), L22, U22, p22, q22)];
    
    K_bd = eye(size(lT_bd, 1), 'like', lT_bd) ...
        - lT_bd * (V_bd' * A_int);
    
    % ---------------------------------------------------------
    % 5. Batched Hutchinson trace estimation
    %
    % Zint = inv(D - U_bd*H_bd*V_bd')
    %
    % Zint*z = D^{-1}*z
    %        + D^{-1}*U_bd*K_bd^{-1}
    %          *H_bd*V_bd'*D^{-1}*z
    % ---------------------------------------------------------
    
    z = 2 * (rand(nInt, m) >= 0.5) - 1;
    z = cast(z, 'like', full(A11));
    
    W = [solve_A(z(1:nLeft, :), ...
                L11, U11, p11, q11); ...
        solve_A(z(nLeft+1:end, :), ...
                L22, U22, p22, q22)];
    
    reduced_rhs = lT_bd * (V_bd' * W);
    GintZ = W + A_int * (K_bd \ reduced_rhs);
    
    quad_forms = sum(conj(z) .* GintZ, 1);
    trZint = mean(quad_forms);
end


function x = solve_A(b, L, U, p, q)
    b_permuted = b(p, :);
    y = L \ b_permuted;
    z = U \ y;
    x = zeros(size(z), 'like', z);
    x(q, :) = z;
end