function trZ00 = limSGF( ...
A1, B1, C1, A0, B0, C0, ...
tol, maxit, m, use_symamd)
% [ A0 B0            ]
% [ C0 A1 B1         ]
% [    C1 A1 B1      ]
% [       C1 A1 B1   ]
% [       ... ... ...]
%
% Inputs:
% A1 : n x n bulk intra-layer matrix
% B1 : n x n forward bulk coupling matrix
% C1 : n x n backward bulk coupling matrix
% A0 : n0 x n0 boundary-layer matrix
% B0 : n0 x n boundary-to-bulk coupling matrix
% C0 : n x n0 bulk-to-boundary coupling matrix
% tol : Riccati tolerance, default 1e-8
% maxit : maximum number of iterations, default 8000
% m : number of Rademacher vectors, default 200
% use_symamd : use symamd ordering, default true
%
% Output:
% trZ00 : Hutchinson estimate of trace(Z00)

    if nargin < 7, tol = 1e-8;end
    if nargin < 8, maxit = 8000;end
    if nargin < 9, m = 200;end
    if nargin < 10, use_symamd = true;end
    
    % ---------------------------------------------------------
    % 1. Low-rank decomposition of B1,C1,B0 and C0
    % ---------------------------------------------------------
    %B1 = U_B1 * V_B1'
    [Q_B1, R_B1, ~] = qr(B1, 0);
    diagR_B1 = abs(diag(R_B1));
    rB1 = sum(diagR_B1 > 1e-12 * max(diagR_B1));
    U_B1 = Q_B1(:, 1:rB1);
    V_B1 = B1' * U_B1;
    
    %C1 = U_C1 * V_C1'
    [Q_C1, R_C1, ~] = qr(C1, 0);
    diagR_C1 = abs(diag(R_C1));
    rC1 = sum(diagR_C1 > 1e-12 * max(diagR_C1));
    U_C1 = Q_C1(:, 1:rC1);
    V_C1 = C1' * U_C1;
    
    %B0 = U_B0 * V_B0'
    [Q_B0, R_B0, ~] = qr(B0, 0);
    diagR_B0 = abs(diag(R_B0));
    rB0 = sum(diagR_B0 > 1e-12 * max(diagR_B0));
    U_B0 = Q_B0(:, 1:rB0);
    V_B0 = B0' * U_B0;
    
    %C0 = U_C0 * V_C0'
    [Q_C0, R_C0, ~] = qr(C0, 0);
    diagR_C0 = abs(diag(R_C0));
    rC0 = sum(diagR_C0 > 1e-12 * max(diagR_C0));
    U_C0 = Q_C0(:, 1:rC0);
    V_C0 = C0' * U_C0;
    
    % ---------------------------------------------------------
    % 2. Reduced Riccati quantities
    %
    % Define
    %   Lt = V_B1' * L
    % and
    % P = V_B1' * A1^{-1} * U_C1; Q = V_B1' * A1^{-1} * U_B1
    % R = V_C1' * A1^{-1} * U_B1; T = V_C1' * A1^{-1} * U_C1
    % ---------------------------------------------------------
    
    % Sparse LU factorization of A0 and A1
    n0 = size(A0, 1);
    if use_symamd
        pA0_1 = symamd(A0);
        pA0_1 = pA0_1(:);
    else
        pA0_1 = (1:n0).';
    end
    
    A0_p = A0(pA0_1, pA0_1);
    [LA0, UA0, pA0_2, qA0_2] = lu(A0_p, 'vector');
    
    % A1
    n1 = size(A1, 1);
    if use_symamd
        pA1_1 = symamd(A1);
        pA1_1 = pA1_1(:);
    else
        pA1_1 = (1:n1).';
    end
    
    A1_p = A1(pA1_1, pA1_1);
    [LA1, UA1, pA1_2, qA1_2] = lu(A1_p, 'vector');
    
    A1_U_B1 = solve_A(U_B1, LA1, UA1, pA1_2, qA1_2, pA1_1);
    A1_U_C1 = solve_A(U_C1, LA1, UA1, pA1_2, qA1_2, pA1_1);
    
    P = V_B1' * A1_U_C1;Q = V_B1' * A1_U_B1;
    R = V_C1' * A1_U_B1;T = V_C1' * A1_U_C1;
    
    % ---------------------------------------------------------
    % 3. Reduced Riccati iteration
    %    Lt_{k+1} = -P + Q * (I + Lt_k * R)^{-1} * Lt_k * T
    % ---------------------------------------------------------
    Lt = zeros(rB1, rC1, 'like', A1);
    I_rB1 = eye(rB1, 'like', A1);
    
    for iter = 1:maxit
        K = I_rB1 + Lt * R;
        M = Lt * T;
        Y = K \ M;
        Lt_new = -P + Q * Y;
    
        diff = norm(Lt_new - Lt, 'fro') / (norm(Lt, 'fro') + eps);
        Lt = Lt_new;
           if diff < tol
                break;
            end
    end
    
    if iter == maxit
        warning('Riccati iteration did not converge within maxit steps');
    end
    
    % ---------------------------------------------------------
    % 4. Bulk surface Green function
    %   z11 = inv(A1 + B1*R)
    %   z11*z = A1^{-1}z
    %       - A1^{-1}U_B1
    %       * (I + Lt*R)^{-1}
    %       * Lt*V_C1'*A1^{-1}z
    % ---------------------------------------------------------
    Kbulk = I_rB1 + Lt * R;
    A1_U_C0 = solve_A(U_C0, LA1, UA1, pA1_2, qA1_2, pA1_1);
    
    rhs_bulk = Lt * (V_C1' * A1_U_C0);
    S_bulk = Kbulk \ rhs_bulk;
    
    z11_U_C0 = A1_U_C0 - A1_U_B1 * S_bulk;
    % ---------------------------------------------------------
    % 5. Boundary effective matrix
    %
    %   H = V_B0' * z11 * U_C0
    %
    %   A0 - B0*z11*C0
    %       = A0 - U_B0*H*V_C0'
    %
    %   Z00 = A0^{-1}
    %       + A0^{-1}U_B0
    %         *(I - H*V_C0'*A0^{-1}U_B0)^{-1}
    %         *H*V_C0'*A0^{-1}.
    % ---------------------------------------------------------
    Lt00 = V_B0' * z11_U_C0;
    A0_U_B0 = solve_A(U_B0, LA0, UA0, pA0_2, qA0_2, pA0_1);
    V_C0_A0_U_B0 = V_C0' * A0_U_B0;
    
    Kbd = eye(rB0, 'like', A0) - Lt00 * V_C0_A0_U_B0; 
    % ---------------------------------------------------------
    % 6. Batched Hutchinson trace estimation
    %
    %   Z00*z =
    %       A0^{-1}z
    %       + A0^{-1}U_B0
    %         Kbd^{-1}
    %         H*V_C0'*A0^{-1}z.
    % ---------------------------------------------------------
    
    z = 2 * (rand(n0, m) >= 0.5);
    W = solve_A(z, LA0, UA0, pA0_2, qA0_2, pA0_1);
    
    rhs_boundary = Lt00 * (V_C0' * W);
    S_boundary = Kbd \ rhs_boundary;
    Z00z = W + A0_U_B0 * S_boundary;
    quad_forms = sum(conj(z) .* Z00z, 1);
    trZ00 = sum(quad_forms) / m;
end


function x = solve_A(b, L, U, p2, q2, p1)
    b1 = b(p1, :);
    b2 = b1(p2, :);
    y = L \ b2;
    z = U \ y;
    y_original = zeros(size(b1), 'like', b1);
    y_original(q2, :) = z;
    x = zeros(size(b), 'like', b);
    x(p1, :) = y_original;
end