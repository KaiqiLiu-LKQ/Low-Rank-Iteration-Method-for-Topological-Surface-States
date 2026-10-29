function trG0 = limSGF(A, B, C, tol, maxit, m, use_symamd)
% dsr_SGF  Compute trace of (1,1) block of inverse of
%    semi-infinite block Toeplitz tridiagonal matrix with blocks A, B, C.
%
%   本版本使用 symamd 排序 + LU 分解，并显式传播所有置换向量，
%   确保 solve_A 求解的是原始 A x = b，而非被置换过的系统。
%
%   Inputs:
%       A     - n×n sparse matrix (assumed invertible, near Hermitian)
%       B, C  - n×n sparse matrices, low-rank
%       tol   - tolerance for Riccati iteration (default 1e-8)
%       maxit - maximum number of iterations (default 500)
%       m     - number of random vectors for trace estimation (default 50)
%       use_symamd - (optional) true/false, 是否使用 symamd 排序 (default true)
%
%   Output:
%       trG00   - trace of (1,1) block of the inverse (complex)

    if nargin < 4, tol = 1e-8; end
    if nargin < 5, maxit = 500; end
    if nargin < 6, m = 50; end
    if nargin < 7, use_symamd = true; end

    n = size(A, 1);

    % ---------------------------------------------------------
    % 1. QR decomposition of B and C
    % ---------------------------------------------------------
    % B = U_B * V_B'
    [Q_B, R_B, ~] = qr(B, 0);
    diagR_B = abs(diag(R_B));
    rB = sum(diagR_B > 1e-12 * max(diagR_B));
    if rB == 0
        error('B is numerically zero');
    end
    U_B = Q_B(:, 1:rB);
    V_B = B' * U_B;

    % C = U_C * V_C'
    [Q_C, R_C, ~] = qr(C, 0);
    diagR_C = abs(diag(R_C));
    rC = sum(diagR_C > 1e-12 * max(diagR_C));
    if rC == 0
        error('C is numerically zero');
    end
    U_C = Q_C(:, 1:rC);
    V_C = C' * U_C;

    % ---------------------------------------------------------
    % 2. Reduced Riccati quantities
    %
    % Define
    %   Lt = V_B' * L
    % and
    %
    % P = V_B' * A^{-1} * U_C; Q = V_B' * A^{-1} * U_B
    % R = V_C' * A^{-1} * U_B; T = V_C' * A^{-1} * U_C
    % ---------------------------------------------------------

    if use_symamd
        p1 = symamd(A);                     % symamd 置换：A(p1,p1) 填充更少
    else
        p1 = (1:n)';                        % 恒等置换
    end

    [L_A, U_A, p2, q2] = lu(A(p1, p1), 'vector');
    solve_A = @(b) solve_symamd_lu(b, L_A, U_A, p2, q2, p1);


    A_U_B = solve_A(U_B);            % n×rB
    A_U_C = solve_A(U_C);            % n×rC

    P = V_B' * A_U_C;                % rB×rC
    Q = V_B' * A_U_B;                % rB×rB
    R = V_C' * A_U_B;                % rC×rB
    T = V_C' * A_U_C;                % rC×rC

    % ---------------------------------------------------------
    % 3. Reduced Riccati iteration
    %    Lt_{k+1} = -P + Q * (I + Lt_k * R)^{-1} * Lt_k * T
    % ---------------------------------------------------------
    Lt = zeros(rB, rC);
    I_rB = eye(rB);

    for iter = 1:maxit
        K = I_rB + Lt * R;
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
    % 4. Interface Woodbury preparation and batched Hutchinson trace estimation
    %    g00 = A^{-1} - A^{-1} U_B (I + Ltilde*R)^{-1} Ltilde V_C' A^{-1}
    % ---------------------------------------------------------
    Vb_L_final = Lt;
    K_final = I_rB + Vb_L_final * R;
    [Lk, Uk] = lu(K_final);

    Z = 2 * (rand(n, m) >= 0.5) - 1;% Generate m Rademacher random vectors
    W = solve_A(Z);  % solve: A*W = Z → W = A^{-1}Z
    
    Vc_W = V_C' * W;                 % rC × m
    rhs  = Vb_L_final * Vc_W;        % rB × m
    S = Uk \ (Lk \ rhs);          % K_final * S = rhs   (rB × m)
    G0Z = W - A_U_B * S;           % X * Z   (n × m)
    quad_forms = sum(conj(Z) .* G0Z, 1); % 每一列计算 z_i' * (G0 z_i)
    trG0 = sum(quad_forms) / m; % Hutchinson 迹估计
    
end

% =====================================================================
% 局部函数：用 symamd + LU 的置换求解原始 A x = b
% =====================================================================
function x = solve_symamd_lu(b, L, U, p2, q2, p1)
% 求解 A x = b，其中：
%   A(p1,p1) 是 symamd 置换后的矩阵
%   A(p1,p1)(p2,q2) = L * U
%
% 步骤：
%   1. b1 = b(p1)                     行置换，对应 A(p1,p1) 的右端项
%   2. 解 A(p1,p1) * y = b1，用 L,U,p2,q2
%   3. x(p1) = y                      逆置换，恢复原始未知数顺序

    % 步骤1：行置换
    b1 = b(p1, :);

    % 步骤2：解 A(p1,p1) * y = b1
    % A(p1,p1)(p2,q2) = L * U  =>  L*U*y(q2) = b1(p2)
    y = zeros(size(b1));
    b1_perm = b1(p2, :);
    tmp = L \ b1_perm;
    tmp2 = U \ tmp;
    y(q2, :) = tmp2;

    % 步骤3：逆置换
    x = zeros(size(b));
    x(p1, :) = y;
end