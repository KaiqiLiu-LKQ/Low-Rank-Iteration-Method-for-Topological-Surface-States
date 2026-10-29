function trZ00 = limSGF( ...
A11, B11, C11, ...
A22, B22, C22, ...
A00, B01, C10, B02, C20, ...
tol, maxit, m)

% [C11 A11 B11                  ]
% [    C11 A11 B01              ]
% [        C10 A00 B02          ]
% [            C20 A22 B22      ]
% [                C22 A22 B22  ]
% [                ... ... ...  ]
%
% Inputs:
%
% A11, B11, C11 : left bulk blocks
% A22, B22, C22 : right bulk blocks
% A00 : interface block
% B01 : left bulk row, interface column
% C10 : interface row, left bulk column
% B02 : interface row, right bulk column
% C20 : right bulk row, interface column
% tol : Riccati tolerance, default 1e-6
% maxit : maximum Riccati iterations, default 8000
% m : number of Hutchinson vectors, default 200
%
% Output:
%
% trG00 : Hutchinson estimate of trace(G00)


% Default parameters
if nargin < 12 ,tol = 1e-6; end
if nargin < 13 ,maxit = 8000; end
if nargin < 14 ,m = 200; end

% ---------------------------------------------------------
% 1. QR decomposition of B01,C10,B02,C20,B11,C11,B22 and C22
% X = U_X*V_X'
% ---------------------------------------------------------
rankTol = 1e-12;
blocks = {B01, C10, B02, C20, B11, C11, B22, C22};
U_blocks = cell(1, 8);
V_blocks = cell(1, 8);
ranks = zeros(1, 8);

for j = 1:8
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

U_B01 = U_blocks{1};V_B01 = V_blocks{1};
U_C10 = U_blocks{2};V_C10 = V_blocks{2};
U_B02 = U_blocks{3};V_B02 = V_blocks{3};
U_C20 = U_blocks{4};V_C20 = V_blocks{4};
U_B11 = U_blocks{5};V_B11 = V_blocks{5};
U_C11 = U_blocks{6};V_C11 = V_blocks{6};
U_B22 = U_blocks{7};V_B22 = V_blocks{7};
U_C22 = U_blocks{8};V_C22 = V_blocks{8};

% r_B01 = ranks(1);r_C10 = ranks(2);
% r_B02 = ranks(3);r_C20 = ranks(4);
r_B11 = ranks(5);r_C11 = ranks(6);
r_B22 = ranks(7);r_C22 = ranks(8);

clear blocks U_blocks V_blocks Q R block

% ---------------------------------------------------------
% 2. Reduced Riccati quantities
%
% Define
%   Lt = V_B' * L
% and
% P = V_B' * A^{-1} * U_C; Q = V_B' * A^{-1} * U_B
% R = V_C' * A^{-1} * U_B; T = V_C' * A^{-1} * U_C
% ---------------------------------------------------------

n0 = size(A00, 1);
[L00, U00, p00, q00] = lu(sparse(A00), 'vector');
[L11, U11, p11, q11] = lu(sparse(A11), 'vector');
[L22, U22, p22, q22] = lu(sparse(A22), 'vector');

A11_U_C11 = solve_A(U_C11, L11, U11, p11, q11);
A11_U_B11 = solve_A(U_B11, L11, U11, p11, q11);
A11_U_C10 = solve_A(U_C10, L11, U11, p11, q11);


A22_U_B22 = solve_A(U_B22, L22, U22, p22, q22);
A22_U_C22 = solve_A(U_C22, L22, U22, p22, q22);
A22_U_C20 = solve_A(U_C20, L22, U22, p22, q22);

P_left = V_B11' * A11_U_C11;P_right = V_B22' * A22_U_C22;
Q_left = V_B11' * A11_U_B11;Q_right = V_B22' * A22_U_B22;
R_left = V_C11' * A11_U_B11;R_right = V_C22' * A22_U_B22;
T_left = V_C11' * A11_U_C11;T_right = V_C22' * A22_U_C22;

Lt_left  = zeros(r_C11, r_B11, 'like', A11);
Lt_right = zeros(r_B22, r_C22, 'like', A22);
I_right = eye(r_B22, 'like', A22);I_left = eye(r_C11, 'like', A11);

% ---------------------------------------------------------
% 3. Reduced Riccati iteration
%    Lt_{k+1} = -P + Q * (I + Lt_k * R)^{-1} * Lt_k * T
% ---------------------------------------------------------

for iter_left = 1:maxit
    K_left = I_left + Lt_left * R_left;
    M_left = Lt_left * T_left;
    Y_left = K_left \ M_left;
    Lt_left_new = -P_left + Q_left * Y_left;

    diff_left = norm(Lt_left_new - Lt_left, 'fro') / (norm(Lt_left, 'fro') + eps);
    Lt_left = Lt_left_new;
       if diff_left < tol
            break;
        end
end
if iter_left == maxit
    warning('Riccati iteration on left did not converge within maxit steps');
end

for iter_right = 1:maxit
    K_right = I_right + Lt_right * R_right;
    M_right = Lt_right * T_right;
    Y_right = K_right \ M_right;
    Lt_right_new = -P_right + Q_right * Y_right;

    diff_right = norm(Lt_right_new - Lt_right, 'fro') / (norm(Lt_right, 'fro') + eps);
    Lt_right = Lt_right_new;
       if diff_right < tol
            break;
        end
end
if iter_right == maxit
    warning('Riccati iteration on right did not converge within maxit steps');
end

% ---------------------------------------------------------
% 4. Interface Woodbury preparation and batched Hutchinson trace estimation
%    G00 = inv(A00 - U_bd*H_bd*V_bd')
%    G00*Z = A00^{-1}*Z
%           +A00^{-1}*U_bd*K_bd^{-1}
%           *H_bd*V_bd'*A00^{-1}*Z
% ---------------------------------------------------------

K_left = I_left + Lt_left * R_left;
K_right = I_right + Lt_right * R_right;

Z_left_U_C10 = A11_U_C10 - A11_U_B11 * (K_left \ (Lt_left * (V_C11' * A11_U_C10)));
ZRight_U_C20 = A22_U_C20 - A22_U_B22 * (K_right \ (Lt_right * (V_C22' * A22_U_C20)));

H_left = V_B01' * Z_left_U_C10;H_right = V_B02' * ZRight_U_C20;

% Combine left and right self-energies

U_bd = [U_B01, U_B02];
V_bd = [V_C10, V_C20];
Ltint = blkdiag(H_left, H_right);
A00_U_bd = solve_A(U_bd, L00, U00, p00, q00);
K_bd = eye(size(Ltint, 1), 'like', Ltint) - Ltint * (V_bd' * A00_U_bd);

% Hutchinson trace estimation
Z = 2 * (rand(n0, m) >= 0.5) - 1;% Generate m Rademacher random vectors
W = solve_A(Z, L00, U00, p00, q00);

reduced_rhs = Ltint * (V_bd' * W);
Z00z = W + A00_U_bd * (K_bd \ reduced_rhs);

quad_forms = sum(conj(Z) .* Z00z, 1);
trZ00 = mean(quad_forms);
end

function x = solve_A(b, L, U, p, q)
    b_permuted = b(p, :);
    y = L \ b_permuted;
    z = U \ y;
    x = zeros(size(z), 'like', z);
    x(q, :) = z;
end