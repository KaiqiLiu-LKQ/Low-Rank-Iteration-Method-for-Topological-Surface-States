% ------------------------------------------------------------------------
% Main script for computing SDOS 
% ------------------------------------------------------------------------
clear all
close all
clc

% Parameters
a = 1;                                % lattice constant
c = 3e8;                              % speed of light
nf = 201;                             % number of frequency
nk = 201;                             % number of wavevector
freq = c/a*linspace(0.45, 0.80, nf);  % normalized frequency
omega = 2*pi*freq;                    % angular frequency
k = linspace(-pi, pi, nk)/a;          % wave vector
phase = exp(1i*k*a);                  % bloch phase
SDOS_crm = zeros(nf, nk);             % surface density of states
SDOS_tmm = zeros(nf, nk);             % surface density of states
SDOS_lim = zeros(nf, nk);             % surface density of states
eta = 1e-3;                           % imaginary part of frequency
posTol = a*1e-3;                      % position tolerance


% Load FEM information and require COMSOL-MATLAB environment
% Extract mass and stiffness matrices
model = example;
mtrInfo = mphmatrix(model,'sol1','out',{'K','E','Kc','Ec'},'initmethod','init');
mtrE = (mtrInfo.E);
mtrK = (mtrInfo.K);

% Extract mesh information
meshInfo = mphxmeshinfo(model);
no2xyz = meshInfo.nodes.coords;

% Find specific nodes indexes
noIdx_cell0 = find( no2xyz(2,:)<=(0.5*a-posTol) & no2xyz(2,:)>=(-0.5*a-posTol) );  % cell0 except its top edge (PEC)
noIdx_cell1 = find( no2xyz(2,:)<=(-0.5*a-posTol) & no2xyz(2,:)>=(-1.5*a-posTol) );  % cell1 except its top edge (PEC)
noIdx_cell0TopEdge = find( abs(no2xyz(2,:)-0.5*a)<posTol );  % cell0 top edge (PEC)

% Reconstruct node indexes
[~, noIdx_tmpCell] = min(pdist2(no2xyz(:,noIdx_cell0)', no2xyz(:,noIdx_cell1)'+[0,1*a])');
noIdx_cell1 = noIdx_cell1(:, noIdx_tmpCell);

% Reconstruct node positions
noIdx_tmpAll = [noIdx_cell0TopEdge, noIdx_cell0, noIdx_cell1];
no2xyz = no2xyz(:, noIdx_tmpAll);  

% Reconstruct mass and stiffness matrices
mtrE = mtrE(noIdx_tmpAll(1:end), noIdx_tmpAll(1:end)); 
mtrK = mtrK(noIdx_tmpAll(1:end), noIdx_tmpAll(1:end));

% Node indexes for left and right edges (PBC)
noIdx_leftEdge = find( abs(no2xyz(1,:)-(-0.5*a))<posTol & no2xyz(2,:)<=0.5*a-posTol );
noIdx_rightEdge = find( abs(no2xyz(1,:)-(0.5*a))<posTol & no2xyz(2,:)<=0.5*a-posTol );
[~, noIdx_tmpEdge] = min( pdist2(no2xyz(:,noIdx_leftEdge)', no2xyz(:,noIdx_rightEdge)'-[a,0])' );
noIdx_rightEdge = noIdx_rightEdge(:, noIdx_tmpEdge);

% Node indexes for cell0 top edge (PEC)
noIdx_cell0TopEdge = 1:length(noIdx_cell0TopEdge);

% Node indexes outside cell0 top edge (PEC)
noIdx_int = length(noIdx_cell0TopEdge)+1:length(no2xyz);

% Define DOFs for truncating matrices
dof_cell0LeftEdge = length( find( abs(no2xyz(1,:)-(-0.5*a))<posTol & no2xyz(2,:)<=0.5*a-posTol & no2xyz(2,:)>=-0.5*a-posTol ) );  % DOFs of cell0 left edge 
dof_cell1LeftEdge = length( find( abs(no2xyz(1,:)-(-0.5*a))<posTol & no2xyz(2,:)<=-0.5*a-posTol & no2xyz(2,:)>=-1.5*a-posTol ) );  % DOFs of cell1 left edge
dof_cell0Int = length(noIdx_cell0) - dof_cell0LeftEdge;  % DOFs of reduced cell0
dof_cell1Int = length(noIdx_cell1) - dof_cell1LeftEdge;  % DOFs of reduced cell1
dof_all = length(no2xyz); % DOFs of cell0 and cell1 

dof_cell0TopEdge = length(noIdx_cell0TopEdge); % DOFs of cell0 top edge (PEC)
dof_leftEdge = length(noIdx_leftEdge); % DOFs of left edge (PBC)

% Compute SDOS 
for mf = 1:nf
    for mk = 1:nk
        
        % Initialize the matrix for imposing boundary conditions
        Pn = eye(dof_all); 
        
        % Impose PEC boundary condition
        for i = 1:dof_cell0TopEdge
            Pn(noIdx_cell0TopEdge(i), noIdx_cell0TopEdge(i)) = 0;
        end
        
        % Impose periodic boundary conditions
        for i = 1:dof_leftEdge
            Pn(noIdx_leftEdge(i), noIdx_leftEdge(i)) = 0;
            Pn(noIdx_leftEdge(i), noIdx_rightEdge(i)) = phase(mk);
        end
        
        % Reduce matrix dimension
        Pn = Pn(any(Pn, 2), any(Pn, 1));
        Pn = sparse(Pn);
        
        % Construct eliminated mass and stiffness matrices
        mtrKc = Pn' * mtrK(noIdx_int, noIdx_int) * Pn;
        mtrEc = Pn' * mtrE(noIdx_int, noIdx_int) * Pn;
        
        % Construct block Z matrix
        Z = (1-1i*eta)^2 * omega(mf)^2 * mtrEc - mtrKc;
        Z00 = Z(1:dof_cell0Int, 1:dof_cell0Int);
        Z01 = Z(1:dof_cell0Int, dof_cell0Int+1:dof_cell0Int+dof_cell1Int);
        Z10 = Z(dof_cell0Int+1:dof_cell0Int+dof_cell1Int, 1:dof_cell0Int);
        
        % Implement CRM TMM and LIM
        G00_crm = crmSGF(Z00,Z01,Z10); 
        G00_tmm = tmmSGF(Z00,Z01,Z10);
        G00_lim = limSGF(Z00,Z01,Z10);
        SDOS_crm(mf,mk) = omega(mf) * imag(trace(G00_crm));
        SDOS_tmm(mf,mk) = omega(mf) * imag(trace(G00_tmm));
        SDOS_lim(mf,mk) = omega(mf) * imag(G00_lim);
    end
    disp(mf);
end