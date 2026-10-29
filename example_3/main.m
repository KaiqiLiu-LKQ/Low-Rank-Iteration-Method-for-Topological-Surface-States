% ------------------------------------------------------------------------
% Main script for computing SDOS 
% ------------------------------------------------------------------------
clear all
close all
clc

% Parameters
a = 1;                                % lattice constant
c = 3e8;                              % speed of light
nf = 201;                              % number of frequency
nk = 201;                              % number of wavevector
freq = c/a*linspace(0.35, 0.60, nf);  % normalized frequency
omega = 2*pi*freq;                    % angular frequency
k = linspace(-pi, pi, nk)/a;          % wave vector
phase = exp(1i*k*a);                  % bloch phase
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
noIdx_cellNeg2=find(no2xyz(2,:)>=-3*a/sqrt(3)/4+a*sqrt(3)/2+posTol);
noIdx_cellNeg1=find(no2xyz(2,:)<=-3*a/sqrt(3)/4+a*sqrt(3)/2+posTol & no2xyz(2,:)>=-3*a/sqrt(3)/4+posTol);
noIdx_cell0=find(abs(no2xyz(2,:)-(-3*a/sqrt(3)/4))<posTol);
noIdx_cellPos1=find(no2xyz(2,:)<=-3*a/sqrt(3)/4-posTol & no2xyz(2,:)>=-3*a/sqrt(3)/4-a*sqrt(3)/2-posTol);
noIdx_cellPos2=find(no2xyz(2,:)<=-3*a/sqrt(3)/4-a*sqrt(3)/2-posTol);

noIdx_cellNeg2TopEdge=find(abs(no2xyz(2,:)-(-3*a/sqrt(3)/4+a*sqrt(3)/2+a*sqrt(3)/2))<posTol);
noIdx_cellPos1BotEdge=find(abs(no2xyz(2,:)-(-3*a/sqrt(3)/4-a*sqrt(3)/2))<posTol);

% Match node indexes in negative direction
[~,noIdx_tmpCellNeg]=min(pdist2(no2xyz(:,noIdx_cellNeg2)',no2xyz(:,noIdx_cellNeg1)'+[-a/2,a*sqrt(3)/2])');
noIdx_cellNeg1=noIdx_cellNeg1(:,noIdx_tmpCellNeg);
noIdx_MapCellNeg=[noIdx_cellNeg2;noIdx_cellNeg1];

% Match node indexes in positive direction
[~,noIdx_tmpCellPos]=min(pdist2(no2xyz(:,noIdx_cellPos1)',no2xyz(:,noIdx_cellPos2)'+[-a/2,a*sqrt(3)/2])');
noIdx_cellPos2=noIdx_cellPos2(:,noIdx_tmpCellPos);
noIdx_MapCellPos=[noIdx_cellPos1;noIdx_cellPos2];

% Reconstruct interior node indexes of cell-1 according to those of cell-2
noIdx_cellNeg2Int=setdiff(noIdx_cellNeg2,noIdx_cellNeg2TopEdge);
[~,tmp]=ismember(noIdx_cellNeg2Int,noIdx_MapCellNeg(1,:));
noIdx_cellNeg1Int=noIdx_MapCellNeg(2,tmp);

% Reconstruct interior node indexes of cell+2 acoording to those of cell+1
noIdx_cellPos1Int=setdiff(noIdx_cellPos1,noIdx_cellPos1BotEdge);
[~,tmp]=ismember(noIdx_cellPos1Int,noIdx_MapCellPos(1,:));
noIdx_cellPos2Int=noIdx_MapCellPos(2,tmp);

% Reconstruct top edge node indexes of cell-1 according to those of cell-2
[~,tmp]=ismember(noIdx_cellNeg2TopEdge,noIdx_MapCellNeg(1,:));
noIdx_cellNeg1TopEdge=noIdx_MapCellNeg(2,tmp);

% Reconstruct bottom edge node indexes of cell+2 acoording to those of cell+1
[~,tmp]=ismember(noIdx_cellPos1BotEdge,noIdx_MapCellPos(1,:));
noIdx_cellPos2BotEdge=noIdx_MapCellPos(2,tmp);

% Reconstruct node positions
noIdx_tmpAll=[noIdx_cellNeg2TopEdge,noIdx_cellNeg2Int,noIdx_cellNeg1TopEdge,noIdx_cellNeg1Int,noIdx_cell0,noIdx_cellPos1Int,noIdx_cellPos1BotEdge,noIdx_cellPos2Int,noIdx_cellPos2BotEdge];
no2xyz=no2xyz(:,noIdx_tmpAll);

% Reconstruct mass and stiffness matrices
mtrE=mtrE(noIdx_tmpAll(1:end),noIdx_tmpAll(1:end));
mtrK=mtrK(noIdx_tmpAll(1:end),noIdx_tmpAll(1:end));

% Node indexes for left and right edges (PBC)
noIdx_leftEdge=find(abs(no2xyz(2,:)-(-sqrt(3)*no2xyz(1,:)-sqrt(3)/2*a))<posTol);
noIdx_rightEdge=find(abs(no2xyz(2,:)-(-sqrt(3)*no2xyz(1,:)+sqrt(3)/2*a))<posTol);
[~,noIdx_tmpEdge]=min(pdist2(no2xyz(:,noIdx_leftEdge)',no2xyz(:,noIdx_rightEdge)'-[1*a,0])');
noIdx_rightEdge=noIdx_rightEdge(:,noIdx_tmpEdge);

% Define DOFs for truncating matrices
dof_cellNeg2LeftEdge=length(find(abs(no2xyz(2,:)-(-sqrt(3)*no2xyz(1,:)-sqrt(3)/2*a))<posTol & no2xyz(2,:)>=-3*a/sqrt(3)/4+a*sqrt(3)/2+posTol));
dof_cellNeg1LeftEdge=dof_cellNeg2LeftEdge;
dof_cellPos2LeftEdge=length(find(abs(no2xyz(2,:)-(-sqrt(3)*no2xyz(1,:)-sqrt(3)/2*a))<posTol & no2xyz(2,:)<=-3*a/sqrt(3)/4-a*sqrt(3)/2-posTol));
dof_cellPos1LeftEdge=dof_cellPos2LeftEdge;

dof_cellNeg2Int=length(noIdx_cellNeg2)-dof_cellNeg2LeftEdge;
dof_cellNeg1Int=length(noIdx_cellNeg1)-dof_cellNeg1LeftEdge;
dof_cell0Int=length(noIdx_cell0)-1;
dof_cellPos1Int=length(noIdx_cellPos1)-dof_cellPos1LeftEdge;
dof_cellPos2Int=length(noIdx_cellPos2)-dof_cellPos2LeftEdge;

dof_all = length(no2xyz);
dof_leftEdge = length(noIdx_leftEdge);

% Compute SDOS 
for mf = 1:nf
    for mk = 1:nk

        % Initialize the matrix for imposing boundary conditions
        Pn = eye(dof_all); 

        % Impose periodic boundary conditions
        for i=1:dof_leftEdge
            Pn(noIdx_leftEdge(i),noIdx_leftEdge(i))=0;
            Pn(noIdx_leftEdge(i),noIdx_rightEdge(i))=phase(mk);
        end

        % Reduce matrix dimension
        Pn = Pn(any(Pn, 2), any(Pn, 1));
        Pn = sparse(Pn);
        
        % Construct eliminated mass and stiffness matrices
        mtrKc=Pn'*mtrK*Pn;
        mtrEc=Pn'*mtrE*Pn;

        % Construct block Z matrix
        Z = (1-1i*eta)^2 * omega(mf)^2 * mtrEc - mtrKc;

        iNeg2 = dof_cell0Int + (1:dof_cellNeg2Int);
        iNeg1 = dof_cell0Int + dof_cellNeg2Int + (1:dof_cellNeg1Int);
        iPos1 = dof_cellNeg2Int + dof_cellNeg1Int + dof_cell0Int + (1:dof_cellPos1Int);
        iPos2 = dof_cellNeg2Int + dof_cellNeg1Int + dof_cell0Int + dof_cellPos1Int + (1:dof_cellPos2Int);
        iA11 = iNeg1;%   A11 corresponds to [Neg1; cell0]
        iA22 = iPos1;%   A22 corresponds to Pos1
        B11 = Z(iNeg2, iNeg1);%   B11: Neg2 row -> Neg1 column
        C11 = Z(iNeg1, iNeg2);%   C11: Neg1 row -> Neg2 column
        A11 = Z(iA11, iA11);%   Left surface block, including the cell0 interface DOFs

        % Right periodic tail coupling:
        A22 = Z(iA22, iA22);
        B22 = Z(iPos1, iPos2);%   B22: Pos1 row -> Pos2 column
        C22 = Z(iPos2, iPos1);%   C22: Pos2 row -> Pos1 column

        % Coupling between the two interface blocks
        B12 = Z(iA11, iA22);
        C21 = Z(iA22, iA11);

        % Implement LIM
        G00 = limSGF(A11, B11, C11,A22, B22, C22,B12, C21);       
        SDOS_lim(mf,mk) = omega(mf) * imag(trace(G00));
    end
    disp(mf);
end