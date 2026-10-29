clear all
close all
clc

%% 基本参数
ax = 3e-3;                           % lattice constant
ay = 3e-3;                           % lattice constant
az = 6e-3;                           % lattice constant
kx0 = 30;                            % kx0 in Comsol
ky0 = 40;                            % ky0 in Comsol
kz0 = 50;                            % kz0 in Comsol
c = 3e8;                             % speed of light
GHz = 1e9;                           % 1GHz
nf = 61;                             % number of frequency
nk = 61;                             % number of wavevector
freq = GHz*linspace(0.1, 0.32, nf);  % normalized frequency
omega = 2*pi*freq;                    % angular frequency
eta = 1e-3;                           % imaginary part of frequency
posTol = ax*1e-3;                     % position tolerance
kType = 0;                            % swtich sweep k direction
[phaseSweep.phasekx,phaseSweep.phaseky,kx,ky] = KsweepDir(kType,nk,ax,ay);

%% Load FEM information and require COMSOL-MATLAB environment
% model = All_PEC_half_inf;
model = example;
% Extract mass and stiffness matrices
mtrInfo = mphmatrix(model,'sol1','out',{'K','E','Kc','Ec','Null','Nullf'},'initmethod','init');
% Extract mesh information
meshInfo = mphxmeshinfo(model);
no2xyz = meshInfo.dofs.coords;

% Find the boundary position from null info
phasekx0 = exp(-1i*kx0*ax);
phaseky0 = exp(-1i*ky0*ay);
[IndexPECInside,IndexBoundaryOutSide] = FindBoundaryByNull3D(mtrInfo.Null,no2xyz,phasekx0,phaseky0,length(mtrInfo.K),1e-7);

% Find the sequence numbers of the periodic boundary edges
[IndexCell,dofCell,mtrE,mtrK,no2xyz] = ReconstruCell3D(mtrInfo.E,mtrInfo.K,no2xyz,IndexPECInside,ax,ay,az,posTol);

% Compute SDOS 
SDOS_lim=zeros(nf,nk);
for mk = 1:nk
    % Initialize the matrix for imposing boundary conditions
    Pn = buildPn(dofCell,IndexCell,phaseSweep, mk);
    % Construct eliminated mass and stiffness matrices
    mtrKc = Pn' * mtrK(IndexCell.total_bulk,IndexCell.total_bulk) * Pn;
    mtrEc = Pn' * mtrE(IndexCell.total_bulk,IndexCell.total_bulk) * Pn;
    for mf = 1:nf
        % Construct block Z matrix
        Z = (1-1i*eta)^2 * omega(mf)^2 * mtrEc - mtrKc;

        noIdx_cellZ0=1:dofCell.cellInt;
        noIdx_cellZ1=dofCell.cellInt+(1:dofCell.cellInt);

        Z00 = Z(noIdx_cellZ0, noIdx_cellZ0);
        Z01 = Z(noIdx_cellZ0, noIdx_cellZ1);
        Z10 = Z(noIdx_cellZ1, noIdx_cellZ0);
        Z11 = Z(noIdx_cellZ1, noIdx_cellZ1);

        G00_LIM  = limSGF(Z00,Z01,Z10);
        SDOS_lim(mf,mk) = omega(mf) * imag(G00_LIM);

    end
    disp(mk);
end

%%
function [phasekx,phaseky,kx,ky] = KsweepDir(kType,nk,ax,ay)
    switch kType
        case 0                                    %sweep kx
            nkx = nk;                               % number of wavevector kx
            nky = 1;                                % number of wavevector ky
            ky = 0*pi/ay;                                 % wave vector ky
            kx = linspace(-pi, pi, nkx)/ax;          % wave vector kx

        case 1                                    %sweep ky
            nkx = 1;                                % number of wavevector kx
            nky = nk;                               % number of wavevector ky
            kx = 0;                                 % wave vector kx
            ky = linspace(-pi, pi, nky)/ay;          % wave vector ky

        case 2                                      %sweep Tao to M
            nkT = nk;
            kx = linspace(-pi, pi, nkT)/ax;
            ky = 0.3*pi/ay;
    
    end
        phasekx = exp(1i*kx*ax)+0*exp(1i*ky*ay);               % bloch phase kx
        phaseky = 0*exp(1i*kx*ax)+exp(1i*ky*ay);               % bloch phase kx
end



