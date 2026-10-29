function [IndexCell,dofCell,mtrE,mtrK,no2xyz] = ReconstruCell3D(mtrE,mtrK,no2xyz,IndexPECInside,ax,ay,az,posTol)
    %% Find specific nodes indexes
    noIdx_cell0 = find(no2xyz(3,:)<=(0.5*az-posTol) & no2xyz(3,:)>=(-0.5*az-posTol));  % cell0 except its top edge (PEC)
    noIdx_cell1 = find(no2xyz(3,:)<=(-0.5*az-posTol) & no2xyz(3,:)>=(-1.5*az-posTol));  % cell1 except its top edge (PEC)
    noIdx_cell0_PECinside = IndexPECInside(no2xyz(3,IndexPECInside)<=(0.5*az-posTol) & no2xyz(3,IndexPECInside)>=(-0.5*az-posTol));
    noIdx_cell1_PECinside = IndexPECInside(no2xyz(3,IndexPECInside)<=(-0.5*az-posTol) & no2xyz(3,IndexPECInside)>=(-1.5*az-posTol));
    noIdx_cell0_bulk = setdiff(noIdx_cell0,noIdx_cell0_PECinside);
    noIdx_cell1_bulk = setdiff(noIdx_cell1,noIdx_cell1_PECinside);
    noIdx_cell0TopEdge = find(abs(no2xyz(3,:)-0.5*az)<posTol);  % cell0 top edge (PEC)

    %% Reorder all the numbers
    [~, noIdx_tmpCell] = min(pdist2(no2xyz(:,noIdx_cell0_bulk)', no2xyz(:,noIdx_cell1_bulk)'+[0,0,1*az])');
    noIdx_cell1_bulk = noIdx_cell1_bulk(:, noIdx_tmpCell);
    [~, noIdx_tmpCell] = min(pdist2(no2xyz(:,noIdx_cell0_PECinside)', no2xyz(:,noIdx_cell1_PECinside)'+[0,0,1*az])');
    noIdx_cell1_PECinside = noIdx_cell1_PECinside(:, noIdx_tmpCell); 

    %% Reconstruct node positions
    noIdx_tmpAll = [noIdx_cell0TopEdge,noIdx_cell0_PECinside,noIdx_cell0_bulk,...  %cell0
                                       noIdx_cell1_PECinside,noIdx_cell1_bulk];    %cell1
    no2xyz = no2xyz(:, noIdx_tmpAll);
    
    %% Reconstruct mass and stiffness matrices
    mtrE = mtrE(noIdx_tmpAll(1:end), noIdx_tmpAll(1:end));
    mtrK = mtrK(noIdx_tmpAll(1:end), noIdx_tmpAll(1:end));
    
    %% Node indexes for cell0 top edge (PEC)
    IndexCell.cell0TopEdge = 1:length(noIdx_cell0TopEdge);
    %% Node indexes outside cell0 top edge (PEC)
    IndexCell.cell0_PECinside = length(noIdx_cell0TopEdge) + (1:length(noIdx_cell0_PECinside));
    IndexCell.cell0_bulk = length(noIdx_cell0TopEdge) + length(noIdx_cell0_PECinside) + (1:length(noIdx_cell0_bulk));
    IndexCell.cell1_PECinside = length(noIdx_cell0TopEdge) + length(noIdx_cell0_PECinside) + length(noIdx_cell0_bulk) + (1:length(noIdx_cell1_PECinside));
    IndexCell.cell1_bulk = length(noIdx_cell0TopEdge) + length(noIdx_cell0_PECinside) + length(noIdx_cell0_bulk) + length(noIdx_cell1_PECinside) + (1:length(noIdx_cell1_bulk));
    IndexCell.total_bulk = [IndexCell.cell0_bulk,IndexCell.cell1_bulk];

    %% Find boundary Index
    [IndexCell.XleftSurf,IndexCell.XrightSurf,...
     IndexCell.YleftSurf,IndexCell.YrightSurf,...
     IndexCell.XY00edge,IndexCell.XY01edge,...
     IndexCell.XY10edge,IndexCell.XY11edge]=FindPBCboundary3D(no2xyz,ax,ay,az,posTol);

    %% number of dof
    dofCell.cell0TopEdge = length(IndexCell.cell0TopEdge);
    dofCell.cell0_PECinside = length(IndexCell.cell0_PECinside);
    dofCell.cell0_bulk = length(IndexCell.cell0_bulk);
    dofCell.cell1_PECinside = length(IndexCell.cell1_PECinside);
    dofCell.cell1_bulk = length(IndexCell.cell1_bulk);
    dofCell.XleftSurf = length(IndexCell.XleftSurf);
    dofCell.XrightSurf = length(IndexCell.XrightSurf);
    dofCell.YleftSurf = length(IndexCell.YleftSurf);
    dofCell.YrightSurf = length(IndexCell.YrightSurf);
    dofCell.XY00edge = length(IndexCell.XY00edge);
    dofCell.XY01edge = length(IndexCell.XY01edge);
    dofCell.XY10edge = length(IndexCell.XY10edge);
    dofCell.XY11edge = length(IndexCell.XY11edge);
    
    dofCell.total_bulk = length(IndexCell.total_bulk);


    dofCell.cellInt = dofCell.cell0_bulk-dofCell.XY00edge/2*3-dofCell.XleftSurf/2-dofCell.YleftSurf/2; %实际计算的单cell自由度，单边加体
    dofCell.All = length(no2xyz);
end