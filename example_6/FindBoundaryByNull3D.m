function [IndexPECInside,IndexBoundaryOutSide] = FindBoundaryByNull3D(Null,no2xyz,phasekx,phaseky,DofTotal,phase0tol)
    %% Removed PEC label
    [row_indices, ~, ~] = find(Null);
    non_zero_rows = unique(row_indices);
    dof_all = 1:DofTotal;
    PEC_nodes = setdiff(dof_all,non_zero_rows);
    %% Pending cycle edge strip
    [row_PBC_all, col_PBC_all, ~] = find(Null==1);

    [row_XY00,col_XY00, ~] = find(Null~=1 & Null~=0 & imag(Null)==0);
    norm_coner = Null(row_XY00,col_XY00);
    [row_XY_bourndary,~, ~] = find(abs(Null-phasekx*phaseky*norm_coner(1))<phase0tol | abs(Null+phasekx*phaseky*norm_coner(1))<phase0tol);
    
    [row_X_bourndaryR,col_X_bourndaryR, ~] = find(abs(Null-phasekx)<phase0tol | abs(Null+phasekx)<phase0tol | abs(Null-phasekx*norm_coner(1))<phase0tol | abs(Null+phasekx*norm_coner(1))<phase0tol);
    [row_Y_bourndaryR,col_Y_bourndaryR, ~] = find(abs(Null-phaseky)<phase0tol | abs(Null+phaseky)<phase0tol | abs(Null-phaseky*norm_coner(1))<phase0tol | abs(Null+phaseky*norm_coner(1))<phase0tol);

    row_X_bourndaryL = row_PBC_all(ismember(col_PBC_all,col_X_bourndaryR));
    row_Y_bourndaryL = row_PBC_all(ismember(col_PBC_all,col_Y_bourndaryR));

    %% Extract data based on the coordinate system's extreme values, default is square grids
    Xmax = max(no2xyz(1,:));Xmin = min(no2xyz(1,:));
    Ymax = max(no2xyz(2,:));Ymin = min(no2xyz(2,:));
    Zmax = max(no2xyz(3,:));Zmin = min(no2xyz(3,:));
    
    IndexPECInside = PEC_nodes(no2xyz(3,PEC_nodes) < (Zmax-Zmax*1e-2) & no2xyz(3,PEC_nodes) > (Zmin+Zmax*1e-2));
    IndexPECOutside = setdiff(PEC_nodes,IndexPECInside);
    IndexBoundaryOutSide = [IndexPECOutside,row_X_bourndaryR',row_X_bourndaryL',row_Y_bourndaryR',row_Y_bourndaryL',row_XY_bourndary'];
end