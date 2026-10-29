function [noIdx_leftSurf,noIdx_rightSurf,...
          noIdy_leftSurf,noIdy_rightSurf,...
          noIdxy_00edge,noIdxy_01edge,...
          noIdxy_10edge,noIdxy_11edge]=FindPBCboundary3D(no2xyz,ax,ay,az,posTol)
    %% Node indexes for edge (PBC)
    noIdxy_00edge = find(abs(no2xyz(1,:)-(-0.5*ax))<=posTol & abs(no2xyz(2,:)-(-0.5*ay))<=posTol & no2xyz(3,:)<=0.5*az-posTol);%00角点
    noIdxy_11edge = find(abs(no2xyz(1,:)-(0.5*ax))<=posTol & abs(no2xyz(2,:)-(0.5*ay))<=posTol & no2xyz(3,:)<=0.5*az-posTol);%11角点
    noIdxy_01edge = find(abs(no2xyz(1,:)-(-0.5*ax))<=posTol & abs(no2xyz(2,:)-(0.5*ay))<=posTol & no2xyz(3,:)<=0.5*az-posTol);%01角点
    noIdxy_10edge = find(abs(no2xyz(1,:)-(0.5*ax))<=posTol & abs(no2xyz(2,:)-(-0.5*ay))<=posTol & no2xyz(3,:)<=0.5*az-posTol);%10角点
    %按照00排序对应
    [~, noIdxy_tmp11edge] = min( pdist2(no2xyz(:,noIdxy_00edge)', no2xyz(:,noIdxy_11edge)'-[ax,ay,0])' );
    [~, noIdxy_tmp01edge] = min( pdist2(no2xyz(:,noIdxy_00edge)', no2xyz(:,noIdxy_01edge)'-[0,ay,0])' );
    [~, noIdxy_tmp10edge] = min( pdist2(no2xyz(:,noIdxy_00edge)', no2xyz(:,noIdxy_10edge)'-[ax,0,0])' );

    noIdxy_11edge = noIdxy_11edge(:,noIdxy_tmp11edge);
    noIdxy_01edge = noIdxy_01edge(:,noIdxy_tmp01edge);
    noIdxy_10edge = noIdxy_10edge(:,noIdxy_tmp10edge);

    %% Node indexes for left and right edges (PBC)
    noIdx_leftSurf = find( abs(no2xyz(1,:)-(-0.5*ax))<posTol & abs(no2xyz(2,:)-(-0.5*ay))>=posTol & abs(no2xyz(2,:)-(0.5*ay))>=posTol & no2xyz(3,:)<=0.5*az-posTol );
    noIdx_rightSurf = find( abs(no2xyz(1,:)-(0.5*ax))<posTol & abs(no2xyz(2,:)-(-0.5*ay))>=posTol & abs(no2xyz(2,:)-(0.5*ay))>=posTol & no2xyz(3,:)<=0.5*az-posTol );
    [~, noIdx_tmpSurf] = min( pdist2(no2xyz(:,noIdx_leftSurf)', no2xyz(:,noIdx_rightSurf)'-[ax,0,0])' );
    noIdx_rightSurf = noIdx_rightSurf(:, noIdx_tmpSurf);


    %% Node indexes for left and right edges (PBC)
    noIdy_leftSurf = find(abs(no2xyz(1,:)-(-0.5*ax))>=posTol & abs(no2xyz(1,:)-(0.5*ay))>=posTol & abs(no2xyz(2,:)-(-0.5*ay))<posTol & no2xyz(3,:)<=0.5*az-posTol );
    noIdy_rightSurf = find(abs(no2xyz(1,:)-(-0.5*ax))>=posTol & abs(no2xyz(1,:)-(0.5*ay))>=posTol & abs(no2xyz(2,:)-(0.5*ay))<posTol & no2xyz(3,:)<=0.5*az-posTol );
    [~, noIdy_tmpSurf] = min( pdist2(no2xyz(:,noIdy_leftSurf)', no2xyz(:,noIdy_rightSurf)'-[0,az,0])' );
    noIdy_rightSurf = noIdy_rightSurf(:, noIdy_tmpSurf);
   
end