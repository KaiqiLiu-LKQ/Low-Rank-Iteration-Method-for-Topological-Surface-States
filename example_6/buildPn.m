function Pn = buildPn(dofCell, IndexCell, phaseSweep, mk)

    N   = dofCell.All;

    nTop = dofCell.cell0TopEdge;
    nPEC = dofCell.cell0_PECinside;   
    nX   = dofCell.XrightSurf;
    nY   = dofCell.YrightSurf;
    nXY  = dofCell.XY00edge;


    px = phaseSweep.phasekx(mk);
    py = phaseSweep.phaseky(mk);

    col = @(v,n) reshape(v(1:n), [], 1);


    zeroedIdx = [
        col(IndexCell.cell0TopEdge,    nTop);  
        col(IndexCell.cell0_PECinside, nPEC);   
        col(IndexCell.cell1_PECinside, nPEC);   
        col(IndexCell.XrightSurf,      nX );    
        col(IndexCell.YrightSurf,      nY );    
        col(IndexCell.XY01edge,        nXY);    
        col(IndexCell.XY10edge,        nXY);
        col(IndexCell.XY11edge,        nXY);
    ];

    keepDiag = setdiff((1:N).', zeroedIdx);

    rows = [
        keepDiag;                          
        col(IndexCell.XrightSurf, nX );   
        col(IndexCell.YrightSurf, nY );  
        col(IndexCell.XY01edge,   nXY);   
        col(IndexCell.XY10edge,   nXY);   
        col(IndexCell.XY11edge,   nXY);   
    ];

    cols = [
        keepDiag;
        col(IndexCell.XleftSurf, nX );
        col(IndexCell.YleftSurf, nY );
        col(IndexCell.XY00edge,  nXY);
        col(IndexCell.XY00edge,  nXY);
        col(IndexCell.XY00edge,  nXY);
    ];

    vals = [
        ones(numel(keepDiag), 1);
        px      * ones(nX , 1);
        py      * ones(nY , 1);
        py      * ones(nXY, 1);
        px      * ones(nXY, 1);
        (px*py) * ones(nXY, 1);
    ];


    Pn = sparse(rows, cols, vals, N, N);
    Pn = Pn(any(Pn, 2), any(Pn, 1));
end