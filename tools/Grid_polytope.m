function xgrid = Grid_polytope(P,n)

    Pbox = interval(P);

    lb = Pbox.inf;
    ub = Pbox.sup;

    xdim = size(lb,1);

    xpoints = cell(1,xdim);

    for j = 1:xdim  
        xpoints{j} = linspace(lb(j),ub(j),n);
    end

    grids = cell(1,xdim);
    [grids{:}] = ndgrid(xpoints{:});

    xgrid = zeros(xdim,n^xdim);
    for j = 1:xdim
        xgrid(j,:) = grids{j}(:)';
    end
    
    inside = contains(P, xgrid, 'exact', 1e-7);

    xgrid = xgrid(:, inside);

end