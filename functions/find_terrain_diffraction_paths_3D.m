% Inputs format:
%   - xTerrain = [row]
%   - yTerrain = [row]
%   - Zterrain = [matrix: roxs = y, columns = x]
% 
%   - S = [xs, 0, zs];
%   - R = [xr, 0, zr];

function paths = find_terrain_diffraction_paths_3D(S, R, xTerrain, yTerrain, Zterrain)

    xTerrain = xTerrain(:).';
    yTerrain = yTerrain(:);

    if size(Zterrain,1) ~= length(yTerrain) || ...
       size(Zterrain,2) ~= length(xTerrain)

        error(['Zterrain must have size ' ...
               'length(yTerrain) x length(xTerrain).']);

    end

    if abs(S(2)) > 1e-9 || abs(R(2)) > 1e-9
        error('S and R must lie on y = 0.');
    end

    if S(1) >= R(1)
        error('The function assumes S is before R along the x-axis.');
    end

    d = norm(R - S);

    paths.top = findTopPath( ...
        S, R, xTerrain, yTerrain, Zterrain, d);

    paths.left = findSidePath( ...
        S, R, xTerrain, yTerrain, Zterrain, d, +1);

    paths.right = findSidePath( ...
        S, R, xTerrain, yTerrain, Zterrain, d, -1);

end


function path = findTopPath(S, R, xTerrain, yTerrain, Zterrain, d)

    inside = xTerrain > S(1) & xTerrain < R(1);

    xProfile = xTerrain(inside);

    zProfile = interp2( ...
        xTerrain, ...
        yTerrain, ...
        Zterrain, ...
        xProfile, ...
        zeros(size(xProfile)), ...
        'linear');

    terrainPoints = [
        xProfile(:), ...
        zProfile(:)
    ];

    S2 = [S(1), S(3)];
    R2 = [R(1), R(3)];

    nodes = [
        S2
        terrainPoints
        R2
    ];

    nNodes = size(nodes,1);

    visible = false(nNodes,nNodes);

    for i = 1:nNodes-1

        for j = i+1:nNodes

            x1 = nodes(i,1);
            z1 = nodes(i,2);

            x2 = nodes(j,1);
            z2 = nodes(j,2);

            between = ...
                xProfile > x1 & ...
                xProfile < x2;

            xCheck = xProfile(between);

            if isempty(xCheck)

                visible(i,j) = true;

            else

                zGround = interp1( ...
                    xProfile, ...
                    zProfile, ...
                    xCheck, ...
                    'linear');

                zLine = z1 + ...
                        (z2-z1) .* ...
                        (xCheck-x1) ./ ...
                        (x2-x1);

                tol = 1e-9;

                visible(i,j) = ...
                    all(zLine >= zGround - tol);

            end

        end

    end

    pathIndices = shortestVisiblePath(nodes, visible);

    path2D = nodes(pathIndices,:);

    pathPoints = [
        path2D(:,1), ...
        zeros(size(path2D,1),1), ...
        path2D(:,2)
    ];

    path = buildPathOutput(pathPoints, d);

end


function path = findSidePath( ...
    S, R, xTerrain, yTerrain, Zterrain, d, sideSign)

    [X,Y] = meshgrid(xTerrain,yTerrain);

    insideX = X > S(1) & X < R(1);

    if sideSign > 0

        sideMask = Y >= 0;

    else

        sideMask = Y <= 0;

    end

    mask = insideX & sideMask;

    terrainPoints = [
        X(mask), ...
        Y(mask), ...
        Zterrain(mask)
    ];

    [~,order] = sortrows(terrainPoints,[1 2]);

    terrainPoints = terrainPoints(order,:);

    nodes = [
        S
        terrainPoints
        R
    ];

    nNodes = size(nodes,1);

    visible = false(nNodes,nNodes);

    for i = 1:nNodes-1

        for j = i+1:nNodes

            if nodes(j,1) <= nodes(i,1)
                continue
            end

            if segmentAboveTerrain( ...
                    nodes(i,:), ...
                    nodes(j,:), ...
                    xTerrain, ...
                    yTerrain, ...
                    Zterrain)

                visible(i,j) = true;

            end

        end

    end

    pathIndices = shortestVisiblePath(nodes, visible);

    pathPoints = nodes(pathIndices,:);

    path = buildPathOutput(pathPoints, d);

end


function isVisible = segmentAboveTerrain( ...
    P1, P2, xTerrain, yTerrain, Zterrain)

    segmentLength = norm(P2 - P1);

    dx = median(diff(xTerrain));

    if length(yTerrain) > 1
        dy = median(abs(diff(yTerrain)));
    else
        dy = dx;
    end

    ds = min(dx,dy) / 2;

    nSamples = max( ...
        ceil(segmentLength / ds), ...
        2);

    t = linspace(0,1,nSamples);

    x = P1(1) + ...
        t .* (P2(1)-P1(1));

    y = P1(2) + ...
        t .* (P2(2)-P1(2));

    z = P1(3) + ...
        t .* (P2(3)-P1(3));

    zGround = interp2( ...
        xTerrain, ...
        yTerrain, ...
        Zterrain, ...
        x, ...
        y, ...
        'linear', ...
        NaN);

    valid = ~isnan(zGround);

    tol = 1e-9;

    isVisible = ...
        all(z(valid) >= zGround(valid) - tol);

end


function pathIndices = shortestVisiblePath(nodes, visible)

    nNodes = size(nodes,1);

    dist = inf(nNodes,1);
    previous = zeros(nNodes,1);

    dist(1) = 0;

    for j = 2:nNodes

        for i = 1:j-1

            if ~visible(i,j)
                continue
            end

            if isinf(dist(i))
                continue
            end

            segmentLength = ...
                norm(nodes(j,:) - nodes(i,:));

            newDistance = ...
                dist(i) + segmentLength;

            if newDistance < dist(j)

                dist(j) = newDistance;
                previous(j) = i;

            end

        end

    end

    if isinf(dist(end))
        error('No valid diffraction path was found.');
    end

    current = nNodes;

    pathIndices = current;

    while current ~= 1

        current = previous(current);

        if current == 0
            error('Could not reconstruct the diffraction path.');
        end

        pathIndices = [
            current, ...
            pathIndices
        ];

    end

end


function path = buildPathOutput(pathPoints, d)

    nPoints = size(pathPoints,1);

    path.points = pathPoints;

    if nPoints <= 2

        path.edges = zeros(0,3);

        path.d_ss = 0;
        path.d_sr = 0;
        path.e = 0;

        path.L_diff = d;
        path.z = 0;

        return

    end

    path.edges = ...
        pathPoints(2:end-1,:);

    path.d_ss = ...
        norm(pathPoints(2,:) - pathPoints(1,:));

    path.d_sr = ...
        norm(pathPoints(end,:) - pathPoints(end-1,:));

    if nPoints == 3

        path.e = 0;

    else

        edgeSegments = ...
            diff(path.edges,1,1);

        path.e = ...
            sum(sqrt(sum(edgeSegments.^2,2)));

    end

    segmentVectors = ...
        diff(pathPoints,1,1);

    segmentLengths = ...
        sqrt(sum(segmentVectors.^2,2));

    path.L_diff = ...
        sum(segmentLengths);

    path.z = ...
        path.L_diff - d;

end