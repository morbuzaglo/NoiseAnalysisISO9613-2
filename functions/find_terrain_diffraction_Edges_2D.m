% Inputs format:
%   - xTerrain = [0 10 20 30 40 50 60 70 80 90 100];
%   - zTerrain = [0  1  3  8 12  7  5 10  6  2   0];
%   - S = [0, 2];
%   - R = [100, 3];

% IMPORTANT: the function assumes SOURCE is to the left from the RECIEVER,
% and the [x,z] points are from left to right.

function [edges, pathPoints, pathLength, visible] = find_terrain_diffraction_Edges_2D(S, R, xTerrain, zTerrain)

    xTerrain = xTerrain(:);
    zTerrain = zTerrain(:);

    if length(xTerrain) ~= length(zTerrain)
        error('xTerrain and zTerrain must have the same length.');
    end

    if any(diff(xTerrain) <= 0)
        error('xTerrain must be strictly increasing.');
    end

    if S(1) >= R(1)
        error('This function assumes S is to the left of R.');
    end

    inside = xTerrain > S(1) & xTerrain < R(1);

    terrainPoints = [xTerrain(inside), zTerrain(inside)];

    nodes = [
        S
        terrainPoints
        R
    ];

    nNodes = size(nodes,1);

    visible = false(nNodes,nNodes);

    for i = 1:nNodes-1

        for j = i+1:nNodes

            x1 = nodes(i,1);
            z1 = nodes(i,2);

            x2 = nodes(j,1);
            z2 = nodes(j,2);

            between = xTerrain > x1 & xTerrain < x2;

            xCheck = xTerrain(between);
            zCheck = zTerrain(between);

            if isempty(xCheck)

                visible(i,j) = true;

            else

                zLine = z1 + ...
                        (z2-z1) .* ...
                        (xCheck-x1) ./ ...
                        (x2-x1);

                tol = 1e-9;

                if all(zLine >= zCheck - tol)
                    visible(i,j) = true;
                end

            end

        end

    end

    dist = inf(nNodes,1);
    previous = zeros(nNodes,1);

    dist(1) = 0;

    for j = 2:nNodes

        for i = 1:j-1

            if visible(i,j)

                segmentLength = norm(nodes(j,:) - nodes(i,:));

                newDistance = dist(i) + segmentLength;

                if newDistance < dist(j)

                    dist(j) = newDistance;
                    previous(j) = i;

                end

            end

        end

    end

    if isinf(dist(end))
        error('No valid propagation path was found above the terrain.');
    end

    pathIndices = nNodes;

    current = nNodes;

    while current ~= 1

        current = previous(current);

        if current == 0
            error('Could not reconstruct propagation path.');
        end

        pathIndices = [current, pathIndices];
    end

    pathPoints = nodes(pathIndices,:);

    if size(pathPoints,1) > 2
        edges = pathPoints(2:end-1,:);
    else
        edges = zeros(0,2);
    end

    pathLength = dist(end);

end