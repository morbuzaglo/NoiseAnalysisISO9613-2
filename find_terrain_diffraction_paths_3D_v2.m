function paths = find_terrain_diffraction_paths_3D_v2( ...
    S,R,xTerrain,yTerrain,Zterrain,varargin)

% ============================================================
% FAST 3D TERRAIN DIFFRACTION PATH SEARCH - V2
%
% INPUTS:
%   S = [x y z] source coordinates
%   R = [x y z] receiver coordinates
%
%   xTerrain = terrain X vector
%   yTerrain = terrain Y vector
%   Zterrain = terrain elevation matrix
%
% OPTIONAL:
%   opts.topSpacing       [m]
%   opts.lateralSpacing   [m]
%   opts.maxLateralOffset [m]
%   opts.nLateralPaths
%   opts.clearanceTol     [m]
%
% OUTPUT:
%   paths.top
%   paths.left
%   paths.right
%
% Each path may contain:
%   .z
%   .d_ss
%   .d_sr
%   .e
%   .L_diff
%   .edges
%   .profile_s
%   .profile_z
%
% ============================================================

%% ------------------------------------------------------------
% OPTIONS
% ------------------------------------------------------------

opts = struct;

opts.topSpacing = [];
opts.lateralSpacing = [];
opts.maxLateralOffset = [];
opts.nLateralPaths = 6;
opts.clearanceTol = 0.05;

if ~isempty(varargin)
    userOpts = varargin{1};

    names = fieldnames(userOpts);

    for k = 1:length(names)
        opts.(names{k}) = userOpts.(names{k});
    end
end

%% ------------------------------------------------------------
% GRID RESOLUTION
% ------------------------------------------------------------

dxGrid = median(abs(diff(xTerrain)));
dyGrid = median(abs(diff(yTerrain)));

gridSpacing = min(dxGrid,dyGrid);

if isempty(opts.topSpacing)
    opts.topSpacing = gridSpacing;
end

if isempty(opts.lateralSpacing)
    opts.lateralSpacing = 2*gridSpacing;
end

%% ------------------------------------------------------------
% BASIC GEOMETRY
% ------------------------------------------------------------

Sxy = S(1:2);
Rxy = R(1:2);

deltaXY = Rxy-Sxy;

dHorizontal = norm(deltaXY);
dDirect = norm(R-S);

if dHorizontal == 0
    error('Source and receiver have identical XY coordinates.');
end

u = deltaXY/dHorizontal;

% Unit vector perpendicular to source-receiver direction
n = [-u(2) u(1)];

%% ------------------------------------------------------------
% DEFAULT MAXIMUM LATERAL SEARCH DISTANCE
% ------------------------------------------------------------

if isempty(opts.maxLateralOffset)

    opts.maxLateralOffset = ...
        min( ...
            0.35*dHorizontal, ...
            250);

    opts.maxLateralOffset = ...
        max( ...
            opts.maxLateralOffset, ...
            3*gridSpacing);

end

%% ------------------------------------------------------------
% INITIAL OUTPUT
% ------------------------------------------------------------

paths = struct;

paths.top = [];
paths.left = [];
paths.right = [];

paths.meta = struct;

paths.meta.dDirect = dDirect;
paths.meta.dHorizontal = dHorizontal;
paths.meta.LOS = true;

%% ============================================================
% TOP PATH
% =============================================================

nTop = ...
    max( ...
        3, ...
        ceil(dHorizontal/opts.topSpacing)+1);

t = linspace(0,1,nTop);

xLine = ...
    S(1) + t*(R(1)-S(1));

yLine = ...
    S(2) + t*(R(2)-S(2));

zTerrainLine = ...
    interp2( ...
        xTerrain, ...
        yTerrain, ...
        Zterrain, ...
        xLine, ...
        yLine, ...
        'linear');

sLine = ...
    t*dHorizontal;

%% LOS height

zLOS = ...
    S(3) + ...
    t*(R(3)-S(3));

blocked = ...
    any( ...
        zTerrainLine(2:end-1) > ...
        zLOS(2:end-1)+opts.clearanceTol);

paths.meta.LOS = ~blocked;

%% ------------------------------------------------------------
% CLEAR LINE OF SIGHT
% ------------------------------------------------------------

if ~blocked

    paths.meta.message = ...
        'Direct source-receiver line is clear.';

    return

end

%% ------------------------------------------------------------
% FIND TOP DIFFRACTION EDGES
% ------------------------------------------------------------

S2 = ...
    [0 S(3)];

R2 = ...
    [dHorizontal R(3)];

edgesTop = ...
    find_terrain_diffraction_Edges_2D( ...
        S2, ...
        R2, ...
        sLine, ...
        zTerrainLine);

if ~isempty(edgesTop)

    [z,d_ss,d_sr,e,~,L_diff] = ...
        z_multiEdge_2D( ...
            S2, ...
            R2, ...
            edgesTop);

    paths.top.z = z;
    paths.top.d_ss = d_ss;
    paths.top.d_sr = d_sr;
    paths.top.e = e;

    paths.top.L_diff = L_diff;

    paths.top.edges = edgesTop;

    paths.top.profile_s = sLine;
    paths.top.profile_z = zTerrainLine;

    paths.top.type = 'top';

end

%% ============================================================
% FAST LATERAL SEARCH
% =============================================================

offsets = ...
    linspace( ...
        opts.lateralSpacing, ...
        opts.maxLateralOffset, ...
        opts.nLateralPaths);

bestLeft = [];
bestRight = [];

bestLeftLength = inf;
bestRightLength = inf;

for sideSign = [-1 1]

    for offset = offsets

        %% ----------------------------------------------------
        % MIDPOINT OF CANDIDATE DETOUR
        % -----------------------------------------------------

        midpoint = ...
            0.5*(Sxy+Rxy);

        Mxy = ...
            midpoint + ...
            sideSign*offset*n;

        %% Check terrain boundaries

        if Mxy(1) < min(xTerrain) || ...
           Mxy(1) > max(xTerrain) || ...
           Mxy(2) < min(yTerrain) || ...
           Mxy(2) > max(yTerrain)

            continue

        end

        %% ----------------------------------------------------
        % BUILD S -> M -> R POLYLINE
        % -----------------------------------------------------

        d1 = norm(Mxy-Sxy);
        d2 = norm(Rxy-Mxy);

        n1 = ...
            max( ...
                2, ...
                ceil(d1/opts.lateralSpacing)+1);

        n2 = ...
            max( ...
                2, ...
                ceil(d2/opts.lateralSpacing)+1);

        t1 = linspace(0,1,n1);
        t2 = linspace(0,1,n2);

        x1 = ...
            S(1) + ...
            t1*(Mxy(1)-S(1));

        y1 = ...
            S(2) + ...
            t1*(Mxy(2)-S(2));

        x2 = ...
            Mxy(1) + ...
            t2*(R(1)-Mxy(1));

        y2 = ...
            Mxy(2) + ...
            t2*(R(2)-Mxy(2));

        % Avoid repeating midpoint

        xPath = ...
            [x1 x2(2:end)];

        yPath = ...
            [y1 y2(2:end)];

        %% ----------------------------------------------------
        % TERRAIN ALONG CANDIDATE PATH
        % -----------------------------------------------------

        zPath = ...
            interp2( ...
                xTerrain, ...
                yTerrain, ...
                Zterrain, ...
                xPath, ...
                yPath, ...
                'linear');

        if any(isnan(zPath))
            continue
        end

        %% ----------------------------------------------------
        % CUMULATIVE HORIZONTAL PATH DISTANCE
        % -----------------------------------------------------

        ds = ...
            sqrt( ...
                diff(xPath).^2 + ...
                diff(yPath).^2);

        sPath = ...
            [0 cumsum(ds)];

        horizontalCandidateLength = ...
            sPath(end);

        %% ----------------------------------------------------
        % MAP TO 2D PROFILE
        % -----------------------------------------------------

        S2lat = ...
            [0 S(3)];

        R2lat = ...
            [horizontalCandidateLength R(3)];

        edgesLat = ...
            find_terrain_diffraction_Edges_2D( ...
                S2lat, ...
                R2lat, ...
                sPath, ...
                zPath);

        %% ----------------------------------------------------
        % IF NO TERRAIN EDGE EXISTS
        %
        % The candidate still contains horizontal detour.
        % -----------------------------------------------------

        if isempty(edgesLat)

            candidateLength = ...
                sqrt( ...
                    horizontalCandidateLength^2 + ...
                    (R(3)-S(3))^2);

            candidate.z = ...
                candidateLength-dDirect;

            candidate.d_ss = candidateLength;
            candidate.d_sr = 0;

            candidate.e = 0;

            candidate.L_diff = candidateLength;

            candidate.edges = [];

        else

            [zLat,dssLat,dsrLat,eLat,~,LdiffLat] = ...
                z_multiEdge_2D( ...
                    S2lat, ...
                    R2lat, ...
                    edgesLat);

            % Include lateral horizontal detour in path difference.

            lateralExtra = ...
                horizontalCandidateLength - ...
                dHorizontal;

            zLat = ...
                zLat + lateralExtra;

            LdiffLat = ...
                dDirect + zLat;

            candidate.z = zLat;
            candidate.d_ss = dssLat;
            candidate.d_sr = dsrLat;
            candidate.e = eLat;

            candidate.L_diff = LdiffLat;

            candidate.edges = edgesLat;

        end

        candidate.profile_s = sPath;
        candidate.profile_z = zPath;

        candidate.x = xPath;
        candidate.y = yPath;

        candidate.offset = offset;

        %% ----------------------------------------------------
        % KEEP SHORTEST CANDIDATE
        % -----------------------------------------------------

        if sideSign > 0

            candidate.type = 'left';

            if candidate.L_diff < bestLeftLength

                bestLeft = candidate;

                bestLeftLength = ...
                    candidate.L_diff;

            end

        else

            candidate.type = 'right';

            if candidate.L_diff < bestRightLength

                bestRight = candidate;

                bestRightLength = ...
                    candidate.L_diff;

            end

        end

    end

end

%% ============================================================
% SAVE BEST LATERAL PATHS
% =============================================================

if ~isempty(bestLeft)
    paths.left = bestLeft;
end

if ~isempty(bestRight)
    paths.right = bestRight;
end

%% ============================================================
% META INFORMATION
% =============================================================

paths.meta.topFound = ...
    ~isempty(paths.top);

paths.meta.leftFound = ...
    ~isempty(paths.left);

paths.meta.rightFound = ...
    ~isempty(paths.right);

paths.meta.nLateralCandidates = ...
    opts.nLateralPaths;

paths.meta.maxLateralOffset = ...
    opts.maxLateralOffset;

end