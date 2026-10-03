% For the barrier attenuation function.

% Inputs format: 
%   - S = [x_s, z_s];
%   - R = [x_r, z_r];
%   - edges = [
%              x1, z1
%              x2, z2
%              x3, z3
%                    ]; direction from SOURCE to RECIEVER

function [z, d_ss, d_sr, e, d, L_diff] = z_multiEdge_2D(S, R, edges)

    d = norm(R - S); % Straight distance

    if isempty(edges) % No edges, no barrrier
        d_ss = 0;
        d_sr = 0;
        e = 0;
        L_diff = d;
        z = 0;
        return
    end

    d_ss = norm(edges(1,:) - S); % Distance from source to the 1st edge
    d_sr = norm(R - edges(end,:)); % Distance from last edge to reciever

    if size(edges,1) == 1 % If there is only one edge
        e = 0;
    else
        edgeSegments = diff(edges,1,1);
        edgeLengths = sqrt(sum(edgeSegments.^2,2));
        e = sum(edgeLengths);
    end

    L_diff = d_ss + e + d_sr; % By ISO definition
    z = L_diff - d; % Difference between straight length to the bending

end