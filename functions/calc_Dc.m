% This function is where the hemi-spherical geometry is taken into account
% D_source = 0 if the source relese sound equally for all directions.

% This value is ADDED to the Lp, not to the attenuations.

function D_c = calc_Dc(D_source, nReflectingSurfaces)

    switch nReflectingSurfaces

        case 0 % free field
            D_surface = 0;
        case 1 % surface, wall
            D_surface = 3;
        case 2 % edge, two intersecting walls
            D_surface = 6;
        case 3 % corner, three intersecting walls
            D_surface = 9;
        otherwise
            error('nReflectingSurfaces must be 0, 1, 2, or 3.');
    end

    D_c = D_source + D_surface;

end