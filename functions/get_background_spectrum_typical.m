function [f,Lp_bg,LpA_bg,description] = ...
    get_background_spectrum_typical(environmentType)

    f = [63 125 250 500 1000 2000 4000 8000];

    switch lower(environmentType)

        case 'rural'

            Lp_bg = ...
                [42 37 32 27 22 18 14 12];

            description = ...
                'Typical quiet rural outdoor background';

        case 'suburban'

            Lp_bg = ...
                [47 42 37 32 27 23 19 17];

            description = ...
                'Typical quiet suburban outdoor background';

        case 'urban'

            Lp_bg = ...
                [52 47 42 37 32 28 24 22];

            description = ...
                'Typical relatively quiet urban outdoor background';

        case 'commercial'

            Lp_bg = ...
                [57 52 47 42 37 33 29 27];

            description = ...
                'Typical commercial outdoor background';

        otherwise

            error([ ...
                'environmentType must be: ', ...
                '''rural'', ''suburban'', ', ...
                '''urban'', or ''commercial''.' ...
            ]);

    end

    A_weight = ...
        [-26.2 -16.1 -8.6 -3.2 0 1.2 1.0 -1.1];

    LpA_bg = ...
        10*log10(sum( ...
        10.^((Lp_bg + A_weight)/10)));

end