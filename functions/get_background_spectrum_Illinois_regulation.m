function [f,Lp_bg,LpA_bg,description] = ...
    get_background_spectrum_Illinois_regulation(category,timeOfDay)

    f = [63 125 250 500 1000 2000 4000 8000];

    day = [
        71 72 70 67 63 57 53 48
        64 64 63 59 54 48 42 36
        57 57 55 51 45 38 30 24
        51 51 49 45 39 33 25 19
        45 45 43 39 33 26 20 13
    ];

    night = [
        61 62 60 57 53 47 43 38
        56 56 55 51 46 40 34 28
        52 52 50 46 40 33 25 19
        46 46 44 40 34 28 20 14
        40 40 38 34 28 21 15 8
    ];

    descriptions = {
        'Noisy commercial and industrial'
        'Moderate commercial/industrial or noisy residential'
        'Quiet commercial/industrial or moderate residential'
        'Quiet residential'
        'Very quiet sparse suburban or rural'
    };

    if category < 1 || category > 5
        error('category must be an integer from 1 to 5.');
    end

    switch lower(timeOfDay)

        case 'day'
            Lp_bg = day(category,:);

        case 'night'
            Lp_bg = night(category,:);

        otherwise
            error('timeOfDay must be ''day'' or ''night''.');

    end

    A_weight = [-26.2 -16.1 -8.6 -3.2 0 1.2 1.0 -1.1];

    LpA_bg = 10*log10( ...
        sum(10.^((Lp_bg + A_weight)/10)));

    description = descriptions{category};

end