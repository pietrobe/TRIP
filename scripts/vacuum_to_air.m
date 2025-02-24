function wave_air = vacuum_to_air(wave, to_air_limit)
    % Default value for to_air_limit if not provided
    if nargin < 2
        to_air_limit = 200.0;
    end
    
    wave2 = (wave .* wave) / 100.0;
    
    fact = 1.0 + 2.735182e-4 + (1.314182e0 + (1.0 / 2.76249e+4) * wave2) ./ wave2;
    fact = fact .* (wave > to_air_limit) + 1.0 * (wave < to_air_limit);
    
    wave_air = wave ./ fact;
end

