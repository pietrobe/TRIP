function wave_air = vacuum_to_air(wave, to_air_limit)
    if nargin < 2
        to_air_limit = 200.0;
    end
    
    wave2 = wave .^ 2;
    
    fact = 1.0 + 2.735182e-4 + (1.314182e0 + 2.76249e+4 ./ wave2) ./ wave2;
    fact = fact .* (wave > to_air_limit) + 1.0 .* (wave < to_air_limit);
    
    wave_air = wave ./ fact;
end