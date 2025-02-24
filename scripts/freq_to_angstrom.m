% Function to convert frequency to angstroms
function angstrom = freq_to_angstrom(freq)
    speed_of_light = 29979245800; % Speed of light in centmeters per second
    cm = speed_of_light ./ freq;
    angstrom = cm_to_angstrom(cm);
end
