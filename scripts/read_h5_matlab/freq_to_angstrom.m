function angstrom = freq_to_angstrom(freq)
    angstrom = cm_to_angstrom(speed_of_light() ./ freq);
end