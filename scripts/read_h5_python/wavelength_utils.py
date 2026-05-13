def speed_of_light():
    return 2.99792458e10


def cm_to_angstrom(cm):
    return cm * 1e8


def angstrom_to_cm(cm):
    return cm * 1e-8


def freq_to_angstrom(freq):
    return cm_to_angstrom(speed_of_light() / freq)


def vacuum_to_air(wave, to_air_limit=200.0):
    """
    https://github.com/ITA-Solar/rh/blob/master/idl/vacuumtoair.pro
    """
    wave2 = wave * wave

    fact = 1.0 + 2.735182e-4 + (1.314182e0 + 2.76249e+4 / wave2) / wave2
    fact = fact * (wave > to_air_limit) + 1.0 * (wave < to_air_limit)

    return wave / fact