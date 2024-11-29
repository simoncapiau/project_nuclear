def molarMass(X):
    """
    Database: http://wwwndc.jaea.go.jp/NuC/
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        i.e.: for the Uranium 235, X = 'U235'; for a neutron X = 'n'
    :return: double
        molar mass of the nucleon, or nuclide X in [kg/mol]
    """

    molar_mass_data = {
        "U233": 233.039636574,
        "U235": 235.043931368,
        "U238": 238.050789466,
        "U239": 239.054294518,
        "Pu239": 239.052164844,
        "Pu240": 240.053815008,
        "Np239": 239.052940487,
        "Xe135": 134.907226844,
        "Th232": 232.038060026,
        "Th233": 233.041586541,
        "Pa233": 233.040248815,
        "n": 1.0086649
    }

    
    if X in molar_mass_data:
        M = molar_mass_data[X] * 1e-3
    else:
        print('\n ---------------- WARNING ----------------- \n No molar mass for species (', X, ')')
        M = 0

    return M

# Try your function:
#Mm = molarMass('U235')
#print(Mm)