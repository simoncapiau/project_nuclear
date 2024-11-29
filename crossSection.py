def crossSection(X, Transfo, E_neutron):
    """
    Computes the cross section for a specific energy level
    ENDF database: https://www-nds.iaea.org/exfor/endf.htm

    ----------------
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        e.g.: for the Uranium 235, X = 'U235'
    :param Transfo: string
        name of the considered transformation, should be one of the following in the list
        ['Fission', 'Capture']
    :param E_neutron: array, list
        value(s) of the cross section for an incident neutron of energy level equal to E_neutron in [eV]
        the value should be in the range [1e-5; 2e7] [eV]
    :return: array, list
        values of the cross section in [barn] corresponding to the input parameters
    """

    # Your work : here
    CrossSections = {
        "Fission": {
            "U235" : 590.66974,
            "Pu239": 751.4277, 
            "Pu241": 1018.6381,
        },
        "Capture": {
            "Kr95" : 0, # NOPE
            "Zr104": 0, # NOPE
            "Sn134": 0, # NOPE
            "Xe135": 2666632.5,
            "Ce135": 35.103786,
            "Xe136": 0, # NOPE
            "U236" : 5.153439,
            "U237" : 455.4131,
            "Np237": 176.58429,
            "U238" : 2.699, # NOPE
            "U239" : 22.66,
            "Np239": 45.29,
            "Pu240": 287,
            "Am241": 689.123,
            "Am242": 255.7,
            "Cm242": 19.1517,
            "Pu243": 99.24,
            "Am243": 0, # NOPE
            "Cm243": 131.531,
            "Am244": 0, # NOPE
            "Cm244": 15.248
        },
    }


    # ==================================  Check arguments  ==================================
    #  Check that the nucleon/nuclide asked, and that the associated transformation exists in the database
    if X != 'Th232' and X != 'Th233' and X != 'Pa233' and X != 'U233' and X != 'U235' and X != 'U236' and X != 'U237' and X != 'U238' \
            and X != 'U239' and X != 'Np239' and X != 'Pu239' and X != 'Pu240' and X != 'Xe135':
        print('\n WARNING : There is no database for element ', X, '. \n Please check function information')
    
    # sigma = 0.

    return sigma

# Try the function:
# cs = crossSection(X='Pu240', Transfo='Fission', E_neutron=np.logspace(-5,6,10000).tolist())