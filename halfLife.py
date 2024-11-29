def halfLife(X, Transfo):
    """
    Return the half life of a nuclide or nucleon for a specific transformation
    Data available on http://wwwndc.jaea.go.jp/NuC/ (more precisely: http://wwwndc.jaea.go.jp/CN14/index.html)

    ----------------
    :param X: string
        nuclide or nucleon (for a neutron), that follows the atomic notation of the element,
        i.e.: for the Uranium 235, X = 'U235'
    :param Transfo: string
        name of the considered transformation, should be one of the following in the list
        ['Alpha', 'BetaMinus', 'BetaPlus', 'Gamma']
    :return: the half life of X for the transformation Transfo in seconds
    """

    # Time units in seconds
    M = 60         # [s]
    H = 3600       # [s]
    D = 86400      # [s]
    Y = 31536000   # [s]

    
    half_life_data = {
        "Alpha": {
            "Th232": 14.05e9 * Y,
            "U233": 159.2e3 * Y,
            "U235": 703.8e6 * Y,
            "U236": 2.342e7 * Y,
            "U238": 4.468e9 * Y,
            "Pu238": 87.7 * Y,
            "Pu239": 2.411e4 * Y,
            "Pu240": 16.561e3 * Y,
            "Pu241": 14.29 * Y,
            "Np237": 2.144e6 * Y,
        },
        "BetaMinus": {
            "Th233": 22.3 * M,
            "Pa233": 26.975 * D,
            "U237": 6.75 * D,
            "U239": 23.45 * M,
            "U240": 14.1 * H,
            "Np238": 2.117 * D,
            "Np239": 2.356 * D,
            "Np240": 61.9 * H,
            "Xe135": 9.14 * H,
            "I135": 6.57 * H,
        },
    }

    # Validate transformation type
    if Transfo not in half_life_data:
        raise ValueError(f" The transformation required isn't matching with a nucleon. ")

    # Validate nuclide
    if X not in half_life_data[Transfo]:
        raise ValueError(f"Nuclide '{X}' is not supported for transformation '{Transfo}'. ")

    
    hl =  half_life_data[Transfo][X]

    # ==================================  Check arguments  ==================================
    #  Check that the nucleon/nuclide asked, and that the associated transformation exists in the database
    if X != 'Th232' and X != 'Th233' and X != 'Pa233' and X != 'U233' and X != 'U235'and X != 'U236'and X != 'U237' and X != 'U238'\
            and X != 'U239' and X != 'Np239' and X != 'Pu239' and X != 'Pu240' and X != 'Xe135':
        print('\n WARNING : There is no database for element ', X, '. \n Please check function information')

     # Check whether the transformation exists or not
    if Transfo != 'Alpha' and Transfo != 'BetaMinus' and Transfo != 'BetaPlus' and Transfo != 'Gamma': # Transfos Gamma très courtes
        print('\n WARNING : These transformation are not implemented :', Transfo, '.\n Please check function information')

    # hl = 0.

    return hl

# Try your function
#print(halfLife(X='Pa233', Transfo='BetaMinus'))
