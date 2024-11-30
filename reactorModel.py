import numpy as np
import sys
import matplotlib.pyplot as plt

import molarMass as mM
import halfLife as hL
import crossSection as cS


def reactorModel(fuelCompo, FPCompo, t_final, n_th_init, n_fa_init, mTot):

    """
    Computes the evolution of the species in the reactor, and compute the power evolution in time
    :param fuelCompo: class
    :param FPCompo: class
            FPCompo.Xe135: percentage of FP that are considered as poison (in this project, we assume only
                Xe135 is produced as a poison)
            FPCompo.FP: percentage of FP that are considered "others fission products" of U235 (assume only one isotope
                different than Xe135 is produced by the fission)
    :param t_final: double
            final time of the simulation in seconds [s]
    :param n_th_init: double
            initial number of thermal neutrons
    :param n_fa_init: double
            initial number of fast neutrons
    :param mTot: double
            total number of kg of fuel at the initial state in [kg]
    :return fuelBurnup: array
            depending on the fuel, returns the burnup
    """
    # units
    eV  = 1.60218e-19 
    amu = 1.66e-27
    
    # parameters
    V = 10           # [m^3] volume of the reactor
    dt = 1e-4        # [s]   time step
    E_fast = 1e6     # [ev]
    E_th = 25e-3     # [ev]
    T_fast = 5e-4    # [s] half-time for a fast neutron
    T_fpro = 1       # [s] half-time for a fission product
    nu = 2.3         # nombre de neutrons produits par une fission
    NA = 6.02214076 * 10**23     # nombre d'Avogadro [1/mol]

    PXe135 = 0.05   # proportion de Xe135 produit par la fission de l'U235 et Pu239
    PSn134 = 0.05   # proportion de Sn134 produit par la fission du Pu241


    # variables d'etat : neutrons rapides, neutrons thermiques, nucleons
    # systeme de la forme : dotx = f(x)

    # [nr, nth, Kr95, Zr104, Sn134, Xe135, ]
    
    # index des differents elements dans le vecteur
    itime = 0    # temps
    infast  = 1  # neutrons rapides
    inth = 2     # neutrons thermiques
    iKr95 = 3    # index du Kr95
    iZr104 = 4
    iSn134 = 5
    iXe135 = 6
    # iCe135 = 7
    iXe136 = 8
    iU235 = 9
    iU236 = 10
    iU237 = 11
    iNp237 = 12
    iU238 = 13
    iU239 = 14
    iNp239 = 15
    iPu239 = 16
    iPu240 = 17
    iPu241 = 18
    iAm241 = 19
    iPu242 = 20
    iAm242 = 21
    iCm242 = 22
    iPu243 = 23
    iAm243 = 24
    iCm243 = 25
    iAm244 = 26
    iCm244 = 27
    iotherFP = 28    # produit de fission quelconque
    iE = 29    # variable d'etat pour l'energie

    # les quantites d'elements sont exprimees en mol

    X = np.zeros((1000, 30))    # vecteur des variables d'etat en fonction du temps, l'index 0 donne les conditions initiales
    F = np.zeros(30)            # vecteur des derivees des variables d'etat

    vfast = 3*1e8 * np.sqrt(2*E_fast  / 939565379)  # vitesse des  neutrons rapides
    vth = 3*1e8 * np.sqrt(2*E_th  / 939565379)
    ln2 = np.ln(2)

    for j in range(1, 1000):
        F = np.zeros(26)

        nfast = X[j-1][inth]
        nth = X[j-1][infast]

        Ith = NA * nth * vth/V   # flux de neutrons thermiques en 1/m**2 s

        Ifast = NA * nfast * vfast/V

        # temps
        F[itime] = 1

        # neutrons rapides
        F[infast] = - nfast * (ln2/T_fast) \
                    + nu * cS.crossSection("U235", "Fission", ) * Ith * X[j-1][iU235]   \
                    + nu * cS.crossSection("Pu239", "Fission", ) * Ith * X[j-1][iPu239] \
                    + nu * cS.crossSection("Pu241", "Fission", ) * Ith * X[j-1][iPu241]
        # il faut ajouter des termes pour la perte, le controle et les neutrons retardes


        # neutrons thermiques
        #      diffusion
        F[inth] += nfast * (ln2/T_fast) * (E_fast - E_th)
        #      fissions
        F[inth] += - cS.crossSection("U235", "Fission", ) * Ith * X[j-1][iU235]    \
                   - cS.crossSection("Pu239", "Fission", ) * Ith * X[j-1][iPu239]  \
                   - cS.crossSection("Pu241", "Fission", ) * Ith * X[j-1][iPu241]
        #      captures
        F[inth] += - cS.crossSection("Xe135", "Capture", ) * Ith * X[j-1][iXe135]  \
                   - cS.crossSection("U235", "Capture", ) * Ith * X[j-1][iU235]    \
                   - cS.crossSection("U236", "Capture", ) * Ith * X[j-1][iU236]    \
                   - cS.crossSection("U238", "Capture", ) * Ith * X[j-1][iU238]    \
                   - cS.crossSection("Pu239", "Capture", ) * Ith * X[j-1][iPu239]  \
                   - cS.crossSection("Pu240", "Capture", ) * Ith * X[j-1][iPu240]  \
                   - cS.crossSection("Pu241", "Capture", ) * Ith * X[j-1][iPu241]  \
                   - cS.crossSection("Am241", "Capture", ) * Ith * X[j-1][iAm241]  \
                   - cS.crossSection("Pu242", "Capture", ) * Ith * X[j-1][iPu242]  \
                   - cS.crossSection("Am242", "Capture", ) * Ith * X[j-1][iAm242]  \
                   - cS.crossSection("Cm242", "Capture", ) * Ith * X[j-1][iCm242]  \
                   - cS.crossSection("Am243", "Capture", ) * Ith * X[j-1][iAm243]  \
                   - cS.crossSection("Cm243", "Capture", ) * Ith * X[j-1][iCm243]  
        
        # il faut ajouter des termes pour la perte, le controle et les neutrons retardes

        # Kr95
        F[iKr95] = - X[j-1][iKr95] * (ln2/hL.halfLife("Kr95", "betaMinus")) \
                    + PXe135 * cS.crossSection("U235", "Fission", ) * Ith * X[j-1][iU235]   \
                    + PXe135 * cS.crossSection("Pu239", "Fission", ) * Ith * X[j-1][iPu239]

        # Zr104
        F[iZr104] = - X[j-1][iZr104] * (ln2/hL.halfLife("Zr104", "betaMinus")) \
                    + PSn134 * cS.crossSection("Pu241", "Fission", ) * Ith * X[j-1][iPu241]

        # Sn134
        F[iSn134] = - X[j-1][iXe135] * (ln2/hL.halfLife("Sn134", "betaMinus")) \
                    - cS.crossSection("Sn134", "Capture", ) * Ith * X[j-1][iXe135]  \
                    + PSn134 * cS.crossSection("Pu241", "Fission", ) * Ith * X[j-1][iPu241]
        
        # Xe135
        F[iXe135] = - X[j-1][iXe135] * (ln2/hL.halfLife("Xe135", "betaMinus")) \
                    - cS.crossSection("Xe135", "Capture", ) * Ith * X[j-1][iXe135]  \
                    + PXe135 * cS.crossSection("U235", "Fission", ) * Ith * X[j-1][iU235]   \
                    + PXe135 * cS.crossSection("Pu239", "Fission", ) * Ith * X[j-1][iPu239]

        # Xe136
        F[iXe136] = + cS.crossSection("Xe135", "Capture", ) * Ith * X[j-1][iXe135]

        # U235
        F[iU235] = - cS.crossSection("U235", "Fission", ) * Ith * X[j-1][iU235]   \
                   - cS.crossSection("U235", "Capture", ) * Ith * X[j-1][iU235]

        # U236
        F[iU236] = + cS.crossSection("U235", "Capture", ) * Ith * X[j-1][iU235]   \
                   - cS.crossSection("U236", "Capture", ) * Ith * X[j-1][iU236]
        
        # U237
        F[iU237] = + cS.crossSection("U236", "Capture", ) * Ith * X[j-1][iU236]   \
                   - X[j-1][iU237] * (ln2/hL.halfLife("U237", "betaMinus"))
        
        # Np237
        F[iNp237] = + X[j-1][iU237] * (ln2/hL.halfLife("U237", "betaMinus"))

        # U238
        F[iU238] = - cS.crossSection("U238", "Fission", ) * Ith * X[j-1][iU238]   \
                   - cS.crossSection("U238", "Capture", ) * Ith * X[j-1][iU238]
        
        # U239
        F[iU239] = + cS.crossSection("U238", "Capture", ) * Ith * X[j-1][iU238]   \
                   - X[j-1][iU239] * (ln2/hL.halfLife("U239", "betaMinus"))
        
        # Np239
        F[iNp239] = + X[j-1][iU239] * (ln2/hL.halfLife("U239", "betaMinus")) \
                    - X[j-1][iNp239] * (ln2/hL.halfLife("Np239", "betaMinus"))
        
        # Pu239
        F[iPu239] = + X[j-1][iNp239] * (ln2/hL.halfLife("Np239", "betaMinus"))     \
                    - cS.crossSection("Pu239", "Fission", ) * Ith * X[j-1][iPu239]   \
                    - cS.crossSection("Pu239", "Capture", ) * Ith * X[j-1][iPu239]
        
        # Pu240
        F[iPu240] = + cS.crossSection("Pu239", "Capture", ) * Ith * X[j-1][iPu239]   \
                    - cS.crossSection("Pu240", "Capture", ) * Ith * X[j-1][iPu240]
        
        # Pu241
        F[iPu241] = + cS.crossSection("Pu240", "Capture", ) * Ith * X[j-1][iPu240]   \
                    - cS.crossSection("Pu241", "Fission", ) * Ith * X[j-1][iPu241]   \
                    - cS.crossSection("Pu241", "Capture", ) * Ith * X[j-1][iPu241]
        
        # Pu242
        F[iPu242] = + cS.crossSection("Pu241", "Capture", ) * Ith * X[j-1][iPu241]

        # autres PF
        F[iotherFP] = + (1-PXe135) * cS.crossSection("U235", "Fission", ) * Ith * X[j-1][iU235]   \
                      + (1-PXe135) * cS.crossSection("Pu239", "Fission", ) * Ith * X[j-1][iPu239] \
                      + (1-PSn134) * cS.crossSection("Pu241", "Fission", ) * Ith * X[j-1][iPu241] \
                      + X[j-1][iXe135] * (ln2/hL.halfLife("Xe135", "betaMinus"))   \
                      + X[j-1][iKr95] * (ln2/hL.halfLife("Kr95", "betaMinus"))   \
                      + X[j-1][iZr104] * (ln2/hL.halfLife("Zr104", "betaMinus"))   \
                      + X[j-1][iSn134] * (ln2/hL.halfLife("Sn134", "betaMinus"))
        



        ######################
        # methode d'Euler explicite
        X[j] = X[j-1] + F * dt

    return fuelBurnup

class Fuel:
    def __init__(self):
        self.U235 = 3
        self.U238 = 97
        self.Pu239 = 0
        self.Th232 = 0

class FP:
    def __init__(self):
        self.Xe135 = 5
        self.FP = 95

fuelCompo = Fuel()
FPCompo = FP()



# Try the function

# Mm = mM.molarMass('U235')
# hl = hL.halfLife(X='Pa233', Transfo='BetaMinus')
# cs = cS.crossSection(X='Pu240', Transfo='Fission', E_neutron=np.logspace(-5,6,10000).tolist())
#
# print(Mm, hl, cs)

# fBurnup = reactorModel(fuelCompo=fuelCompo, FPCompo=FPCompo, t_final=200., n_th_init=1e+10, n_fa_init=0., mTot=25.)
