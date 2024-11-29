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

    V = 10           # [m^3] volume of the reactor
    dt = 1e-4        # [s]   time step
    E_fast = 1e6     # [ev]
    E_th = 25e-3     # [ev]
    T_fast = 5e-4    # [s] half-time for a fast neutron
    T_fpro = 1       # [s] half-time for a fission product

    # variables d'etat : neutrons rapides, neutrons thermiques, nucleons
    # systeme de la forme : dotx = f(x)

    # [nr, nth, Kr95, Zr104, Sn134, Xe135, ]
    
    # index des differents elements dans la matrice
    inth = 0     # neutrons thermiques
    infast  = 1  # n0 rapides
    iKr95 = 2   # index du Kr85
    iZr104 = 3
    iSn134 = 4
    iXe135 = 5
    iCe135 = 6
    iXe136 = 7
    iU235 = 8
    iU236 = 9
    iU237 = 10
    iNp237 = 11
    iU238 = 12
    iU239 = 13
    iNp239 = 14
    iPu239 = 15
    iPu240 = 16
    iPu241 = 17
    iAm241 = 18
    iAm242 = 19
    iCm242 = 20
    iPu243 = 21
    iAm243 = 22
    iCm243 = 23
    iAm244 = 24
    iCm244 = 25
    iE = 26    # variable d'etat pour l'energie

    X = np.zeros((26, 1000))    # vecteur des variables d'etat au depart
    F = np.zeros(26)            # vecteur des derivees des variables d'etat

    vfast = 3*1e8 * np.sqrt(2*E_fast  / 939565379)
    vth = 3*1e8 * np.sqrt(2*E_th  / 939565379)
    ln2 = np.ln(2)

    for j in range(1, 1000):
        F = np.zeros(26)

        nfast = X[j-1][inth]
        nth = X[j-1][infast]

        Ith = nth * vth/V
        Ifast = nfast * vfast/V

        # neutrons rapides
        F[infast] = - nfast * (T_fast)/ln2    # il faut ajouter un terme pour la perte et le controle

        F[iE] += nfast * (T_fast/ln2) * (E_fast - E_th)

        # neutrons thermiques
        F[inth] = + nfast * (T_fast/ln2)



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
