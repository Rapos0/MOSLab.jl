using MOSLab
using Plots

N_d = 6e16
tox = 4e-7
T = 300.0
Nb = 5.7e17
l = 1e-4
h = 1e-4
MI = MOSFETInputDeck(N_d,Nb,T,l,h,0.5e-4,0.5e-4,tox,:NMOS)
G = Grid(MI,101,101)