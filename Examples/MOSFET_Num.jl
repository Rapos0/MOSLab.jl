using MOSLab
using Plots
using LaTeXStrings
Vg = 0.1 # Gate Voltage
Vd = 1.0 # Dain Voltage
Nx = 101 # Number of grid points on the x direction
Ny = 101 # Number of grid points on the y direction
N_d = 6e16 # Dopants density at the drain and source side  in cm ^-3
T = 300.0 # Temperature in Kelvin
Nb = 5.7e17 # Dopants ensity at the bulk in cm ^-3
l = .180e-4 # Channel length
h = 0.5e-4 # Silicon depth
t_ox = 4e-7 # Oxide thickens in cm
Ms = MOSFETInputDeck(N_d,Nb,T,l,h,h/2.0,l,t_ox,:NMOS) # Create the input deck for an NMOS transistor
Vs = ψs_PSP(Vg,Ms.Gate) # Calculate the surface potential at the gate
G,psi,nn,pp = MOSFETSimulation(Ms,Vs,Vd,Nx,Ny;verbose=false,ξ₀=0.01) # run 2D MOSFET simulation
gr()
X = unique(G[1,:])
Y = unique(G[2,:])
contour(X,Y,psi,fill=true)
surface(G[1,:],G[2,:],psi,fill=true)
xlabel!(L"x \ [\mu m]")
ylabel!(L"y \ [\mu m]")
plot!(zlabel = L"\psi \ [V]")
savefig("Psi_Surface.png")
yvals = unique(G[2,:])
ySel = G[2,:] .== yvals[end-1]
sum(ySel)
minimum(psi[ySel])
plot(G[1,ySel],psi[ySel],label="\$ E_C (x,$(yvals[end-1]) ) \$",legend=:bottomleft)
savefig("Psi_channel.png")
plot(G[1,ySel],nn[ySel],label="\$ E_C (x,$(yvals[end-1]) ) \$",legend=:topright)
savefig("nn_channel.png")