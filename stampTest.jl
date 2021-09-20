### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ ca4ccad0-121b-11ec-2869-3b8f6afc2528
import Pkg

# ╔═╡ 7a7e21fa-4d13-4a06-8eba-5d857c726fe4
Pkg.activate(".")

# ╔═╡ 430a0642-4a75-4629-9aed-1c2552acb7af
using Revise

# ╔═╡ 15390d18-055d-4397-a73a-3f615215be38
using MOSLab

# ╔═╡ a84cbe79-86f5-483d-8af0-3ad7a8e3189f
begin

	NpolyDoping = 5e17 # Gate Npoly doping concentration in cm^-3
	N_a = 5e17 # Bulk doping concentration in cm^-3
	E_a = 0.044 # The ground state energy of the acceptor in eV
	t_ox = 4e-7 # oxide thickness in cm
	t_si = 1e-6
	MOSf(T) = MOSStructure(NPoly(NpolyDoping),SiO2(),SemiconductorData(T,BoltzmanDist(),PSilicon(N_a,E_a)),t_ox,t_si) ## Calculate Parameters of a MOS Structure the given parameters using Boltzman Distribution at temperature T 
	
end

# ╔═╡ be32370a-39ca-4c37-a3a3-1d2b263e5980
net = Netlist()

# ╔═╡ d7a825c6-1cae-44fa-97fb-7d9d9f9b3cc7
begin
	I1 = CircuitComponent("I1",CurrentSource(1))
	R1 = CircuitComponent("R1",Resistor(1))
	R2 = CircuitComponent("R2",Resistor(1))
	R3 = CircuitComponent("R3",Resistor(1))
	G = CircuitComponent("G",TransConductance(1))
end

# ╔═╡ 21d4f9bf-0625-4c03-93a6-0789081b402c
begin
	addComponent(net,R1,Dict("p"=>1,"m"=>2))
	addComponent(net,R2,Dict("p"=>2,"m"=>3))
	addComponent(net,R3,Dict("p"=>3,"m"=>0))
	#addComponent(net,I1,Dict("p"=>1,"m"=>0))
	#addComponent(net,G,Dict("ip"=>1,"im"=>2,"op"=>0,"om"=>2))
end

# ╔═╡ 3163e019-e4ba-4e81-a272-32910bc399a9
net

# ╔═╡ 0d805a03-40d1-4d7a-9ff9-c30abe120bf1
 

# ╔═╡ Cell order:
# ╠═ca4ccad0-121b-11ec-2869-3b8f6afc2528
# ╠═430a0642-4a75-4629-9aed-1c2552acb7af
# ╠═7a7e21fa-4d13-4a06-8eba-5d857c726fe4
# ╠═15390d18-055d-4397-a73a-3f615215be38
# ╠═a84cbe79-86f5-483d-8af0-3ad7a8e3189f
# ╠═be32370a-39ca-4c37-a3a3-1d2b263e5980
# ╠═d7a825c6-1cae-44fa-97fb-7d9d9f9b3cc7
# ╠═21d4f9bf-0625-4c03-93a6-0789081b402c
# ╠═3163e019-e4ba-4e81-a272-32910bc399a9
# ╠═0d805a03-40d1-4d7a-9ff9-c30abe120bf1
