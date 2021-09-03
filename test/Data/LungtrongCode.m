% Surface potential vs. Gate voltage for MOS Capacitor 
% Date: Oct. 23, 2012 % Author: Xingshu Sun and Mark S. Lundstrom(Purdue University) % 
% References % [1] Mark Lundstrom and Xingshu Sun (2012), "Notes on the Solution of the Poisson-Boltzmann Equation for MOS Capacitors and MOSFETs, 2nd Ed." http://nanohub.org/resources/5338. 
%Initialize the range of the surface potential 
psi = -0.2:0.01:1.2; %[V] 
%Specify the physical constants 
epsilon_siO2 = 4*8.854*1e-14; %Permittivity of SiO2 [F/cm]
epsilon_si = 11.68*8.854*1e-14; %Permittivity of Si [F/cm] 
k_b = 1.380e-23; %Boltzmann constant [J/K] 
q = 1.6e-19; %Elementary charge [C]

TL = 300; %Room temperature[K] 
ni = 6.321247991453874e9; %Intrinsic semiconductor carrier density [cm-3] 
N_A = 5e17; %Acceptor concentration [cm-3]
tox = 2e-7; %The thickness of SiO2 [cm]
V_FB = 0; %The flat-band voltage [V]

%Calculate the capacitance 
cox = epsilon_siO2/tox; %[F/cm2] 
%Calculate the F function (to calculate the total charge) with respect to different surface potentials.  See: eqn. (10) in [1] 
F_psi1 = exp(-psi*q/k_b/TL) + psi*q/k_b/TL -1; 
F_psi2 = ni^2/N_A^2*(exp(psi*q/k_b/TL) - psi*q/k_b/TL -1);
F_psi = sqrt(F_psi1 + F_psi2); 
%Initialize the vectors 
psi_length = length(psi); 
Qs = zeros(psi_length,1);  
V_G = zeros(psi_length,1); 
psi_B = zeros(psi_length,1); 

for i = 1:psi_length 
    if psi(i) <= 0 
        Qs(i) = sqrt(2*epsilon_si*k_b*TL*N_A)*F_psi(i); 
        %Calculate charge when the surface potential is negative
    else
        Qs(i) = -sqrt(2*epsilon_si*k_b*TL*N_A)*F_psi(i);
    %Calculate charge when the surface potential is positive    
    end  
    psi_B(i) = k_b*TL/q*log(N_A/ni); 
    V_G(i) = V_FB+psi(i)-Qs(i)/cox; 
end

%Plot total charge vs. surface potential (psi_S)
figure(1)
semilogy(psi,abs(Qs)/q,'r','linewidth',3); 
hold on 
set(gca, 'xlim', [-0.4 1.4], 'ylim', [1e11 1e14]); 
set(gca,'fontsize',13); xlabel('\psi_S (V)'); 
ylabel('|Qs|/q (cm^-2)');