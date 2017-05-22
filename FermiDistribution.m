%Calculate the Fermi distribution function

%function FermiDistribution = FermiDistribution()

%Define some constants
hbar = 1.054571628e-34; %Js
Te = 10; %K electron temperature
me = 9.1094e-31; %kg -electron mass
mw = 0.063; %kg -effective electron mass in well
mb = 0.15; %kg -effective electron mass in the barrier
ne = 2e20; %m^(-1) electron denisty in the wells

%Calculate the Fermi level in the well
R = (ne*pi*hbar^2) / (me*mw*Te);
Ef = Te * (log(exp(R) - 1));

%create the energy range over which to work
E = [1e-24:1e-20:2.064e-19]';

%calculate the Fermi distribution
U = (E - Ef) / Te;
f = 1 / (1 + exp(U));

%f = FermiDistribution

%plot the distribution
figure;
hold on
plot(E,U)
%plot(E,f);





