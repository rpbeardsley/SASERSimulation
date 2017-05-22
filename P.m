%calcaulate P(epsilon) for the calculation of the incrament

%variables
%E1 -deformation potential constant
%rho1 -density
%d -SL period
%vs -sound velocity
%mw -effective mass well
%mb -effective mass barrier
%Ef -Fermi energy
%Te -electron temperature
%En -energy
%omega -phonon frequency (angular)
%q -phonon wavevector
%theta -angle between q and qz
%d -SL period

function P = P(epsilon,Te,Ef,omega,rhow,vsw,q,theta,E1,mw,mb,me,hbar,d)

qpar = q*sin(theta);

%calculate Estar
Estar1 = me*mw / (2*(hbar^2)*(qpar^2));
Estar2 = (epsilon - ( ( (hbar^2)*(qpar^2) ) / (2*mw*me) ) )^2;
Estar = Estar1*Estar2;

%create the range of E over which to integrate
Eupper = Estar*50; %upper integral limit
step = (Eupper-Estar)/10000;
En = [Estar:step:Eupper];

%calculate the number preceding the integral
Fac1 = E1^2*omega;
Fac2 = sqrt(2)*pi*rhow*d*(vsw^2)*qpar;
Fac3 = (me*(0.9*mw + 0.1*mb) /(hbar^2))^(3/2); %effective mass should be slightly into barrier (0.9*mw + 0.1*mb)
Fac = (Fac1/Fac2)*Fac3;

%inside the integral
In1 = 1 + exp((En - Ef) / Te);
In2 = 1 + exp((Ef - En - epsilon) / Te);
In3 = sqrt((En+step) - Estar);
In = 1 ./ (In1.*In2.*In3);

%figure;
%plot (En,In,'x');

%integrate the In with respect to E between infinity and Estar
PInt = trapz(En,In); 
P = PInt * Fac;





