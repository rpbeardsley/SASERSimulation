%electron dispersion using the Kronig-Penney model

%constants
me=5.68572e-32; %kg electron mass
hbar=1.054e-34; %hbar
db=5.9e-9 ; %m barrier thickness
dw=3.9e-9 ; %m well thickness
a=db+dw; %m SL period
Emax=1200e-3*1.6e-19; %J maximum energy for calculation
Egw = 1.42*1.6e-19; %J -energy gap of the well at 10K
Egb = 2.16*1.6e-19; %J -energy gap of the barrier at 10K
V0 = Egb - Egw; %J barrier height

%create the energy range over which to calculate
%E = 0.1e-24:0.1e-24:Emax;
E = 0.1e-25:0.1e-25:4e-20;

%calculate the wavevectors in the barrier and well for each E
D = V0 - E(:,:);
Kb = ((D.*2*me)/hbar^2).^0.5;
Kw = (E.*2*me/hbar^2).^0.5;

%calculate the dispersion
Q = ((Kb.^2)-(Kw.^2))./(Kb.*Kw.*2) ;
W = sinh(Kb.*db).*sin(Kw.*dw);
T = cosh(Kb.*db).*cos(Kw.*dw);
U = (Q.*W + T);


I_band = find(U <=1 & U >= -1);
I_gap = find(U > 1 | U < -1);
U_band = U(I_band);
E_band = E(I_band);
U_gap = U(I_gap);
E_gap = E(I_gap);

ka_band = acos(U_band);

%change units for plotting
k_band = ka_band/pi;
E_band = E_band./1.6e-22;


%plot the result
figure
plot(k_band,E_band,'x')
xlabel('k (\pi/a)');
ylabel('energy (meV)');



