%Find the single quantum well wave function by solving for the boundary
%conditions numerically to obtain the energy and put this energy into the
%wave function.

%define some constants
hbar = 1.054571628e-34; %Js
me = 9.1094e-31; %kg -electron mass
mw = 0.063; % rel effective electron mass in well
mb = 0.15; % rel effective electron mass in the barrier
Egw = 1.52*1.6e-19; %J -energy gap of the well at 10K
Egb = 2.81*1.6e-19; %J -energy gap of the barrier at 10K
vs = 4730; %m/s -sound velocity in well
a = 5.9e-9; %m -well width
V0 = Egb - Egw; %J the barrier height (find and apply the band offset)

%create a vector to hold the energy values
E = [2.4e-20:1e-25:2.6e-19]';

%method computing substitution
%k = ((2*me*mw.*E).^0.5) / hbar;
%A = ((mw/mb)*((2*me*mw*V0) ./ (hbar^2 .* k.^2)) - 1).^0.5;
%B = tan( (a/2) .* k );
%Res = A - B;

%method with substituted eqs
%A = ( (mw/mb)*((V0./E)-1) ).^0.5;
%B = tan( ( (E.*(a^2)*me*mw) ./ (2*(hbar^2)) ) );
%Res = A-B; %should be 0 for the energy I am looking for

%method straight from book
k = ((2*me*mw.*E).^0.5) / hbar;
theta = (k.*a)./2;
theta0sqrd = (me*mw*V0*a^2) / 2*hbar^2;
A = tan(theta);
B = sqrt( (mw/mb) .* (theta0sqrd./theta.^2) - 1 );
Res = A-B;

figure
hold on
%plot(E,A,'b')
%plot(E,B,'g')
plot(E,Res,'r')


%b = 1
%while (a>0)
%    a = Res(b,1)
%    b = b + 1
%end
%c = ( E(b,1) + E(b-1,1) ) / 2
    
    
    