%Calculate the screening factor kappa

%variables
%e -elementary charge
%ne -electron density
%q -SL wavevector
%theta -angle between q and qz
%d -SL period
%RelPerm -GaAs relative permittivity
%AbsPerm -Absolute permittivity freespace
%me -electron mass
%mw -effective mass in well
%a -well width


function ScreFacKappa = SreeningFactor(e,ne,q,theta,d,RelPerm,AbsPerm,me,mw,a,omega);

mstar = me*mw;
qparallel = q*sin(theta);

F1 = (e^2)*ne*(qparallel^2)*d / (8*RelPerm*AbsPerm*mstar);

F2 = ( 4 / (1-cos(q*d)) ) - ( (a/d) * ( (4/3) - (5/pi^2) ) );

omegaPLSquared = F1 * F2;

ScreFacKappa = 1 - omegaPLSquared / omega^2;