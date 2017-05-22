%calculate the intra-well form factor Jintra

%variables
%TWWF - the two well model wave function for a coupled well
%q - phonon wave vector
%theta - small theta means small qpar

function IntraFormFac = JintraFunc(TWWF,q,theta)

qz = q*cos(theta);
qpar = q*sin(theta);

TWWFintra = TWWF;

%seperate axes
zax = TWWFintra(:,1);
TWWFintra(:,1) = [];

%integrate and multiply by complex conjugate
Intfun = (TWWFintra.^2) .* exp(i.*qz.*zax);
Ingrtd = trapz(zax,Intfun);
Ingrtd = Ingrtd.*conj(Ingrtd);

IntraFormFac = Ingrtd;
