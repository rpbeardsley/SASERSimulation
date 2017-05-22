%calculate the inter-well form factor Jinter

%variables
%TWWF - the two well model wave function for a coupled well
%q - phonon wave vector
%theta - small theta means small qpar
%d - SL period

function InterFormFac = JinterFunc(TWWF,q,theta,d,Jintra,CoupleFactor)

TWWFinter = TWWF;

qz = q*cos(theta);
qpar = q*sin(theta);

%put the TWWF in neighbouring wells and on the same z axis
TWWFinternpls1 = TWWFinter;
TWWFinternpls1(:,1) = TWWFinternpls1(:,1) + d;

zstep = abs(TWWFinter(1,1) - TWWFinter(2,1));
resize = round(d/zstep);

TWWFinter(1:(resize+1),:) = [];
TWWFinternpls1((size(TWWFinternpls1,1)-(resize-1)):end,:) = [];
TWWFinternpls1(end,:) = []; %It puts a point way out on the end for some reason. I delete it with this

zax = TWWFinter(:,1);
TWWFinter(:,1) = [];
TWWFinternpls1(:,1) = [];

%figure
%hold on
%plot(zax,TWWFinter,'-');
%plot(zax,TWWFinternpls1,'rx');

%%%%%%%%%%%%%%%
%alternative method (get Jinter from Jintra)
%AB1 = 4*(CoupleFactor^2);
%AB2 = (sin( (q*d)/2 ))^2;
%AB = AB1*AB2*Jintra;
%%%%%%%%%%%%%%%

%the function to be integrated
PhonDisp = exp(i.*qz.*zax);
Intfun = TWWFinter.*TWWFinternpls1.*PhonDisp;

%integrate and multiply by complex conjugate
Ingrtd = trapz(zax,Intfun);
Ingrtd = Ingrtd.*conj(Ingrtd);

InterFormFac = Ingrtd;



