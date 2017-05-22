%calculate the increment for the verticle SASER

%Define some constants for the calculation
global hbar
global me
global mw
global mb
global Egw
global Egb
global a
global b
global d
global V0
global E
global e
global Delta
global q
global qs
global theta
global omega
global Te
global En
global E1
global rhow
global rhob
global vsw
global vsb
global ne
global Ef

%constants
hbar = 1.054571628e-34; %Js
E1 = 9*1.6e-19; %J -deformation potential constant (9eV)
AbsPerm = 8.85e-12; %Fm^-1 permittivity of freespace
RelPerm  = 12.9; %GaAs relative permittivity
e = 1.602176487e-19; %-elementary charge
me = 9.1094e-31; %kg -electron mass
kb = 1.3806503e-23; %m2 kg s^-2 K^-1 -Boltzmann constant
%electrons
mw = 0.067; % -effective electron mass in well
mb = 0.146; % -effective electron mass in the barrier
Te = 0.0007295*1.6e-19; %J -electron temperature in energy units (15K = 0.001295meV)
ne = 118; %m^-2 %e10; %m^(-3) -electron density 2e10
%phonons
%theta = 0.002; %rad -angle between q and qz
%freq = 620e9; %Hz -phonon frerquency
%omega = freq*2*pi; % -angular phonon frequency
Th = 15; %K -heater temperature
%geometry and mechanical
vsw = 4730; %m^(-1) -sound velocity in the well
vsb = 5650; %m^(-1) -sound velocity in the barrier
rhow = 5320; %kg/m^3 -density of well material
rhob = 3760; %kg/m^3 -density of the barrier material
a = 5.9e-9; %m -well width
b = 4.9e-9; %m -barrier width
d = a + b; %m -SL -period
n = 50; %number of SL periods
%energy levels
Egw = 1.52*1.6e-19; %J -energy gap of the well at 10K
Egb = 2.81*1.6e-19; %J -energy gap of the barrier at 10K
V0 = Egb - Egw; %J -the barrier height
E = 2.747195e-20; %J -energy of the first confined state (from bottom of well)
%Delta = 3e-3*1.6e-19; %J -Stark splitting

%SL sound velocity (weighted by time spent in a given layer)
tw = a / vsw;
tb = b / vsb;
ttot = tw + tb;
fractw = tw / ttot;
fractb = tb / ttot;
vsSL = (fractw*vsw) + (fractb*vsb);

%inner loop parameter (Delta in this case)
lower2 = 2.55e-3*1.6e-19; %2.55meV opt
upper2 = 4.5e-3*1.6e-19; %4.5meV opt
step2 = 0.01e-3*1.6e-19; %0.01 opt
freepara2 = [lower2:step2:upper2];
IncrVls2 = zeros(size(freepara2,1),1);
JintraVls2 = zeros(size(freepara2,1),1);
JinterVls2 = zeros(size(freepara2,1),1);
TWWFVls = zeros(size(freepara2,1),4313);

%outer loop parameter (theta in this case)
lower = 0.0028; %0.0028 opt
upper = 0.005; %0.005 opt
step = 0.0002; %0.0002 opt
freepara = [lower:step:upper];
IncrVls = zeros(size(freepara,2),size(freepara2,2));
JintraVls = zeros(size(freepara,2),size(freepara2,2));
JinterVls = zeros(size(freepara,2),size(freepara2,2));





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%calculate the gain coefficient for a range of phonon energies
lower3 = 1e9;
upper3 = 1000e9;
step3 = 10e9;
phononFreq = [lower3:step3:upper3];
gainVlsOmega = zeros(size(phononFreq,2),size(IncrVls,2));

index3 = 1;
for phononFreqloop = lower3:step3:upper3;
    
omega = phononFreq(1,index3)*2*pi;    
       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    




NumberOfSteps = (upper-lower)/step
l = 1;
for theta = lower:step:upper;
    l
    l2 = 1;
    for Delta = lower2:step2:upper2;
    
        %get the wavefunction for a single quantum well and plot it
        SWWF = Single_well_wavefunction(hbar,me,mw,mb,Egw,Egb,a,b,d,V0,E);

        %get the overlap integral, t, to use in the two well model and calculate
        %the coupling factor
        OVLP = OverlapIntegral(SWWF,Egw,Egb,V0,a,b,d,e);
        t = OVLP; %in J
        CoupleFactor = t/Delta;

        %get the two well model (coupled well) wavefunction
        TWWF = TwoWellModelWF(SWWF,Egw,Egb,V0,a,b,d,e,Delta,CoupleFactor);

        %figure;
        %hold on;
        %plot(TWWF(:,1),TWWF(:,2),'bx');
        %plot(SWWF(:,1),SWWF(:,2),'r');
        %title('Single well wavefunction (red) and the coupled well wavefunction (blue)')

        %the phonon wave vector
        q = omega / vsSL;

        %use the TWWF to obtain the form factor Jintra
        Jintra = JintraFunc(TWWF,q,theta);
  
        %calculate the form factor Jinter
        Jinter = JinterFunc(TWWF,q,theta,d,Jintra,CoupleFactor);

        %calculate the quasi Fermi level
        ef1 = exp(((hbar^2)*pi*ne) / (mw*me*Te));
        Ef = Te*log(ef1 - 1);
        Ef = abs(Ef);
    
        %get P for Delta-(hbar omega)
        epsilon1 = Delta - (hbar*omega);
        P1 = P(epsilon1,Te,Ef,omega,rhow,vsw,q,theta,E1,mw,mb,me,hbar,d);
    
        %get P for Delta+(hbar omega)
        epsilon2 = Delta + (hbar*omega);
        P2 = P(epsilon2,Te,Ef,omega,rhow,vsw,q,theta,E1,mw,mb,me,hbar,d);
    
        %get P for -Delta-(hbar omega)
        epsilon3 = -Delta-(hbar*omega);
        P3 = P(epsilon3,Te,Ef,omega,rhow,vsw,q,theta,E1,mw,mb,me,hbar,d);
    
        %get P for (hbar omega)-Delta
        epsilon4 = (hbar*omega)-Delta;
        P4 = P(epsilon4,Te,Ef,omega,rhow,vsw,q,theta,E1,mw,mb,me,hbar,d);
    
        %get P for -(hbar omega)
        epsilon5 = -hbar*omega;
        P5 = P(epsilon5,Te,Ef,omega,rhow,vsw,q,theta,E1,mw,mb,me,hbar,d);
    
        %get P for (hbar omega)
        epsilon6 = hbar*omega;
        P6 = P(epsilon6,Te,Ef,omega,rhow,vsw,q,theta,E1,mw,mb,me,hbar,d);
        
        %get the screening factor
        kappa = SreeningFactorFunc(e,ne,q,theta,d,RelPerm,AbsPerm,me,mw,a,omega);
        
        %calculate the increment
        Incrmnt =( (Jinter*(P1-P2+P3-P4)) + (Jintra*(P5-P6)) ) / kappa^2; %Jintra non-zero
        %Incrmnt = Jinter*(P1-P4) / kappa^2; %Jintra = 0
    
        %put the values to be plotted in the pre defined vectors
        IncrVls2(l2,1) = Incrmnt;
        %JinterVls2(l2,1) = Jinter;
        %JintraVls2(l2,1) = Jintra;
        %TWWFMatrix2(:,l2) = TWWF(:,2);
        l2 = l2 + 1;
              
    end
   
    IncrVls(l,:) = IncrVls2(:,:);
    l = l + 1;
        
end

%take the values for increment and calculate the gain
%for a single point directly under the center of the detector
thetaVls = freepara;
qparaVls = q*sin(thetaVls);

%account for angular spread on the flat detector
Spreadtheta = ((cos(thetaVls)).^2.*2.*pi);
Spreadtheta = Spreadtheta';

%normalize theta to the full theta space
thetaVlsNorm = thetaVls(:,:) ./ (upper-lower);

%gain values (the number of phonons after 1sec seeding with 1 phonon)
gainVls  = zeros(size(IncrVls,2),1);
counter = 1;
for counter = 1:size(IncrVls,2)
    IncrVls(:,counter) = 2*pi*IncrVls(:,counter) ./Spreadtheta; %IncrVls*2*pi for full circle
    gain = trapz(thetaVlsNorm,IncrVls(:,counter)); %with respect to normalized thetaVls
    gainVls(counter,1) = gainVls(counter,1) + gain;
    counter = counter + 1;
end

%calculate the gain coefficient for all Delta
gainVls(:,:) = gainVls(:,:) .* (n*d) ./ vsSL; %from gain/sec to gain over the SL
gainVls(:,:) = gainVls(:,:) + 1; %add the origonal seed phonon to the gain
gainVls(:,:) = log(gainVls(:,:)) ./ (n*d) ./ 100; %gain coefficent (cm^{-1})






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%setup the plank spectrum
PlankSpec = 1 / ( exp( (hbar*omega) / (kb*Th) ) - 1 );
%apply it to the gain values and put them into a matrix of Delta against
%omega for all angles
gainVlsOmega(index3,:) = gainVls(:,:) * PlankSpec;
index3 = index3 + 1;

end

figure;
surf(DeltameV,phononFreq,gainVlsOmega);
shading 'interp'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%Plot and output results%%%%

%convert parameters to the correct units
DeltameV = freepara2 ./ 1.6e-19 ./1e-3; %Stark-splitting in meV
qpara = q.*sin(thetaVls); %theta as wavevector

%plot the gain against Delta
figure
plot(DeltameV,gainVls)
xlabel 'Stark-splitting (meV)'
ylabel 'gain coefficient (cm^{-1})'

%3D plot of the increment
figure;
hold on;
surf(DeltameV,qpara,IncrVls)
shading 'interp'
xlabel 'Stark splitting (meV)'
ylabel 'q parrallel (m^{-1})'
zlabel 'Increment (s^{-1})'

%plot the two well wavefunction as a function of the Stark splitting
%TWWFMatrix(1:1500,:) = [];
%TWWFMatrix(1500:end,:) = [];
%figure;
%surf(TWWFMatrix);
%shading 'interp';
%clear TWWFMatrix

%output to a file (Increment Matrix)
Outfname = ['PATH_HERE\Increment.txt'];
fid = fopen(Outfname, 'w'); 
fprintf(fid, '%d\t', DeltameV);
in1 = 1;
for in1 = 1:size(IncrVls,1);
    in2 = 1;
    for in2 = 1:size(IncrVls,2)
        fprintf(fid, '%d\t', IncrVls(in1,in2));
        in2 = in2 + 1;
    end
    fprintf(fid, '\n');
end
status = fclose(fid);

%output to a file (gain)
Outfname = ['PATH_HERE\Gain.dat'];
fid = fopen(Outfname, 'w'); 
in = 1;
    for in = 1:size(DeltameV,2);
        fprintf(fid, '%d\t', DeltameV(1,in));
        fprintf(fid, '%d\n', gainVls(in,1));
        in = in + 1;
    end
status = fclose(fid);



