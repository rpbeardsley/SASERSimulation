%Calculate the wave function in a single quantumm well

%variables
%me -electron mass
%mw -effective electron mass in well
%mb -effective electron mass in the barrier
%Egw -energy gap of the well at 10K
%Egb -energy gap of the barrier at 10K
%a -well width
%b -barrier width
%d -SL period
%V0 -the barrier height
%E -energy of the first confined state

function SingleWellWF = Single_well_wavefunction(hbar,me,mw,mb,Egw,Egb,a,b,d,V0,E);

k = ((E*(2*me*mw))^0.5)/hbar; %k in the well
K = ((E*(2*me*mb))^0.5)/hbar; %K outside the well

%evaluating the normalised wavefunction at the boundaries
alpha = cos(k*a/2)/exp(-K*a/2);

den = ( (alpha^2) * 2*(exp(-K*a)/K) ) + a/2 + sin(k*a)/(4*k) - sin(-k*a)/(4*k);

%calculate the coefs for the wave functions
C = (1/den)^0.5;
D = -alpha*C;

%output the wave functions over the well
Z1 = [-1e-6:0.05e-9:-2.95e-9]';
wf1 = D*exp(Z1.*K);

Z2 = [-2.75e-9:0.05e-9:2.75e-9]';
wf2 = C*cos(Z2.*k);

Z3 = [2.95e-9:0.05e-9:1e-6]';
wf3 = D*exp(-Z3.*K);

wf = vertcat(wf1,wf2,wf3);
Z = vertcat(Z1,Z2,Z3);

SingleWellWF = horzcat(Z,wf);
%SingleWellWF(:,:) = abs(SingleWellWF(:,:));
%check normalisation
ProbDen = abs(SingleWellWF(:,2)).^2;
norm = trapz(Z,ProbDen) ;% = 1

%figure;
%hold on;
%plot(Z,SingleWellWF(:,2),'x');
%plot(Z,ProbDen);






















