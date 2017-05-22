%obtain the overlap integral

%variables
%Egw -energy gap of the well at 10K
%Egb -energy gap of the barrier at 10K
%V0 -the barrier height
%a -well width
%b -barrier width
%d -SL period
%e -elementary charge
%wfw1 -single well wavefunction

function tJ = OverlapIntegral(wfw1,Egw,Egb,V0,a,b,d,e);

%put the wavefunction into the neighbouring well
Zshift = wfw1(:,1) + d;
wfw2 = horzcat(Zshift,wfw1(:,2));  %well2

%select the region over which to integrate and resize the matricies
%wavefunction well1
[trash, array_position1] = min(abs(wfw1 - -2*d)); %middle of next nearest neighbour of well1
wfw1(1:array_position1(:,1),:) = [];
[trash, array_position2] = min(abs(wfw1 - 3*d)); %middle of next nearest neighbour of well2
wfw1(array_position2(:,1):end,:) = [];

%wavefunction well2
[trash, array_position3] = min(abs(wfw2 - -2*d)); %middle of neighbour of well1
wfw2(1:array_position3(:,1),:) = [];
[trash, array_position4] = min(abs(wfw2 - 3*d)); %middle of neighbour of well2
wfw2(array_position4(:,1):end,:) = [];

%take Z and multiply the two wfs
Z = wfw1(:,1);

ProdWF = wfw1(:,2) .* wfw2(:,2);

%plot the two wave functions in the neigbouring wells
%figure;
%hold on;
%plot(Z,wfw1(:,2),'b');
%plot(Z,wfw2(:,2),'r');

%figure;
%plot(Z,ProdWF,'x');

%try integrating over different regions
ind = find(Z >= 2.95e-9  & Z <= 6.85e-9);

%Integrate over the product of the WFs
tInt = trapz(Z,ProdWF);

%multiply by the barrier height to find the overlap integral t
tJ = V0*tInt; % in J
tmeV = tJ*1000/e; % in meV




