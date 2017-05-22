%calculate the wave function using the two well model

%variables
%t -Overlap integral
%Delta -Stark splitting
%CoupleFactor -t/Delta
%wfw1 -single well wavefunction
%Egw -energy gap of the well at 10K
%Egb -energy gap of the barrier at 10K
%V0 -the barrier height
%a -well width
%b -barrier width
%d -SL period
%e -elementary charge

function Wavefunction = TwoWellModelWF(wfw1,Egw,Egb,V0,a,b,d,e,Delta,CoupleFactor) 

%get the wavefunctions of the neighbouring wells (adjust for the
%translation along Z)
Zshiftw2 = wfw1(:,1) + d; %well 2
wfw2 = horzcat(Zshiftw2,wfw1(:,2)); %well 2
Zshiftw0 = wfw1(:,1) - d; %well 0
wfw0 = horzcat(Zshiftw0,wfw1(:,2)); %well 0

%figure;
%hold on;
%plot(wfw1(:,1),wfw1(:,2),'b');
%plot(wfw2(:,1),wfw2(:,2),'r');
%plot(wfw0(:,1),wfw0(:,2),'g');

%put them all on the same Z axis
%well 1
[trash, array_position1] = min(abs(wfw1 - -10*d)); %middle of next nearest neighbour of well1
wfw1(1:array_position1(:,1),:) = [];
[trash, array_position2] = min(abs(wfw1 - 10*d)); %middle of next nearest neighbour of well2
wfw1(array_position2(:,1):end,:) = [];

%well 2
[trash, array_position3] = min(abs(wfw2 - -10*d)); %middle of neighbour of well1
wfw2(1:array_position3(:,1),:) = [];
[trash, array_position4] = min(abs(wfw2 - 10*d)); %middle of neighbour of well2
wfw2(array_position4(:,1):end,:) = [];

%well 0
[trash, array_position5] = min(abs(wfw0 - -10*d)); %middle of next nearest neighbour of well1
wfw0(1:array_position5(:,1),:) = [];
[trash, array_position6] = min(abs(wfw0 - 10*d)); %middle of next nearest neighbour of well2
wfw0(array_position6(:,1):end,:) = [];

%take the Z values
Z = wfw1(:,1);

%plot the wavefunctions in the neigbouring wells to be sure they are correct
%figure;
%hold on;
%plot(Z,wfw1(:,2),'b');
%plot(Z,wfw2(:,2),'r');
%plot(Z,wfw0(:,2),'g');

%form the two well model wavefunction
A = wfw1(:,2);
B = CoupleFactor; %t/Delta;
C = wfw2(:,2) - wfw0(:,2);

WavefunctionY = A - (B.*C);
Wavefunction = horzcat(Z,WavefunctionY);

%figure
%hold on
%plot(Z,Wavefunction)
























