
% Patrobas Adewumi
% 100963608
% ELEC 4700
% assignment4

disp('Question 2(a) It can be seen from the frequency response plot that this is a low pass filter circuit.')
disp('Increase in the frequency increases the inductance making the capacitor to cause a parallel resistance')
disp('with the resistor')
 
disp('Question 2(b) The expected frequency response would be a pass band during the low frequencies')
disp(' and a 2nd order drop off at its respective frequency')

clear all

runTime = 1; 
timecuts = 1000;
dt = runTime/timecuts;

R1 = 1;
C1 = 0.25;
R2 = 2;
L1 = 0.2;
R3 = 10;
a = 100;
R4 = 0.1;
Ro = 1000;     

C = [  0, 0,   0, 0, 0, 0, 0; ...
     -C1,C1,   0, 0, 0, 0, 0; ...
       0, 0, -L1, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0; ...
       0, 0,   0, 0, 0, 0, 0];
   
G = [     1,             0,  0,    0,  0,     0,            0; ...
      -1/R1, (1/R2 + 1/R1), -1,    0,  0,     0,            0; ...
          0,             1,  0,   -1,  0,     0,            0; ...
          0,             0, -1, 1/R3,  0,     0,            0; ...
          0,             0,  0,    0, -a,     1,            0; ...
          0,             0,  0, 1/R3, -1,     0,            0; ...
          0,             0,  0,    0,  0, -1/R4, (1/R4 +1/Ro)];

disp('Question 2: G and C matrices:')
disp(G)
disp(C)
 
V1 = 0;      
F = zeros(7,1);

Flist = zeros(7,1,timecuts);
Flist(1,1,30:timecuts) = 1;
Vlist = zeros(7,1,timecuts);

for count  = 2:1:timecuts
    A = C/dt +G;
    Vlist(:,:,count) = A\(C*Vlist(:,:,count-1)/dt +Flist(:,:,count));
end
      
V1list(1,:) = Vlist(1,1,:);
V2list(1,:) = Vlist(2,1,:);
ILlist(1,:) = Vlist(3,1,:);
I3list(1,:) = Vlist(4,1,:);
V4list(1,:) = Vlist(5,1,:);
Volist(1,:) = Vlist(7,1,:);


figure(1)
plot((1:timecuts).*dt,Volist(1,:))
xlabel('Time(secs)')
ylabel('Voltage')
title('Vout of Step voltage')
hold on
plot((1:timecuts).*dt,V1list(1,:))
hold off

figure(2)
g = abs(fftshift(fft(Volist(1,:))));
plot(((1:length(g))/timecuts)-0.5,g)
xlim([-0.05 0.05])
xlabel('Freq (Hz)')
ylabel('Magnitude')
title('Fourier Transform of Output')
      
% Sine function
f = 1/0.03;
Flist = zeros(7,1,timecuts);
for count = 1:1:timecuts
    Flist(1,1,count) = sin(2*pi*f*count*dt);
end
Vlist = zeros(7,1,timecuts);

for count  = 2:1:timecuts
    A = C/dt +G;
    Vlist(:,:,count) = A\(C*Vlist(:,:,count-1)/dt +Flist(:,:,count));
end
      
V1list(1,:) = Vlist(1,1,:);
V2list(1,:) = Vlist(2,1,:);
ILlist(1,:) = Vlist(3,1,:);
I3list(1,:) = Vlist(4,1,:);
V4list(1,:) = Vlist(5,1,:);
Volist(1,:) = Vlist(7,1,:);



figure(3)
plot((1:timecuts).*dt,Volist(1,:))
xlabel('Time(secs)')
ylabel('Voltage')
title('Vout of Sine wave')
hold on
plot((1:timecuts).*dt,V1list(1,:))
hold off

figure(4)
g = abs(fftshift(fft(Volist(1,:))));
plot(((1:length(g))/timecuts)-0.5,g)
xlim([-0.05 0.05])
xlabel('Freq (Hz)')
ylabel('Magnitude')
title('Fourier Transform of Output')

% Gaussian pulse
% Vgauss = exp(-1/2*((k/ts-0.06)/(0.03))^2)
mag = 1;
dev = 0.03;
delay = 0.06;
Flist = zeros(7,1,timecuts);
for count = 1:1:timecuts
    Flist(1,1,count) = exp(-((count*dt-0.06)/0.03)^2);
end
Vlist = zeros(7,1,timecuts);

for count  = 2:1:timecuts
    A = C/dt +G;
    
    Vlist(:,:,count) = A\(C*Vlist(:,:,count-1)/dt +Flist(:,:,count));
end
      
V1list(1,:) = Vlist(1,1,:);
V2list(1,:) = Vlist(2,1,:);
ILlist(1,:) = Vlist(3,1,:);
I3list(1,:) = Vlist(4,1,:);
V4list(1,:) = Vlist(5,1,:);
Volist(1,:) = Vlist(7,1,:);



figure(5)
plot((1:timecuts).*dt,Volist(1,:))
xlabel('Time(secs)')
ylabel('Voltage')
title('Vout of Gaussian Pulse')
hold on
plot((1:timecuts).*dt,V1list(1,:))
hold off

figure(6)
g = abs(fftshift(fft(Volist(1,:))));
plot(((1:length(g))/timecuts)-0.5,g)
xlim([-0.05 0.05])
xlabel('Freq (Hz)')
ylabel('Magnitude')
title('Q2: Fourier Transform of Output')

disp('A change in the timecut changes the accuracy of the simulation. Hence, a larger timecut results in less accurate sims')
