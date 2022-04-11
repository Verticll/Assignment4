%Capacitor
C=0.25;
%Resistors
G0 = 1/1000;
G1 = 1;
G2 = 1/2;
G3 = 1/0.213;
G4 = 1/0.1;
%Inductor
L = 0.2;
%Ratio
alpha = 100;

V = alpha*G3;

%% BUILDING MATRICIES

% v1, v2, v3, v4, v5, vin, i4, i1
%  0,  0,  0,  0,  0,  0,  0,  0; v1
%  0,  0,  0,  0,  0,  0,  0,  0; v2
%  0,  0,  0,  0,  0,  0,  0,  0; v3
%  0,  0,  0,  0,  0,  0,  0,  0; v4
%  0,  0,  0,  0,  0,  0,  0,  0; v5
%  0,  0,  0,  0,  0,  0,  0,  0; vin
%  0,  0,  0,  0,  0,  0,  0,  0; i4
%  0,  0,  0,  0,  0,  0,  0,  0; i1

CMatrix = [-C, C, 0, 0, 0, 0, 0, 0; 
            C,-C, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, 0;
            0, 0, 0, 0, 0, 0, 0, -L;];

GMatrix = [-G1, G1, 0, 0, 0,-1, 0, 0;
            G1,-G1-G2, 0, 0, 0, 0, 0,-1;
            0, 0,-G3, 0, 0, 0, 0, 1;
            0, 0, 0,-G4, G4, 0, 1, 0; 
            0, 0, 0, G4,-G4-G0, 0, 0, 0;
            1, 0, 0, 0, 0, 0, 0, 0; 
            0, 0,-V, 1, 0, 0, 0, 0;
            0, 1,-1, 0, 0, 0, 0, 0;];

F = [0;0;0;0;0;1;0;0];

%% DC SWEEP
dcV = -10:0.1:10;
outputV3 = zeros(size(dcV));
outputV5 = zeros(size(dcV));

for i = 1:size(dcV,2)
    V = (GMatrix + CMatrix)\(F.*dcV(i));
    outputV3(i) = V(3);
    outputV5(i) = V(5);
end

figure(1)
plot(dcV,outputV3);
hold on
plot(dcV,outputV5)
title('DC Sweep');
xlabel('Input Voltage(V)');
ylabel('Output Voltage(V)');

%% AC SWEEP
ACfreq = linspace(0.01, 100);
outputV3 = zeros(size(ACfreq));
outputV5 = zeros(size(ACfreq));
for i = 1:size(ACfreq,2)
    w = 2*pi*ACfreq(i);
    jw = 1i * w;
    V = (GMatrix + jw * CMatrix)\F;
    outputV3(i) = V(3);
    outputV5(i) = V(5);
end

figure(2)
plot(2*pi*ACfreq,db(real(outputV3),'voltage'));
hold on
plot(2*pi*ACfreq,db(real(outputV5),'voltage'))
title('AC Sweep');
xlabel('Input Frequency(rad/s)');
ylabel('Output Gain(db)');

%% CAPACITOR
RandomC = C + 0.05 * randn(10000,1);
outputV5 = zeros(size(RandomC));
for i = 1:size(RandomC)
    CMatrix = [-RandomC(i), RandomC(i), 0, 0, 0, 0, 0, 0; 
                RandomC(i),-RandomC(i), 0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0, 0, 0;
                0, 0, 0, 0, 0, 0, 0, -L;];
    V = (GMatrix + (1i*pi)*CMatrix)\F;
    outputV5(i) = V(5);
end
figure(3)
histogram(RandomC);
title('Histogram of Capacitance');
xlabel('Capacitance(F)');
ylabel('Bins');

figure(4)
histogram(real(outputV5));
title('Histogram of Output');
xlabel('Vout(V)');
ylabel('Bins');
