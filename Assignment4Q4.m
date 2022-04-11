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


%% Transient step response
t = 0:0.01:1;
Fnew = [0;0;0;0;0;0;0;0];
inputV = zeros(1,size(t,2));
outputV5 = zeros(size(GMatrix,1),size(t,2)); 

for i = 1:size(t,2)-1

    if i<3
        inputV(i+1) = 0;
    else
        inputV(i+1) = 1;
    end

    outputV5(:,i) = Fnew;
    Fnew = ((CMatrix./0.001) + GMatrix) \ ((CMatrix/0.001) * Fnew + F * inputV(i+1));
    outputV5(:,i+1) = Fnew;
end

figure(1)
plot(t,inputV)
hold on
plot(t,outputV5(5,:))
title('Step Response');
xlabel('Time (s)');
ylabel('Voltage (V)');

figure(2)
semilogy(abs(fftshift(fft(outputV5(5,:)))))
title('Step Response Fourier');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% Transient sine function
t = 0:0.001:1;
Fnew = [0;0;0;0;0;0;0;0];
inputV = zeros(1,size(t,2));
w = 2 * pi * (1/0.03);
inputV(1,1) = 0;
outputV5 = zeros(size(GMatrix,1),size(t,2)); 

for i = 1:size(t,2)-1

    inputV(1,i+1) = sin(w*t(i+1));

    outputV5(:,i) = Fnew;
    Fnew = ((CMatrix./0.001) + GMatrix) \ ((CMatrix/0.001) * Fnew + F * inputV(i+1));
end

figure(3)
plot(t,inputV)
hold on
plot(t,outputV5(5,:))
title('Sine Wave Response');
xlabel('Time (s)');
ylabel('Voltage (V)');

figure(4)
semilogy(abs(fftshift(fft(outputV5(5,:)))))
title('Sine Wave Response Fourier');
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%% Transient gaussian pulse

t = 0:0.001:1;
Fnew = [0;0;0;0;0;0;0;0];
inputV = zeros(1,size(t,2));
std = 0.03;
w = 2*pi*std;
inputV(1,1) = 0;
outputV5 = zeros(size(GMatrix,1),size(t,2)); 

for i = 1:size(t,2)-1
    if i<60
        inputV(i+1) = 0;
    else
        inputV(1,i+1) = 1/sqrt(w) * exp(-(t(i+1).^2)/(2*std^2));
    end
    
    outputV5(:,i) = Fnew;
    Fnew = ((CMatrix./0.001) + GMatrix) \ ((CMatrix/0.001) * Fnew + F * inputV(i+1));
end

figure(5)
plot(t,inputV)
hold on
plot(t,outputV5(5,:))
title('Gaussian Wave Response');
xlabel('Time (s)');
ylabel('Voltage (V)');

figure(6)
semilogy(abs(fftshift(fft(outputV5(5,:)))))
title('Gaussian Wave Response Fourier');
xlabel('Frequency (Hz)');
ylabel('Amplitude');