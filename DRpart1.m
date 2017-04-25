
%setting the values of R and C
R = 20;
C = 10e-6;

%transfer function of the low pass filter
fn = @(f) ((1i*2*pi*R*C*f)+1).^(-1);


%bandwidth of the circuit
bw=(2*pi*R*C)^(-1);

%plotting the transfer function
f = linspace(0,1e4,1000); %frequency vector
Hpoints = fn(f);
figure (1)
loglog(f,real(Hpoints));
grid on;

%finding the differential equation 
syms t s;
Vi = 1;
VoMain = 0;
F = Vi * ((s*R*C)+1)^(-1);
VoMain = @(t) Vi * (1 - exp(-t/(R*C)));

%plotting Vo vs time
timeMain = linspace(0,1e-3);
figure (2)
plot(timeMain, VoMain(timeMain));
grid on;

%three different values of timesteps (delta t) and their number of steps
numA = length(0:1e-6:1e-3);
numB = length(0:1e-5:1e-3);
numC = length(0:1e-4:1e-3);
t1=linspace(0,1e-3,numA);
t2=linspace(0,1e-3,numB);
t3=linspace(0,1e-3,numC);
VoMainD1 = zeros(1,numA);
VoMainD2 = zeros(1,numB);
VoMainD3 = zeros(1,numC);

%the equation is implicit because it uses both Vo(n) and Vi, ie, both the
%dependent and indeprendent variables
%approximating the derivative using first principles
for n=1:numA-1
%approximating the derivitative using a first order difference equation
VoMainD1(n+1) = ((Vi-VoMainD1(n))*(1e-6)/(R*C))+VoMainD1(n);
end
for n=1:numB-1
%approximating the derivitative using a first order difference equation
VoMainD2(n+1) = ((Vi-VoMainD2(n))*(1e-5)/(R*C))+VoMainD2(n);
end
for n=1:numC-1
%approximating the derivitative using a first order difference equation
VoMainD3(n+1) = ((Vi-VoMainD3(n))*(1e-4)/(R*C))+VoMainD3(n);
end
figure(3)
plot(t1, VoMainD1, t2, VoMainD2, t3,VoMainD3);
grid on;
%Vo becomes unstable when delta is 1e-4 or smaller


%Testing with a sinusoidal input
freq1 = 1000; % Enter value of f, and experiment by changing it
freq2 = 10000;
freq3 = 100000;
Lsignal = 1e-2;
time = linspace(0,Lsignal,100);
VoSine1 = zeros(1,numB);
VoSine2 = zeros(1,numB);
VoSine3 = zeros(1,numB);
for n=1:numB-1
VoSine1=((sin(2*freq1*time)-VoMainD2(n))*(1e-5)/(R*C))+VoMainD2(n);
VoSine2=((sin(2*freq2*time)-VoMainD2(n))*(1e-5)/(R*C))+VoMainD2(n);
VoSine3=((sin(2*freq3*time)-VoMainD2(n))*(1e-5)/(R*C))+VoMainD2(n);
end
figure (4)
plot(time,VoSine1, time,VoSine2, time,VoSine3);

grid on;

%frequency responce plot using fft
%using freq1 1000, thus Fs has to be at least 2000
Fs1 = 2e3;
SamplingP = 1/Fs1;
%lng = length(0:SamplingP:freq1-SamplingP);

fresp = fft(VoSine1);
magnitude = abs(fresp);
fresp2 = fft(VoSine2);
magnitude2 = abs(fresp2);
figure (5)
plot(fresp,magnitude,fresp2,magnitude2);
grid on;
%%
%part 2, finding the mean square voltage
K = 1.38064852e-23;
T = 298;
VRMS = (4*K*T*R*bw)^(0.5);
%%
%part 3;
Imax = Vi/5*R;
rng(0,'twister');
Inoise = Imax*randn;
VoutN = zeros(1,numA);
Vin = 1;
for n=1:numA-1
%approximating the derivitative using a first order difference equation
Inoise = Imax*randn;
VoutN = ((((Vin-VoutN(n))/(R*C))+Inoise/C)*(1e-6))+VoutN(n);
end
figure (6)
plot(t2,VoutN);
Vin = 0;
Vrms = sqrt(mean(VoutN)^2);


