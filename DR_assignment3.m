close all
clc

%% 
%Constants ofassignment 1
mass=9.1e-31;
K=1.38*10e-23;  %Boltzmann
m=0.26*(mass); %effective mass of electrons
w=200e-9; %nominal size of region
l=100e-9;
%Calculating the thermal velocity
T=300; %temperature
vth=sqrt((T*K)/m); 
%Mean free path calculation
meanT=2e-13;
mfp=meanT*vth;
%Modeling the random motion of the electrons
npar=1000; %number of particles/electrons
%setting up matrix for location and velocity (i,1) for x and (i,2) for y
locs = zeros(npar,2);
vels = zeros(npar,2);
%Assignning particles random locations in the x-y plane
locs(:,1)=rand(npar,1)*w;
locs(:,2)=rand(npar,1)*l;
%Assignning particles to move in random directions but with the thermal
 %velocity
vels(:,1)=randn(npar,1)*vth;
vels(:,2)=randn(npar,1)*vth;
%setting up the spacial step
deltaT = l/100/vth;
timesteps = 1000;
olocs = locs;
%modelling the scattering of the electrons using an exponential
%scattering probability
pscat= (1-(exp(-deltaT/meanT))); 
temp = zeros(npar,timesteps);

%This assignment presents the electrons in the silicon as particles with the 
%effective mass, using a simplistic Monte-Carlo model just like in assignment 1. 
%However, a voltage source of 0.1V is applied and causes an electric field of  V/m 
%to be generated. The force on an electron is N due to the electric field, and it 
%also accelerates the motion of the electrons. The trajectory of the particles below 
%shows how the electric field applied causes it to curve at certain occasion

Vx=0.8;
q=1.60217662e-19;
elecF=Vx/w;
F=elecF*q;
accx=F/m;
accy = 0;

eleCon=1e15;
I=zeros(npar,1);

for n = 1:timesteps
        if pscat>rand(n)
            vels(:,1)=randn(npar,1)*vth;
            vels(:,2)=randn(npar,1)*vth;
        end
        %updating the particle location using Newton's law of motion
        vels(:,1)= vels(:,1)+(accx*deltaT);
        vels(:,2)= vels(:,2)+(accy*deltaT);
        locs(:,1) = locs(:,1)+ vels(:,1)*deltaT;
        locs(:,2) = locs(:,2)+ vels(:,2)*deltaT;
      
        %setting the boundary condition in the y direction for velocity and
        %position    
        shiftInY=locs(:,2)<0;
        locs(shiftInY,2)= -locs(shiftInY,2); 
        vels(shiftInY,2)= -vels(shiftInY,2);

        shiftInY=locs(:,2)>l;
        locs(shiftInY,2)=2*l-locs(shiftInY,2); 
        vels(shiftInY,2)= -vels(shiftInY,2);   

        %setting the boundary condition in the x direction
        shiftInX=locs(:,1)<0;
        parsLeaving = sum(shiftInX);
        locs(shiftInX,1)= locs(shiftInX,1)+w; 
        shiftInX=locs(:,1)>w;
        parsEntering = sum(shiftInX);
        locs(shiftInX,1)= locs(shiftInX,1)-w; 
        
        I(n) = (parsLeaving-parsEntering);
        

        s = rand(npar,1) < pscat;
        vels(s,1)=randn(sum(s),1)*vth;
        vels(s,2)=randn(sum(s),1)*vth;
        
        temp(:,n) = mean((1/(2*K))*m*((vels(:,1).^2)+(vels(:,2).^2)));
        
        X=[olocs(:,1) locs(:,1)];
        Y=[olocs(:,2) locs(:,2)];
        pause(0.000001)
        plot (X,Y, 'x'); hold on
        title ('trajectory part A');
        olocs = locs;

       
end

%Current is computed through finding the number of electrons crossing a line
%with a fixed x-coordinate value. The number of electrons is the difference 
%between the ones entering and leaving. To get the current, the number of 
%electrons was multiplied by the charge of an electron ?q? and divided by the 
%length of the line being crossed. 

Time = (0:deltaT:deltaT*(timesteps-1));
figure
plot (Time, q*I(:,1)/l);
title('Current')
xlabel('Time (s)')
ylabel('Current (A)')
 
%density and temperature
%density
delta = 1e-8;
ndx = w/delta;
ndy = l/delta;
density=zeros(ndx,ndy);
tempMap = zeros(ndx,ndy);

for n=1:ndx
    for k=1:ndy
        ni=locs(:,2)<(delta*k) &...
                 locs(:,2)>(delta*(k-1)) & locs(:,1)<(delta*n)...
                 & locs(:,1)>(delta*(n-1)); 
             density(n,k)=sum(ni);
             tempMap(n,k) = mean(((1/(2*K))*m*(vels(ni,1).^2 + vels(ni,2).^2)));       
             plot(locs(ni,1),locs(ni,2),'x');hold on
    end            
end

figure
surf(density');pause(0.001)
title('Electron Density Part A');
grid on

figure
surf(tempMap');pause(0.001)
title('Temperature Map Part A');
colorbar
grid on
%%
% Constants of assignment 2
w=200e-9; %length of x
l=100e-9;%length of y
nx=200; %number of steps in the x-direction
ny=100; %number of steps in the y-direction
deltaX=w/(nx-1); %interval between steps
deltaY=l/(ny-1);
V0=1; %inital value of V0

% set up of assignment 2
x = linspace(0,w,nx);
y = linspace(0,l,ny);
sig = zeros(nx,ny); %conductivity matrix
V = zeros(nx,ny);

%code of assignment 2

for i=1:1:nx
    for j=1:1:ny 
        %measurements adjusted to the dimensions of the bottleneck in
        %assg1/3
        if((x(i)>=2*w/5 && x(i)<=3*w/5) && (y(j)<=2*l/5 || y(j)>=3*l/5))
            sig(i,j)=1e-2;
        else
            sig(i,j)=1;
        end
    end
end

e=1; %max error
while(e>=2e-6) 
    V_old=V;
    for i=1:1:nx
        for j=1:1:ny
            if(i==nx)
                V(i,j)=0;
            elseif (i==1)
                V(i,j)=V0;
            elseif (j==ny)
                V(i,j)=V(i,j-1);
            elseif (j==1)
                V(i,j)=V(i,j+1);
            else
                rw1=((1/sig(i,j))+(1/sig(i,j+1)))/2;
                rw2=((1/sig(i,j))+(1/sig(i,j-1)))/2;
                rw3=((1/sig(i,j))+(1/sig(i+1,j)))/2;
                rw4=((1/sig(i,j))+(1/sig(i-1,j)))/2;
                
                V(i,j)=1/((1/rw1)+(1/rw2)+(1/rw3)+(1/rw4))*...
                    ((V(i,j+1)/rw1)+(V(i,j-1)/rw2)+(V(i-1,j)/rw3)+(V(i+1,j)/rw4));
                
            end
        end
    end       
    e = max(max(abs(V - V_old)));
end

%Plotting Potential
figure
surf(V');pause(0.001)
title('Voltage plot Part B');
ylabel('Y Position');
xlabel('X Position');
zlabel('Voltage');
grid on
%%
%Calculating and plotting Electric field
[EyNg,ExNg] = gradient(V);			
Ex = -ExNg;							
Ey = -EyNg;

% figure
% quiver(Ex',Ey');
% pause(0.001)
% title('Voltage plot');
% ylabel('Y Position');
% xlabel('X Position');
% zlabel('Electrical field');
% grid on

%calculating the acceleration using the electric fields
mass=9.1e-31;
m=0.26*(mass); %effective mass of electrons
q=1.60217662e-19;

accx = q*Ex/m;
accy = q*Ey/m;


for n = 1:timesteps
    
        if pscat>rand(n)
            vels(:,1)=randn(npar,1)*vth;
            vels(:,2)=randn(npar,1)*vth;
        end
        for j=1:npar
            axp=abs(floor(locs(j,1)/deltaX))+1;
            ayp=abs(floor(locs(j,2)/deltaY))+1;
            axxx=accx(axp,ayp);
            ayyy=accy(axp,ayp);

            %updating the particle location using Newton's law of motion
            vels(:,1)= vels(:,1)+(accx(axp,ayp)*deltaT);
            vels(:,2)= vels(:,2)+(accy(axp,ayp)*deltaT);
        end
        locs(:,1) = locs(:,1)+ vels(:,1)*deltaT;
        locs(:,2) = locs(:,2)+ vels(:,2)*deltaT;
      
        %setting the boundary condition in the y direction for velocity and
        %position    
        shiftInY=locs(:,2)<0;
        locs(shiftInY,2)= -locs(shiftInY,2); 
        vels(shiftInY,2)= -vels(shiftInY,2);

        shiftInY=locs(:,2)>l;
        locs(shiftInY,2)=2*l-locs(shiftInY,2); 
        vels(shiftInY,2)= -vels(shiftInY,2);   

        %setting the boundary condition in the x direction
        shiftInX=locs(:,1)<0;
        locs(shiftInX,1)= locs(shiftInX,1)+w; 
        shiftInX=locs(:,1)>w;
        locs(shiftInX,1)= locs(shiftInX,1)-w; 
        
        s = rand(npar,1) < pscat;
        vels(s,1)=randn(sum(s),1)*vth;
        vels(s,2)=randn(sum(s),1)*vth;
        
        temp(:,n) = mean((1/(2*K))*m*((vels(:,1).^2)+(vels(:,2).^2)));
        
        X=[olocs(:,1) locs(:,1)];
        Y=[olocs(:,2) locs(:,2)];
        pause(0.000001)
        plot (X,Y, 'x'); hold on
        title ('trajectory Part C');
        olocs = locs;
    
end

for n=1:ndx
    for k=1:ndy
        ni=locs(:,2)<(delta*k) &...
                 locs(:,2)>(delta*(k-1)) & locs(:,1)<(delta*n)...
                 & locs(:,1)>(delta*(n-1)); 
             density(n,k)=sum(ni);
             tempMap(n,k) = mean(((1/(2*K))*m*(vels(ni,1).^2 + vels(ni,2).^2)));       
             plot(locs(ni,1),locs(ni,2),'x');hold on
    end            
end

figure
surf(density');pause(0.001)
title('Electron Density Part C');
grid on

figure
surf(tempMap');pause(0.001)
title('Temperature Map Part C');
colorbar
grid on
