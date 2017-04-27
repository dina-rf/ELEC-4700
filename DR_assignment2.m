%%
Lx= 1.5; %length of x
Wy=1; %length of y
nx=150; %number of steps in the x-direction
ny=100; %number of steps in the y-direction
deltaX=Lx/(nx-1); %interval between steps
deltaY=Wy/(ny-1);
V0=1; %inital value of V0

%Part 1A, 1-D case where dV/dy = 0
V=zeros(nx,ny); %inital conditions
V(1,:)=V0;
V(nx,:)=0;
e=1; %max error
while(e>=0.00001) %keeps looping until enough measurements are taken and the result is accurate enough to be lower than 0.1
    V_old=V;
    for i=2:1:nx-1
        for j = 1:ny
            V(i,j) = (V(i+1,j)+V(i-1,j))/2;
        end
    end
    e = max(max(abs(V - V_old)));
%     surface(V','linestyle','none');pause(0.01)
end
X = linspace(0,Lx,nx);
figure
surface(V','linestyle','none');pause(0.01)
title('1D solution of Laplace equation');
ylabel('Voltage');
xlabel('Position');
grid on

%-------------
%%
%Part 1B, dV/dy is included into the computations

V=zeros(nx,ny); %inital conditions
V(1,:)=V0;
V(nx,:)=V0;
V(:,1)=0;
V(:,ny)=0;
e=1; %max error
while(e>=0.000001) 
    V_old=V;
    for i=2:1:nx-1
        for j=2:1:ny-1
            V(i,j) = (V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1))/4;
        end
    end
    e = max(max(abs(V - V_old)));
end
figure 
surf(V')   
ylabel('Y Position');
xlabel('X Position');
zlabel('Voltage');
title('Matching surface plot');
grid on

%-------------
%%
%solving it using the analytical situation
figure
V = zeros(nx,ny);
b = Lx/2;
a = Wy;
x = linspace(-Lx/2,Lx/2,nx);
y = linspace(0,Wy,ny);

for n=1:2:113 
    for i=2:1:nx-1
        xp = x(i);
        for j=2:1:ny-1           
            yp = y(j);
            addterm = (cosh(n*pi*xp/a)*sin(n*pi*yp/a)/(n*cosh(n*pi*b/a)));
            V(i,j) = V(i,j) + addterm;
        end
    end
    V1 = V.*(4*V0/pi); 
    
    surf(V1');pause(0.001)
    title('Matching surface plot - Analytical');
    ylabel('Y Position');
    xlabel('X Position');
    zlabel('Voltage');
    grid on
end
%%
% The analytical solution approach uses a formula and computes for a certain 
% value, however to get an accurate plot (that has perfect square edges unlike 
% the one we get which is smooth) an infinte number of data points is
% reuqired. The numerical solution approach used earlier make
% approximations while solving and ploting the waverforms, hence the
% waveform results with a lot of ripples.

%Part 2 - Effect on current flow due to the "bottle-neck"

x = linspace(0,Lx,nx);
y = linspace(0,Wy,ny);
sig = zeros(nx,ny); %conductivity matrix
V = zeros(nx,ny);

for i=1:1:nx
    for n=1:1:ny
        if((x(i)>=Lx/3 && x(i)<=2*Lx/3) && (y(n)<=3*Wy/8 || y(n)>=Wy-3*Wy/8))
            sig(i,n)=1e-2;
        else
            sig(i,n)=1;
        end
    end
end

e=1; %max error
while(e>=1e-6) 
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
                
                if (V(i,j) > 1)
                    pause(0.1)
                end
            end
        end
    end       
%     surf(V');pause(0.001)
    e = max(max(abs(V - V_old)));
end

%Plotting Sigma
figure
surf(sig');pause(0.001)
title('Resistivity plot');
ylabel('Y Position');
xlabel('X Position');
zlabel('sigma');
grid on

%Plotting Potential
figure
surf(V');pause(0.001)
title('Voltage plot');
ylabel('Y Position');
xlabel('X Position');
zlabel('Voltage');
grid on

%Calculating and plotting Electric field
[EyNg,ExNg] = gradient(V);			
Ex = -ExNg;							
Ey = -EyNg;

figure
surf(Ex');pause(0.001)
title('Electric field plot in the x-direction');
xlabel('X Position');
zlabel('Electrical field');
grid on

figure
surf(Ey');pause(0.001)
title('Electric field plot in the y-direction');
ylabel('Y Position');
zlabel('Electrical field');
grid on

%Calculating and plotting current density

Jx=sig.*Ex;

Jy=sig.*Ey;

figure
quiver(Jx',Jy');
pause(0.001)
title('Current Density plot');
ylabel('Y Position');
xlabel('X Position');
zlabel('Current Density');
grid on
    

    
    
    
    