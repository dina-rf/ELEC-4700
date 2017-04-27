
close all
clc
% PART 1
clear temp

%Constants
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

%Creating boxes
specular=1; % 0 means diffusive

xm=0.8e-7;
xp=1.2e-7;
yl=0;
yMlw=0.4e-7;
yMup=0.6e-7;
yUp=1e-7;

xline=[xm xp xp xm xm];
topybox=[yl yl yMlw yMlw yl];
bottomybox=[yMup yMup yUp yUp yMup];

figure
plot(xline, topybox,'b-', xline, bottomybox,'b-')
axis([0 w 0 l])
hold on

for j = 1:timesteps
      
        if pscat>rand(j)
            vels(:,1)=randn(npar,1)*vth;
            vels(:,2)=randn(npar,1)*vth;
        end
        %updating the particle location using Newton's law of motion
        locs(:,1) = locs(:,1)+ vels(:,1)*deltaT;
        locs(:,2) = locs(:,2)+ vels(:,2)*deltaT;

        % check to what inside
        ix = locs(:,1) > xm & locs(:,1) < xp;
        iy = locs(:,2) > yMup | locs(:,2) < yMlw;
        it = ix & iy;

        for k=1:length(it)

            if (it(k)==1 && specular ==1)
                olocs(k,1)=locs(k,1)-(vels(k,1)*deltaT);
                olocs(k,2)=locs(k,2)-(vels(k,2)*deltaT);
                if (olocs(k,1)<xm && locs(k,1)>xm) || (locs(k,1)<xp && olocs(k,1)>xp)   
                    vels(k,1)=-vels(k,1);
                elseif ((olocs(k,2)<yMup && olocs(k,2)>yMlw) && (locs(k,2)>yMup || locs(k,2)<yMlw))
                    vels(k,2)=-vels(k,2);
                end
                locs(k,1) = locs(k,1)+ vels(k,1)*deltaT;
                locs(k,2) = locs(k,2)+ vels(k,2)*deltaT;
             end

%             if (it(k)==1 && specular ==0)
%                 olocs(:,1)=locs(:,1)-(vels(:,1)*deltaT);
%                 olocs(:,2)=locs(:,2)-(vels(:,2)*deltaT);
%                 if (olocs(k,1)<xm && locs(k,1)>xm) || (locs(k,1)<xp && olocs(k,1)>xp)   
%                     vels(k,1)=randn*vth;
%                 elseif ((olocs(k,2)<yMlw && olocs(k,2)>yMlw) && (locs(k,2)>yMup || locs(k,2)<yMlw))
%                     vels(j,2)=randn*vth;
%                 end
%                 locs(:,1) = locs(:,1)+ vels(:,1)*deltaT;
%                 locs(:,2) = locs(:,2)+ vels(:,2)*deltaT;
%             end
            X=[olocs(k,1) locs(k,1)];
            Y=[olocs(k,2) locs(k,2)];
            plot (X,Y, 'x');
       end

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

        temp(:,j) = mean((1/(2*K))*m*((vels(:,1).^2)+(vels(:,2).^2)));
       
        X=[olocs(:,1) locs(:,1)];
        Y=[olocs(:,2) locs(:,2)];
        pause(0.00001)
        plot (X,Y, 'x');
        title ('trajectory');
        olocs = locs;

end

%calculating and plotting the temperature of the semiconductor

figure 
plot(1:timesteps,temp); 
title('Semiconductor Temperature over Time')
xlabel('Time')
ylabel('Temperature')
grid on;
 
    
%average of the speeds will be vth using a histogram
%function
vavg= sqrt((vels(:,1).^2) + (vels(:,2).^2));
figure
hist(vavg);

%density
delta = 1e-8;
ndx = w/delta;
ndy = l/delta;
density=zeros(ndx,ndy);
tempMap = zeros(ndx,ndy);

figure
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
title('Electron Density');
grid on

figure
surf(tempMap');pause(0.001)
title('Temperature Map');
colorbar
grid on



