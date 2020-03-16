%{
ELEC 4700 Assignment 3
William Fraser
101001393
%}

colours = ['b','g','r','c','m','y','k'];
maxTraces = length(colours);
numTraces = maxTraces;
numparticles = 30000;

density = 10^15; %cm^-2
kb = 1.38064852e-23;
m0 = 9.11e-31;
m = 0.26*m0;
T = 300;
vth = sqrt((2*kb*T)/m);
vth = vth*1e9;
vth = vth*1e-12; %nm/ps
vth = vth/100; %Timestep = 10fs
vth2 = vth^2;

eField = 0.01/200e-9; %V/m;
q = 1.602e-19;
eForce = q*eField; %N
accel = eForce/m0; %m/s^2
accel = accel * 1e9 * ((1e-12)^2); %nm/ps^2
accel = accel/100; %Timestep = 10fs

grid1 = zeros(100,200);
particles = zeros(numparticles,6); %Linear Index; row; column; X velocity; Y velocity; X Accel
particles(:,6) = accel;

particles(1:100,2) = randperm(100,100);
particles(1:100,3) = randperm(200,100);

start = 1;
while(start <= 29901)
    stop = start + 99;
    particles(start:stop,2) = randperm(100,100);
    particles(start:stop,3) = randperm(200,100);
    start = start + 100;
end

particles(:,1) = sub2ind(size(grid1),particles(:,2),particles(:,3));
tracesY = zeros(numTraces,2);
tracesY(:,1) = particles(1:numTraces,2);
tracesX = zeros(numTraces,2);
tracesX(:,1) = particles(1:numTraces,3);
J = zeros(1000,1);

for i = 1:length(particles(:,1))
    xRat = rand;
    xDir = (-1)^(round(rand));
    yDir = (-1)^(round(rand));
    yRat = 1 - xRat;
    xVel = (sqrt(xRat*vth2));
    yVel = (sqrt(yRat*vth2));
    particles(i,4) = xVel*xDir;
    particles(i,5) = yVel*yDir;
end

grid1(particles(:,1)) = 1;
gridSize = size(grid1);
kbMax = kb*1e18*1e-24; %Fixing units
tempVector = zeros(1,1000);

figure(1)
xlim([0 200])
ylim([0 100])
set(gca, 'YDir','reverse')
title('Electron Path Tracing')
xlabel('X Position (nm)')
ylabel('Y Position (nm)')
hold on

for time = 1:1000
    squaredVel = (particles(:,4).^2) + (particles(:,5).^2);
    meanVel = mean(squaredVel);
    temperature = (m*meanVel)/kbMax;
    tempVector(time) = temperature;
    
    fixUnitsVel = mean(particles(:,4))*(1e-7)*(1e12); %cm/s;
    eField = eField/100; %V/cm
    miu = fixUnitsVel/eField; %cm^2/Vs
    J(time) = fixUnitsVel*q*density*30000; %A/cm^2
 
    display('Time = ',num2str(time));
    if(time~=1)
        tracesX(:,1) = tracesX(:,2);
        tracesY(:,1) = tracesY(:,2);
    end
    particles(:,3) = particles(:,3) + particles(:,4);
    particles(:,2) = particles(:,2) + particles(:,5);
    particles(:,4) = particles(:,4) + particles(:,6);
    yBoundMaxed = particles(:,2) > gridSize(1);
    yBoundMined = particles(:,2) < 0;
    yIsZero = particles(:,2)==0;
    yBounds = yBoundMaxed + yBoundMined + yIsZero;
    keepY = yBounds==0;
    particles(:,2) = (particles(:,2).*keepY) - (yBoundMined.*particles(:,2)) + (yBoundMaxed.*(gridSize(1)-(particles(:,2)-gridSize(1)))) - round(yIsZero.*particles(:,5));
    particles(:,5) = (particles(:,5).*keepY) - (yBoundMined.*particles(:,5)) - (yBoundMaxed.*particles(:,5)) - (particles(:,5).*yIsZero);
    xBoundMaxed = particles(:,3) > gridSize(2);
    xBoundMined = particles(:,3) < 0;
    xIsZero = particles(:,3)==0;
    xBounds = xBoundMaxed + xBoundMined + xIsZero;
    keepX = xBounds==0;
    particles(:,3) = (particles(:,3).*keepX) + (xBoundMaxed.*(particles(:,3)-gridSize(2))) + (xBoundMined.*(particles(:,3)+gridSize(2))) + (xIsZero.*round(gridSize(2)+particles(:,4)));
    fix1 = floor(particles(:,2)) == 0;
    fix2 = floor(particles(:,3)) == 0;
    particles(:,2) = particles(:,2) + fix1;
    particles(:,3) = particles(:,3) + fix2;
    particles(:,1) = sub2ind(gridSize,floor(particles(:,2)),floor(particles(:,3)));
    tracesX(:,2) = particles(1:numTraces,3);
    tracesY(:,2) = particles(1:numTraces,2);
    holder = tracesX(:,2);
    jump = abs(tracesX(:,1) - tracesX(:,2)) > 100;
    noJump = jump == 0;
    temp = tracesX(:,1).*noJump;
    tracesX(:,1) = temp + (jump.*tracesX(:,2));
    holder = tracesY(:,2);
    jump = abs(tracesY(:,1) - tracesY(:,2)) > 70;
    noJump = jump == 0;
    temp = tracesY(:,1).*noJump;
    tracesY(:,1) = temp + (jump.*tracesY(:,2));
    %{
    if(time~=1)
        x(1,:) = tracesX(:,1);
        x(2,:) = tracesX(:,2);
        y(1,:) = tracesY(:,1);
        y(2,:) = tracesY(:,2);

        plot(x(:,1),y(:,1),colours(1),x(:,2),y(:,2),colours(2),x(:,3),y(:,3),colours(3),x(:,4),y(:,4),colours(4),x(:,5),y(:,5),colours(5),x(:,6),y(:,6),colours(6),x(:,7),y(:,7),colours(7))
    end
    %}
    pause(0.1)
end

timeVector = zeros(1,1000);
for i = 1:length(timeVector)
    timeVector(i) = i/100;
end

J = J * 1e6 * (1e-7) * (1e-7); %uA/nm^2;

figure(2)
plot(timeVector,J,'b')
title('Current Density Over Time')
xlabel('Time (ps)')
ylabel('Current Density (uA/nm^2)')

figure(3)
plot(timeVector,tempVector,'b')
title('Material Temperature Over Time')
xlabel('Time (ps)')
ylabel('Temperature (K)')

densityGrid = zeros(100,200);
densityGrid(particles(:,1)) = 1;
%h = ones(10,10);
%densityGrid = imfilter(densityGrid,h);
figure(4)
imagesc(densityGrid)
colorbar
colormap cool
title('Electron Position Density (Electrons/nm^2)')

electronTemperatures = (squaredVel.*m)./kbMax;
tempGrid = zeros(100,200);
tempGrid(particles(:,1)) = electronTemperatures;
%tempGrid = imfilter(tempGrid,h);
figure(5)
imagesc(tempGrid)
colorbar
colormap cool
title('Temperature Map (K)')
