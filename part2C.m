%{
ELEC 4700 Assignment 3
William Fraser
101001393
%}
load Evector
load x
fieldX = x;
clear x
load y
fieldY = y;
clear y;

for i = 1:length(fieldX)
    eFieldX(fieldY(i),fieldX(i)) = Evector(i,1);
    eFieldY(fieldY(i),fieldX(i)) = Evector(i,2);
end    

eFieldX = eFieldX * 1e9; %V/m
eFieldY = eFieldY * 1e9; %V/m

colours = ['b','g','r','c','m','y','k'];
maxTraces = length(colours);
numparticles = 1000;

numTraces = maxTraces;

kb = 1.38064852e-23;
m0 = 9.11e-31;
m = 0.26*m0;
T = 300;
vth = sqrt((2*kb*T)/m);
vth = vth*1e9;
vth = vth*1e-12; %nm/ps
vth = vth/100; %Timestep = 10fs
vth2 = vth^2;

q = 1.602e-19;
eForceX = q*eFieldX; %N
accelX = eForceX/m0; %m/s^2
accelX = accelX * 1e9 * ((1e-12)^2); %nm/ps^2
accelX = accelX/100; %Timestep = 10fs
eForceY = q*eFieldY; %N
accelY = eForceY/m0; %m/s^2
accelY = accelY * 1e9 * ((1e-12)^2); %nm/ps^2
accelY = accelY/100; %Timestep = 10f

velocities = linspace(0,50,102);
v2 = velocities.^2;
coef1 = 4*pi;
kbMax = kb*1e18*1e-26; %Fixing units
coef2 = (m0/(2*pi*kbMax*T));
exponent = exp(-(m0*v2)/(2*kbMax*T));
maxBol = coef1*coef2*v2.*exponent*numparticles;
velAssigner = zeros(1,round(sum(maxBol)));
currentIndex = 1;

for i = 1:length(maxBol)
    if(round(maxBol(i))~=0)
        stopPoint = currentIndex + round(maxBol(i))-1;
        velAssigner(currentIndex:stopPoint) = velocities(i);
    else
        stopPoint = currentIndex;
    end
    currentIndex = stopPoint + 1;
end

velAssigner = velAssigner(randperm(length(velAssigner)));

grid1 = zeros(100,200);
particles = zeros(numparticles,7); %Linear Index; row; column; X velocity; Y velocity; X Accel; Y Accel
if(numparticles <= 100)
    particles(:,2) = randperm(100,numparticles);
    particles(:,3) = randperm(200,numparticles);
else
    start = 1;
    stop = 100;
    endReached = 0;
    while(endReached==0)
        if(stop<=numparticles)
            particles(start:stop,2) = randperm(100,100);
            particles(start:stop,3) = randperm(200,100);
            start = start + 100;
            stop = stop + 100;
        elseif(stop>=numparticles)
            particles(start:numparticles,2) = randperm(100,(numparticles-start+1));
            particles(start:numparticles,3) = randperm(100,(numparticles-start+1));
            endReached = 1;
        end 
    end
end
inBoxX1 = particles(:,3) >= 75;
inBoxX2 = particles(:,3) <= 125;
inBoxX = inBoxX1.*inBoxX2;
inBoxY1 = particles(:,2) <= 40;
inBoxY2 = particles(:,2) >= 60;
inBox1 = inBoxY1.*inBoxX;
inBox2 = inBoxY2.*inBoxX;
inBox = inBox1 + inBox2;
notInBox = inBox==0;
particles(:,2) = (particles(:,2).*notInBox) + (inBox*50);
particles(:,3) = (particles(:,3).*notInBox) + (inBox*100);
particles(:,1) = sub2ind(size(grid1),particles(:,2),particles(:,3));
tracesY = zeros(numTraces,2);
tracesY(:,1) = particles(1:numTraces,2);
tracesX = zeros(numTraces,2);
tracesX(:,1) = particles(1:numTraces,3);
maxIndex = length(velAssigner);

tempInd = sub2ind(size(accelX),particles(:,2),particles(:,3));
particles(:,6) = accelX(tempInd);
particles(:,7) = accelY(tempInd);

for i = 1:length(particles(:,1))
    choice = randperm(maxIndex,1);
    tempVel = velAssigner(choice);
    tempVel2 = tempVel^2;
    xRat = rand;
    xDir = (-1)^(round(rand));
    yDir = (-1)^(round(rand));
    yRat = 1 - xRat;
    xVel = (sqrt(xRat*tempVel2));
    yVel = (sqrt(yRat*tempVel2));
    particles(i,4) = xVel*xDir;
    particles(i,5) = yVel*yDir;
end

grid1(particles(:,1)) = 1;
gridSize = size(grid1);

pScat = 1 - exp(-0.1/0.2);

clear choice coef1 coef2 currentIndex exponent i kb m0 maxBol normalize
clear stopPoint T tempVel tempVel2 v2 vth vth2 velocities histVel maxBolNorm

scatterMatrix = zeros(numparticles,1000);
averageVel = zeros(1,1000);

boxX = [75 75 125 125];
box1y = [0 40 40 0];
box2y = [100 60 60 100];

figure(1)
title('Electron Path Tracing')
ylabel('Y Position (nm)')
xlabel('X Position (nm)')
set(gca, 'YDir','reverse')
plot(boxX,box1y,'k',boxX,box2y,'k')
xlim([0 200])
ylim([0 100])
hold on

for time = 1:1000
    squaredVel = (particles(:,4).^2) + (particles(:,5).^2);
    meanVel = mean(squaredVel);
    averageVel(time) = sqrt(meanVel);
    temperature = (m*meanVel)/kbMax;
    %display('Time = ',num2str(time));
    if(time~=1)
        tracesX(:,1) = tracesX(:,2);
        tracesY(:,1) = tracesY(:,2);
    end
    grid2 = zeros(size(grid1));
    inBoxX1 = particles(:,3) >= 75;
    inBoxX2 = particles(:,3) <= 125;
    inBoxX = inBoxX1.*inBoxX2;
    inBoxY1 = particles(:,2) <= 40;
    inBoxY2 = particles(:,2) >= 60;
    inBoxY = inBoxY1.*inBoxY2;
    enter1Left = particles(:,3)+particles(:,4)>=75;
    enter1Left = enter1Left.*inBoxX2.*inBoxY1;
    enter1Right = particles(:,3)+particles(:,4)<=125;
    enter1Right = enter1Right.*inBoxX1.*inBoxY1;
    enter1Sides = enter1Right + enter1Left;
    not1Sides = enter1Sides==0;
    particles(:,4) = (particles(:,4).*not1Sides) - (particles(:,4).*enter1Sides);
    enter2Left = particles(:,3)+particles(:,4)>=75;
    enter2Left = enter2Left.*inBoxX2.*inBoxY2;
    enter2Right = particles(:,3)+particles(:,4)<=125;
    enter2Right = enter2Right.*inBoxX1.*inBoxY2;
    enter2Sides = enter2Right + enter2Left;
    not2Sides = enter2Sides==0;
    particles(:,4) = (particles(:,4).*not2Sides) - (particles(:,4).*enter2Sides);
    enter1Bottom = particles(:,2)+particles(:,5)<=40;
    enter1Bottom = enter1Bottom.*inBoxX1.*inBoxX2;
    enter2Top = particles(:,2)+particles(:,5)>=60;
    enter2Top = enter2Top.*inBoxX1.*inBoxX2;
    enterHorizontal = enter1Bottom+enter2Top;
    noEnterHorizontal = enterHorizontal==0;
    particles(:,5) = (particles(:,5).*noEnterHorizontal) - (particles(:,5).*enterHorizontal);
    particles(:,3) = round(particles(:,3) + particles(:,4));
    particles(:,2) = round(particles(:,2) + particles(:,5));
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
    particles(:,1) = sub2ind(gridSize,particles(:,2),particles(:,3));
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
    if(time~=1)
        x(1,:) = tracesX(:,1);
        x(2,:) = tracesX(:,2);
        y(1,:) = tracesY(:,1);
        y(2,:) = tracesY(:,2);
        plot(x(:,1),y(:,1),colours(1),x(:,2),y(:,2),colours(2),x(:,3),y(:,3),colours(3),x(:,4),y(:,4),colours(4),x(:,5),y(:,5),colours(5),x(:,6),y(:,6),colours(6),x(:,7),y(:,7),colours(7))
        pause(0.1)
    end
    particles(:,4) = particles(:,4) + particles(:,6);
    particles(:,5) = particles(:,5) + particles(:,7);
    tempInd = sub2ind(size(accelX),particles(:,2),particles(:,3));
    particles(:,6) = accelX(tempInd);
    particles(:,7) = accelY(tempInd);
end

timeVector = zeros(1,1000);
for i = 1:length(timeVector)
    timeVector(i) = i/100;
end

densityGrid = zeros(100,200);
densityGrid(particles(:,1)) = 1;
figure(2)
imagesc(densityGrid)
colorbar
colormap cool
title('Electron Position Density (Electrons/nm^2)')