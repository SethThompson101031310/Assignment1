% Name: Seth Thompson
% Student Number: 101031310

close all
clear
clc

% ELEC 4700 Assignment 1 | Part 2: Collisions with Mean Free Path

% Defining Constants to be used in this part of the assignment.
mRest = 9.109e-31; % kilograms
mEffective = 0.26*mRest; % kilograms
regionLength = 200e-9; % meters
regionWidth = 100e-9; % meters

%-------------------------------------------------------------------------
% Question 1: What is the thermal velocity? Assume T = 300K.
%-------------------------------------------------------------------------

% Defining the given constant
Temperature = 300; % Kelvin
kb = 1.380649e-23; % J*K^-1

% Equation for thermal velocity is sqrt((2*kb*T)/m)(RMS)

thermalVel = sqrt((2*kb*Temperature)/mEffective);
thermalVeldisplay = thermalVel * 1e-3;
% Outputting the thermal velocity to the command window
fprintf('The thermal velocity, assuming T = 300K, is %f km/s. \n', thermalVeldisplay)

%-------------------------------------------------------------------------
% Question 2: If the mean time between collisions is Tmn = 0.2ps, what is
% the mean free path?
%-------------------------------------------------------------------------

% Defining the given constant
meanTime = 0.2e-12; % seconds

% Mean free path is defined as the average distance travelled by a moving
% particle (in this case an electron), calculating this is as simple as d =
% v*t

MFP = thermalVel*meanTime;
MFPdisplay = MFP *1e9;

% Outputting the MFP to the command window.

fprintf('The mean free path, assuming a mean time of 0.2ps, is %f nm.\n', MFPdisplay)
%-------------------------------------------------------------------------
% Question 3: Write a program that will model the random motion of the
% electrons.
%-------------------------------------------------------------------------

% Assigning each particle a random starting position in the 100x200 plane.

electronAmount = 20; % Amount of electrons being simulated

% Assigning a particle a random position on the XY plane within the maximum
% and minimum width and length.
particleXPosition = regionLength.*rand(electronAmount,2);
particleYPosition = regionWidth.*rand(electronAmount,2);

% First column is the initial position, second column is the next position,
% think xn and xn+1
particleXPosition(:,1) = particleXPosition(:,2);
particleYPosition(:,1) = particleYPosition(:,2);

%-------------------------------------------------------------------------
% Question 3: Write a program that will model the random motion of the
% electrons.
%-------------------------------------------------------------------------
% For part 2 of the assignment, each particle must be assigned a random
% velocity using the maxwell boltzman distribution. This is shown next.

% Assigning a particle a random velocity using the randn function (normally distributed random values).

randXV = randn(electronAmount,2)*sqrt((kb*Temperature)/mEffective);
randYV = randn(electronAmount,2)*sqrt((kb*Temperature)/mEffective);

% Calculating all random initial velocities to plot a histogram.

randV = (randXV(:,1).^2 + randYV(:,1).^2).^(0.5);
figure(1)
hist(randV,electronAmount)
title({['Histogram of Particle Speeds in Simulation'],['Seth Thompson | 101031310']})
xlabel('Speed of Particle (m/s)')
ylabel('Number of Occurences')
% Timestep calculated from sqrt(100nm^2 + 200nm^2)/1000, which is approx
% 0.22nm, and then 0.2nm/thermalVel
timeStep = 0.22e-9/thermalVel; % seconds

pause(4)

% Creating a loop that will calculate each particles position as time goes
% on. Also making a counter that will keep track of the amount of
% iterations done. 
counter = 0;

% Calculating the displacement of the electron for each time step
particleXDisplacement = timeStep*randXV(:,1);
particleYDisplacement = timeStep*randYV(:,1);

% To calculate the semiconductor temperature, take the average of the
% kinetic energy of each particle.

% Initializing the kinetic energy as 0.
KineticE = 0;

% Creating a vector of random RGB values so each particle will have its own
% colour
colourVector = rand(electronAmount,3);

% Creating avariable to hold the probabillity of an electron scattering
PScat = 1 - exp(-(timeStep/meanTime));

% Creating two vectors that will hold the initial position before a particle
% has deflected. One for the X position and another for the Y. There will
% also be a collide counter to determine the aount of collisions made
% during the system

beforeCollideX = particleXPosition(:,1);
beforeCollideY = particleYPosition(:,1);

% Making another variable that will hold the timestep before the next collision
timeBeforeCollision = 0;

% Creating a collision counter variable to track the amount of collisions
collisionCount = 0;

% pre-defining semitempkelvin and semitempcelcius
SemiTempKelvin = zeros(1,10000);
SemiTempCelcius = zeros(1,10000);
for m = 1:1000
    
    % Updating the position of the particle
    for n = 1:electronAmount
        % First checking to see if the electron scatters, if it does then a
        % new velocity must be generated and from there a new displacement
        if (rand(1) < PScat)
            % Increment collisionCount
            collisionCount = collisionCount + 1;
            thetaScatter = (2*pi).*rand(1);
            randXV(n,:) = randn(1)*sqrt((kb*Temperature)/mEffective);
            randYV(n,:) = randn(1)*sqrt((kb*Temperature)/mEffective);
            particleXDisplacement(n,:) = timeStep*randXV(n,1);
            particleYDisplacement(n,:) = timeStep*randYV(n,1);
            
            % Calculating the time it took for a particle to deflected and
            % the distance it travelld before doing so. These results will
            % be stored in a vector.
            distanceX = particleXPosition(n,2) - beforeCollideX(n);
            distanceY = particleYPosition(n,2) - beforeCollideY(n);
            distance(collisionCount) = sqrt(distanceX^2 + distanceY^2);
            
            % Calculating and saving time increments between collisions as
            % well.
            
            timeTook(collisionCount) = abs(timeBeforeCollision - timeStep*n);
            timeBeforeCollision = timeStep*n;
            % Saving the new positions before the next collision
            beforeCollideX = particleXPosition(:,1);
            beforeCollideY = particleYPosition(:,1);
        end
        % Checking the boundary conditions for the left and right sides of
        % the plot. If the particle is going to pass a border on the left
        % or right side then simply move it to the left or right.
        if (particleXPosition(n,1) + particleXDisplacement(n) > 2e-7)
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1) - 2e-7;
        elseif (particleXPosition(n,1) + particleXDisplacement(n) < 0)
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1) + 2e-7;
        else
            particleXPosition(n,2) = particleXDisplacement(n) + particleXPosition(n,1);
        end
        
        % Checking boundary conditions for the top and bottom of the plot.
        % If the particle is about to cross the top of the plot, flip the
        % y, displacement value and have it move in the opposite direction;
        % the opposite applies for the bottom of the plot.
        if ((particleYPosition(n,1) + particleYDisplacement(n) > 1e-7) || (particleYPosition(n,1) + particleYDisplacement(n) < 0))
            particleYDisplacement(n) = -particleYDisplacement(n);
            particleYPosition(n,2) = particleYDisplacement(n) + particleYPosition(n,1);
        else
            particleYPosition(n,2) = particleYDisplacement(n) + particleYPosition(n,1);
        end
    end
    
    % Calculating the sum of the kinetic energys
    KineticE = 0.5*(mEffective*(randXV.^2 + randYV.^2));
    KineticE = sum(KineticE(:,1));

    % Calculating the average kinetic energy
    KineticEAVG(m) = KineticE;

    % Calculating the temperature of the semiconductor using KEAVG = 3/2*(k*T)
    SemiTempKelvin(m) = SemiTempKelvin(m) + KineticEAVG(m)/kb;
    SemiTempCelcius(m) = SemiTempKelvin(m)  + (SemiTempCelcius(m) - 273);
    
    % Averaging out both temperatures to approach a final answer.
    SemiTempKelvinAVG(m) = SemiTempKelvin(m)/m;
    SemiTempCelciusAVG(m) = SemiTempCelcius(m)/m;
    
    % Creating a vector that will hold all timesteps over the simulation
    timeVector(m) = timeStep*m;
    
    if (counter == 0)
        % Plotting the particles over time.
        figure(2)
        scatter(particleXPosition(:,2),particleYPosition(:,2),1,colourVector(:,1))
        hold on
        title({['2-D Plot of Particle Trajectories | Semi Temp = ' num2str(SemiTempCelciusAVG(m))],['Seth Thompson | 101031310']})
        xlabel('X-Axis (m)')
        ylabel('Y-Axis (m)')
        xlim([0 200e-9])
        ylim([0 100e-9])
    elseif (counter > 0 && counter < 1000)
        title({['2-D Plot of Particle Trajectories | Semi Temp = ' num2str(SemiTempCelciusAVG(m))],['Seth Thompson | 101031310']})
        scatter(particleXPosition(:,2),particleYPosition(:,2),1,colourVector(:,1))
        hold on
    else
        title({['2-D Plot of Particle Trajectories | Semi Temp = ' num2str(SemiTempCelciusAVG(m))],['Seth Thompson | 101031310']})
        scatter(particleXPosition(:,2),particleYPosition(:,2),1,colourVector(:,1))
        hold off
    end
    
    pause(0.001)
    
    % Updating the first column of each position vector so the position of
    % each particle can be incremented
    particleXPosition(:,1) = particleXPosition(:,2);
    particleYPosition(:,1) = particleYPosition(:,2);
    
    counter = counter + 1;
end

% Plotting the temperature of the simulation over time
figure(3)
plot(timeVector,SemiTempCelciusAVG)
title({['Semiconductor Temperature over time'],['Seth Thompson | 101031310']})
xlabel('Time (s)')
ylabel('Temperature (C)')
ylim([-10,100])

% Calculating and outputting the calculated MFP
distanceSUM = sum(distance);
MFPCalc = distanceSUM/collisionCount;
MFPCalcDisp = MFPCalc*1e9;
fprintf('The calculated MFP of the simulation is %f nm.\n', MFPCalcDisp)
% Calculating and outputting the mean time between collisions
timeSUM = sum(timeTook);
meanTimeCalc = timeSUM/collisionCount;
meanTimeCalcDisp = meanTimeCalc * 1e12;
fprintf('The calculated mean time between collisions is %f fs.\n', meanTimeCalcDisp)