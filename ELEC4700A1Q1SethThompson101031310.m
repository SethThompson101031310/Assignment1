% Name: Seth Thompson
% Student Number: 101031310

close all
clear
clc

% ELEC 4700 Assignment 1 | Part 1: Electron Modelling

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

% Assigning a particle a fixed velocity but with a random direction.
theta = (2*pi).*rand(electronAmount,2);
particleVelocityXDirection = thermalVel * cos(theta);
particleVelocityYDirection = thermalVel * sin(theta);

% First column is the initial velocity, second column is the next velocity,
% think vn and vn+1
particleVelocityXDirection(:,1) = particleVelocityXDirection(:,2);
particleVelocityYDirection(:,1) = particleVelocityYDirection(:,2);

% Timestep calculated from sqrt(100nm^2 + 200nm^2)/1000, which is approx
% 0.22nm, and then 0.2nm/thermalVel
timeStep = 0.22e-9/thermalVel; % seconds

% Creating a loop that will calculate each particles position as time goes
% on. Also making a counter that will keep track of the amount of
% iterations done. 
counter = 0;

% Calculating the displacement of the electron for each time step
particleXDisplacement = timeStep*particleVelocityXDirection(:,1);
particleYDisplacement = timeStep*particleVelocityYDirection(:,1);

% To calculate the semiconductor temperature, take the average of the
% kinetic energy of each particle.

% Initializing the kinetic energy as 0 along with celcius and kelvin temps
KineticE = 0;

% Creating a vector of colours for each particle being simulated
colourVector = rand(electronAmount,3);

for m = 1:1000
    
    % Updating the position of the particle
    for n = 1:electronAmount
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
    KineticE = 0.5*(mEffective*(particleVelocityXDirection.^2 + particleVelocityYDirection.^2));
    KineticE = sum(KineticE(:,1));

    % Calculating the average kinetic energy
    KineticEAVG(m) = KineticE/electronAmount;

    % Calculating the temperature of the semiconductor using KEAVG = 3/2*(k*T)
    SemiTempKelvin(m) = KineticEAVG(m)/kb;
    SemiTempCelcius(m) = SemiTempKelvin(m) - 273;
    
    % Creating a vector that will hold all timesteps over the simulation
    timeVector(m) = timeStep*m;
    
    if (counter == 0)
        % Plotting the particles over time.
        figure(1)
        scatter(particleXPosition(:,2),particleYPosition(:,2),1,colourVector(:,1))
        hold on
        title({['2-D Plot of Particle Trajectories | Semi Temp = ' num2str(SemiTempCelcius(m))],['Seth Thompson | 101031310']})
        xlabel('X-Axis (m)')
        ylabel('Y-Axis (m)')
        xlim([0 200e-9])
        ylim([0 100e-9])
    elseif (counter > 0 && counter < 1000)
        title({['2-D Plot of Particle Trajectories | Semi Temp = ' num2str(SemiTempCelcius(m))],['Seth Thompson | 101031310']})
        scatter(particleXPosition(:,2),particleYPosition(:,2),1,colourVector(:,1))
    else
        title({['2-D Plot of Particle Trajectories | Semi Temp = ' num2str(SemiTempCelcius(m))],['Seth Thompson | 101031310']})
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
figure(2)
plot(timeVector,SemiTempCelcius)
title({['Semiconductor Temperature over time'],['Seth Thompson | 101031310']})
xlabel('Time (s)')
ylabel('Temperature (C)')
