%% ENEL 627 Assignment 3
% Analysis of Van Atta Array using different lengths
% code built off of ENEL 627 assignment 
%
%% Van Atta array parameters
% Frequency 79 GHz
% Transmission line length of first pair = x
% Transmission line length of second pair = y
% Element spacing = z
% Assume Gain and Array factor = 1
% Relative Phase difference calculation
    % First start with same lengths
        % Each pair gets signal 
        % go through phase conjugation calculation
        % see that emitted wave is the same as received
    %
clear all 
close all
clc
% Define parameters of incident wave
fo = 79e9;
c = 2.99792458e8;
lamda = 1/fo;
theta = linspace(0,2*pi,500);
phi = linspace(0,2*pi,500);
d = 0.7*lamda;
incident_phi = 2*pi/lamda*d*sin(theta);

% define parameters of array
l1 = 100;                               % dimensions in mm
l2 = linspace(100-lamda,100+lamda,500);     % l2 over different lengths       
N = l2/(1000*lamda);
l2_electrical = 2*pi*N; 

A1_relative_phase = zeros(1,500);
A2_relative_phase = -1*incident_phi;
A3_relative_phase = -2*incident_phi;
A4_relative_phase = -3*incident_phi;

% use dipole element pattern
% array factor will be generated based on phase of input signal
% array factor is simply superposition of the different relative phases
% (assuming amplitudes of 1)
% so for actual van atta array
% AF_VA = Van_atta_phase(1) + Van_atta_phase(2) + Van_atta_phase(3) + Van_atta_phase(4);


AF_VA = exp(1i*A1_relative_phase) + exp(1i*A2_relative_phase) + exp(1i*A3_relative_phase) + exp(1i*A4_relative_phase);
AF_VA_mag = abs(AF_VA)/4;

for m = 1:500;
    for n = 1:500;
        ElementPattern_HW(m,n) = cos(pi/2*sin(theta(m))*cos(phi(n))) / sqrt((1 - (sin(theta(m))^2*cos(phi(n))^2)));
    end 
end;

% Note we have to plot the element factpr in both planes
figure;
polar(theta, ElementPattern_HW(1,:));
title('Figure 1: Half-Wave Dipole Element Pattern in Y-Z plane');

ElementPatternT = ElementPattern_HW';
figure;
polar(theta, ElementPatternT(1,:));
title('Figure 2: Half-Wave Dipole Element Pattern in X-Z plane');

% Plot Array Pattern
ArrayPattern_HW = zeros(1,500);
%ArrayPattern for X-Z
for m = 1:500;
    for n = 1:500;
        ArrayPattern_HW(m,n) = ElementPattern_HW(m,n) * AF_VA(n);
    end 
end;

%ArrayPattern for Y-Z
for m = 1:500;
    for n = 1:500;
        ArrayPattern_HW_YZ(m,n) = ElementPattern_HW(m,n) * AF_VA(m);
    end 
end;

figure;
polar(theta, ArrayPattern_HW(1,:));
title('Figure 4: Half-Wave Dipole Array Pattern Y-Z');

ArrayPatternT = ArrayPattern_HW_YZ';
figure;
polar(theta, ArrayPatternT(1,:));
title('Figure 5: Half-Wave Dipole Array Pattern X-Z');


%% Question 8.2-5
beta = 2 * pi;
d = 0.001;
Directivity = zeros(1,2000);

for n = 1:2000;
    
    Directivity(n) = 2 / (1 + (sin(beta * d * n) / (beta * d * n)));

end;

figure;
plot(Directivity);
title('Figure 1: Directivity of two-element Broadside Array, alpha = 0');
xlabel('d (function of Lamda)');
ylabel('Directivity');

%% Question 8.4-3
close all
% Plot the element pattern, array factor and array pattern 
% For half-wave dipole spaced lamda apart

% Plot Element Pattern
% create theta and phi array with values ranging 0 to 2pi
theta = pi/250:pi/250:2*pi;
phi = pi/250:pi/250:2*pi;

% Element pattern array of zeros, size 1000
ElementPattern_HW = zeros(1,500);

for m = 1:500;
    for n = 1:500;
        ElementPattern_HW(m,n) = cos(pi/2*sin(theta(m))*cos(phi(n))) / sqrt((1 - (sin(theta(m))^2*cos(phi(n))^2)));
    end 
end;

% Note we have to plot the element factpr in both planes

figure;
polar(theta, ElementPattern_HW(1,:));
title('Figure 1: Half-Wave Dipole Element Pattern in Y-Z plane');

ElementPatternT = ElementPattern_HW';
figure;
polar(theta, ElementPatternT(1,:));
title('Figure 2: Half-Wave Dipole Element Pattern in X-Z plane');

% Plot Array Factor
ArrayFactor_HW = zeros(1,500);

for m = 1:500;
    
    ArrayFactor_HW(m) = cos((pi) * cos(theta(m)));

end;

figure;
polar(theta, ArrayFactor_HW);
title('Figure 3: Half-Wave Dipole Array Factor');

% Plot Array Pattern
ArrayPattern_HW = zeros(1,500);

%ArrayPattern for X-Z
for m = 1:500;
    for n = 1:500;
        ArrayPattern_HW(m,n) = ElementPattern_HW(m,n) * ArrayFactor_HW(n);
    end 
end;

%ArrayPattern for Y-Z
for m = 1:500;
    for n = 1:500;
        ArrayPattern_HW_YZ(m,n) = ElementPattern_HW(m,n) * ArrayFactor_HW(m);
    end 
end;

figure;
polar(theta, ArrayPattern_HW(1,:));
title('Figure 4: Half-Wave Dipole Array Pattern Y-Z');

ArrayPatternT = ArrayPattern_HW_YZ';
figure;
polar(theta, ArrayPatternT(1,:));
title('Figure 5: Half-Wave Dipole Array Pattern X-Z');


%% Question 8.4-6
% Plot the array pattern for four small loop antennas 
% xz plane phi = 0
N = 4;
phi = 0;
% create theta array with values ranging 0 to 2pi
theta = pi/250:pi/250:2*pi;
phi = pi/250:pi/250:2*pi;
% Define psi using general formula and setup for ordinary end-fire
for m = 1:500;
    psi(m) = 0.8*pi*cos(theta(m)) + (0.8*pi);
end;
 
% Plot Element Pattern

ElementPattern_loop_xz = zeros(1,500);

for m = 1:500;
    ElementPattern_loop_xz(m) = sqrt(1 - (sin(theta(m))^2));
    
end;

figure;
polar(theta, ElementPattern_loop_xz);
title('Figure 1: Small Loop Array, Element Pattern');

% Plot Array Factor
% Array Factor array of zeros, size 1000
ArrayFactor_loop_xz = zeros(1,500);

for m = 1:500;
    %ArrayFactor_loop_xz(m) = sin(psi(m)*N/2) / (N*sin(psi(m)/2));
    ArrayFactor_loop_xz(m) = (sin(1.6*pi*cos(theta(m))-1.6*pi))/(4*sin(0.4*pi*cos(theta(m))-0.4*pi));
end;

figure;
AF_Loop_Mag = abs(ArrayFactor_loop_xz);
plot(psi, AF_Loop_Mag) 

figure;
polar(theta, AF_Loop_Mag);
title('Figure 2: Small Loop Array Factor, N = 4');

% Plot Array Pattern
% Array Pattern array of zeros, size 1000
ArrayPattern_loop_xz = zeros(1,500);

for m = 1:500;
    
    ArrayPattern_loop_xz(m) = ElementPattern_loop_xz(m) * ArrayFactor_loop_xz(m);

end;
APattern_Loop_Mag = abs(ArrayPattern_loop_xz);

figure;
polar(theta, APattern_Loop_Mag);
title('Figure 3: Small Loop Array Pattern, N = 4');
%% 8.4-8
% Holy shit it's getting late

for m = 1:500;
    psi(m) = 0.83*2*pi*cos(theta(m));
end;

for m = 1:500;
    ArrayFactorDipoleDesign(m) = abs(sin(2*psi(m))/(4*sin(psi(m)/2)));
end

figure;
polar(theta, ArrayFactorDipoleDesign);
title('Figure 1: Design With 4 Parallel Dipoles, N = 4');


for m = 1:500;
    ArrayPatternEPlane(m) = ArrayFactorDipoleDesign(m)*ElementPattern_loop_xz(m);
end

figure;
polar(theta, ArrayPatternEPlane);
title('Figure 2: Array Pattern in E plane, N = 4');

