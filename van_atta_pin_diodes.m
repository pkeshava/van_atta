%% ENEL 627 Assignment 3
% Analysis of Van Atta Array using different lengths
% code built off of ENEL 627 assignment 
%
%% Van Atta array parameters
% Frequency 79 GHz
% Transmission line length of first pair = l1
% Transmission line length of second pair = l2
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
d = 0.5*lamda;
incident_phi = 2*pi/lamda*d*cos(theta);

% Define parameters of 4 element array
A1_relative_phase = zeros(1,500);
A2_relative_phase = -1*incident_phi;
A3_relative_phase = -2*incident_phi;
A4_relative_phase = -3*incident_phi;

AF_VA = exp(1i*A1_relative_phase) + exp(1i*A2_relative_phase) + exp(1i*A3_relative_phase) + exp(1i*A4_relative_phase);
AF_VA_mag = abs(AF_VA)/4;

for m = 1:500;
    for n = 1:500;
        ElementPattern_HW(m,n) = cos(pi/2*sin(theta(m))*cos(phi(n))) / sqrt((1 - (sin(theta(m))^2*cos(phi(n))^2)));
    end 
end;

% Note we have to plot the element factor in both planes
% Note we are postly interested in YZ (H) plane

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
title('Figure 4: Half-Wave Dipole Array Pattern E');

ArrayPatternT = ArrayPattern_HW_YZ';
figure;
polar(theta, ArrayPatternT(1,:));
title('Figure 5: Half-Wave Dipole Array Pattern H');
%% Now let's start changing the array factor based on varying delta_phi
    % I've proven using hand calculations that we can compensate each element
    % in subarrays of Van Atta Arrays to create a larger retrodirective array

    % We want to now know how Beam Parameters change with slight erros in
    % compensation values

    % I've chose a range of -0.5 to 0.5 delta as the errors to be evaluated.
    % This is an assumtion and can be changed later

close all
clc

for i = 1:500
    delta_phi_temp = linspace(-0.5*incident_phi(i),0.5*incident_phi(i), 500);
    delta_phi(i,:) = delta_phi_temp;
end

% for each value of theta and delta_phi calculate the new relative phase 
for i = 1:500
    for j = 1:500
        A1_rp_new(i,j) = delta_phi(i,j);
        A2_rp_new(i,j) = -1*(2*pi/lamda*d*sin(theta(i))) - delta_phi(i,j);
        A3_rp_new(i,j) = -2*(2*pi/lamda*d*sin(theta(i))) - delta_phi(i,j);
        A4_rp_new(i,j) = -3*(2*pi/lamda*d*sin(theta(i))) - delta_phi(i,j);
    end
end
% Calculate the new array factors as a function of theta and delta_L
for i = 1:500
    for j = 1:500
        AF_pseudo_VA(i,j) = exp(1i*A1_rp_new(i,j)) + exp(1i*A2_rp_new(i,j)) + exp(1i*A3_rp_new(i,j)) + exp(1i*A4_rp_new(i,j));
    end
end

% So at some point in AF_pseudo_VA.. I think 250,: the pattern should be
% the same as the original Van Atta Array. Let's prove this out before
% seeing how beam parameters change with delta_phi

% So let's plot H-Plane changes with values of 0.5 0.4 0.3 0.2 0.1
% delta_phi
% These should correspond to index 100 200 300 400 500
AF_normal = AF_pseudo_VA(:,250); 
AF_1 = AF_pseudo_VA(:,100);
AF_2 = AF_pseudo_VA(:,200);
AF_3 = AF_pseudo_VA(:,300);
AF_4 = AF_pseudo_VA(:,400);
delta_AF = [AF_normal AF_1 AF_3 AF_4];

for d = delta_AF
for m = 1:500;
    for n = 1:500;
        ArrayPattern_HW_YZ(m,n) = ElementPattern_HW(m,n) * d(m);
    end 
end
    ArrayPatternT = ArrayPattern_HW_YZ';
    figure(1)
    polar(theta, ArrayPatternT(1,:));
    legend('Corrected Van Atta','0.1 phi error in correction','0.2 phi error in correction','0.3 phi error in correction', '0.4 phi error in correction')
    title('Figure 5: 4 Element Half-Wave Dipole Array Pattern H Plane');
    hold on
end
hold off
%% Well here's a thought... what if I leave the element patttern out and
% just plot array factor....
close all 
clc

for d = delta_AF
    AF_VA_new = abs(d);
    figure(1)
    polar(theta, AF_VA_new');
    legend('Corrected Van Atta','0.1 phi error in correction','0.2 phi error in correction','0.3 phi error in correction', '0.4 phi error in correction')
    title('Figure 5: 4 Element Half-Wave Dipole Array Pattern H Plane');
    hold on
end

hold off