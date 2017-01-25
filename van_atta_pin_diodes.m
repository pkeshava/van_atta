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
incident_phi = 2*pi/lamda*d*sin(theta);

% define parameters of array
%l1 = 100;                               % dimensions in mm
A1_relative_phase = zeros(1,500);
A2_relative_phase = -1*incident_phi;
A3_relative_phase = -2*incident_phi;
A4_relative_phase = -3*incident_phi;
delta_phi = zeros(500,500);
for i = 1:500
    delta_phi_temp = linspace(-0.5*incident_phi(i),0.5*incident_phi(i), 500);
    delta_phi = [delta_phi; delta_phi_temp];
end

AF_VA = exp(1i*A1_relative_phase) + exp(1i*A2_relative_phase) + exp(1i*A3_relative_phase) + exp(1i*A4_relative_phase);
AF_VA_mag = abs(AF_VA)/4;

for m = 1:500;
    for n = 1:500;
        ElementPattern_HW(m,n) = cos(pi/2*sin(theta(m))*cos(phi(n))) / sqrt((1 - (sin(theta(m))^2*cos(phi(n))^2)));
    end 
end;

% Note we have to plot the element factor in both planes
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
%% Now let's start changing the array factor based on varying length
close all
clc
% for each value of theta and delta_l calculate the new relative phase 
A1_rp_new = zeros(500,500);
for i = 1:500
    for j = 1:500
        A2_rp_new(i,j) = -1*(2*pi/lamda*d*sin(theta(i))) - l2_electrical(j);
        A3_rp_new(i,j) = -2*(2*pi/lamda*d*sin(theta(i))) - l2_electrical(j);
        A4_rp_new(i,j) = -3*(2*pi/lamda*d*sin(theta(i))) - l2_electrical(j);
    end
end
% Calculate the new array factors as a function of theta and delta_L
for i = 1:500
    for j = 1:500
        AF_pseudo_VA(i,j) = exp(1i*A1_rp_new(i,j)) + exp(1i*A2_rp_new(i,j)) + exp(1i*A3_rp_new(i,j)) + exp(1i*A4_rp_new(i,j));
    end
end


