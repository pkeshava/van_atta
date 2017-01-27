%% Background
% I have come up with an algorithm to compensate subarrays of a large Van
% Atta array to create retrodirectivity. I have come up with a method to
% analyze beam parameters of the large array based on error in
% compenstation value. This code is an attempt to demonstrate that

clear all 
close all
clc
% Define angles for plotting of patterns
theta = linspace(0,2*pi,500);
phi = linspace(0,2*pi,500);
% Define parameters of incident wave
fo = 79e9;
c = 2.99792458e8;                                           % frequency of 
lamda = c/fo;                                               % free space wavelength
d = 0.8*lamda;                                              % element seperation
theta_i = pi/4;                                             % incident angle on array
M = 4;                                                      % Number of subarray
N = 4;                                                      % Number of elements per sub array
%delta_phi_1 = N*(M-1)*i_phi;                               % Compensation values for farthest array
%delta_phi = linspace(N*i_phi,delta_phi_1,M/2);             % Vector of compensation values for each subarray  
%i_phi_vector = 0:-i_phi:-(M*N-1)*i_phi;                    % relative incident phases for all elements

%% Prior to compensating analysis
% Plot AF for each array
% Plot total AF assuming van atta
close all 
clc
[i_phi] = IncidentPhases(M,N,theta,lamda,d);
[i_phi_sub] = elementphases(M,N,theta,i_phi);
subAF1 = SubArrayFactor(1,N,theta,i_phi_sub);
figure(1)
polar(theta, subAF1);

%% Plot compensated total AF
% Note the direction should be the same as array 1 but with a thinner
% beamwidth

%% Plot compensated total AF vs Different error values
delta_phi_e = 0.3*i_phi;                                    % relative error in compensation.. picked semi arbitraly








