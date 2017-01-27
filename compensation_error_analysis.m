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
%i_phi = 2*pi/lamda*d*cos(theta);                            % relative incident phases
M = 4;                                                      % Number of subarray
N = 4;                                                      % Number of elements per sub array
% for every elemt

%delta_phi_1 = N*(M-1)*i_phi;                                % Compensation values for farthest array
%delta_phi = linspace(N*i_phi,delta_phi_1,M/2);              % Vector of compensation values for each subarray  
%i_phi_vector = 0:-i_phi:-(M*N-1)*i_phi;                     % relative incident phases for all elements
%i_phi_vector = wrapTo2Pi(i_phi_vector);                         % normalize to 2*pi rad 

%% Prior to compensating analysis
% Plot AF for each array
% Plot total AF
close all
clc
for n = 0:N*M-1
    for i = 1:size(theta,2)
        AF_values(n+1,i) = exp(-1i*n*2*pi/lamda*d*cos(theta(i)));
    end
end

%AF_values_mag = abs(AF_values);
AF_total = sum(AF_values(1:16,:));
AF_element1 = sum(AF_values(1:4,:));
AF_element2 = sum(AF_values(5:8,:));
AF_element3 = sum(AF_values(9:12,:));
AF_element4 = sum(AF_values(13:16,:));
figure(1)
polar(theta, AF_total);

%%
AF_total = [AF_element1; AF_element2; AF_element3; AF_element4];
for i = 1:M
    figure(1)
    polar(theta, AF_total(i,:));
    legend('Subarray 1','Subarray 2','Subarray 3','Subarray 4')
    title('Figure 5: 4 Element Half-Wave Dipole Array Pattern H Plane');
    hold on
end

hold off

%% Plot compensated total AF
% Note the direction should be the same as array 1 but with a thinner
% beamwidth

%% Plot compensated total AF vs Different error values
delta_phi_e = 0.3*i_phi;                                    % relative error in compensation.. picked semi arbitraly








