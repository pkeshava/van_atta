%% Background
% I have come up with an algorithm to compensate subarrays of a large Van
% Atta array to create retrodirectivity. I have come up with a method to
% analyze beam parameters of the large array based on error in
% compenstation value. This code is an attempt to demonstrate that

clear all 
close all
clc
% Define angles for plotting of patterns
theta = linspace(0,2*pi,5000);
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

%% Prior to compensating analysis: Plot AF for each array
close all 
clc
[i_phi] = IncidentPhases(M,N,theta,lamda,d);    % calculate incident phases for all elements
[i_phi_sub] = elementphases(M,N,theta,i_phi);   % segment results into appropriate cells of subarrays
for i = 1:M
    figure;
    [subAFi] = SubArrayFactor(i,N,theta,i_phi_sub);   % calculate AF for subarray indicated by index                                       
    polar(theta, subAFi);
    title(['Array Factor for Sub-Array ' num2str(i) ' of' num2str(M)]);
end
%% Prior to compensating analysis: Plot total AF Assuming CorrectVan Atta Compensation
[i_phi] = IncidentPhases(M,N,theta,lamda,d);
[AF_T] = TotalArrayFactor(N,M,theta,i_phi);
figure;
polar(theta, AF_T);
title(['Total Array Factor of ' num2str(M) ' Sub Van Atta Arrays Containing ' num2str(N) ' Elements and Correctly Compensation']);
%% COMPENSATION ERROR ANALYSIS
% Need to design scrambling algorithm for i_phi
[i_phi] = IncidentPhases(M,N,theta,lamda,d);    % calculate incident phases for all elements
[i_phi_sub] = elementphases(M,N,theta,i_phi);   % segment results into appropriate cells of subarrays
i_phi1 = cell2mat(i_phi_sub(1));
i_phi1_temp = flipud(i_phi1);                   % flips relative phases of elements within array to transmit element i.e. simulating Van Atta behaviour 


%% Plot compensated total AF vs Different error values
delta_phi_e = 0.3*i_phi;                                    % relative error in compensation.. picked semi arbitraly








