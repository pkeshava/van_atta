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
%% Prior to compensating analysis: Plot total AF Assuming Correct Van Atta Compensation ie using linear excitation
[i_phi] = IncidentPhases(M,N,theta,lamda,d);
[AF_T] = TotalArrayFactor(N,M,theta,i_phi);
figure;
polar(theta, AF_T);
title(['Array Factor of ' num2str(M) ' Van Atta Sub-Arrays Containing ' num2str(N) ' Elements and Correctly Compensation']);
%% COMPENSATION ANALYSIS
% First rearrange phases correctly to simulate Van Atta behaviour then use
% compensation algorithm 
clc
% Need to design scrambling algorithm for i_phi: [phi_tx] = RxToTxPhase(N,M,theta,lamda,d);
error_percent = 0;
[compensated_phases] = Compensation(M,N,theta,lamda,d,error_percent);
% Plot
[AF_c] = TotalArrayFactor(N,M,theta,compensated_phases);
figure;
polar(theta, AF_c);
title(['Array Factor of ' num2str(M) ' Van Atta Sub-Arrays Containing ' num2str(N) ' Elements and Correctly Compensated']);
%% COMPENSATION ERROR ANALYSIS
% Let's start approaching this by doing it as a percentage in compensation
% error

error_percent = 0.50;
[compensated_phases] = Compensation(M,N,theta,lamda,d,error_percent);
% Plot
[AF_c] = TotalArrayFactor(N,M,theta,compensated_phases);
figure;
polar(theta, AF_c);
title(['Array Factor of ' num2str(M) ' Van Atta Sub-Arrays Containing ' num2str(N) ' Elements and With ' num2str(error_percent*100) '% Error In Compensation Vector']);

%% COMPENSATING AS FRACTION OF PHI - ERROR ANALYSIS

% Basically capable of shifting with some sort of resolution
% Define resolution variable in degrees, convert to radians
phase_shift_resolution = 10;
phase_shift_resolution = phase_shift_resolution*pi/180;
error_percent = 0.50;
[compensated_phases] = Compensation(M,N,theta,lamda,d,error_percent, phase_shift_resolution);

% Round the required correction to the closest bin defined by resolution








