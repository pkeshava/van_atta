%% ENEL 627 Assignment 3
% Analysis of Van Atta Array using different lengths
% code built off of ENEL 627 assignment 
%
%% Van Atta array parameters
% Frequency 79 GHz
% Transmission line length of first pair = l1
% Transmission line length of second pair = l2
% Element spacing = d
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
% Define angles for plotting of patterns
theta = linspace(0,2*pi,500);
phi = linspace(0,2*pi,500);
% Define parameters of incident wave
fo = 79e9;
c = 2.99792458e8;
lamda = 1/fo;
d = 0.5*lamda;
theta_i = linspace(0,pi,500);
incident_phi = 2*pi/lamda*d*cos(theta_i);

% Define parameters of 4 element array based on incident phase
phase = [zeros(1,500); -1*incident_phi; -2*incident_phi; -3*incident_phi];

%%%%%%%%%%

% These lines of code below are a little bit useless but they demostrate an
% interesting concept which is it shows that the array factor is actually
% indpendent of the actual phase value but rather dependent on the
% difference in relative phases. I've essentially shown that a van atta
% array is a van atta array

    % for m = 1:500       % For every value of theta
    %     for n = 1:500   % For every value of incident phi
    %         AF_VA(m,n) = exp(1i*phase(1,n)) + exp(1i*phase(2,n)) + exp(1i*phase(2,n)) + exp(1i*phase(3,n));
    %     end
    % end
    
%%%%%%%%%% 

% REMEMBER I'm supposed to be thinking of an array of arrays as well!!!
% Let's handle this now



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
%% Now I have to correct some stuff... I don't think theta incident 
% and theta for the plot should be the same variable... 