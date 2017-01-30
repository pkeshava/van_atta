function [compensated_phases] = Compensation(M,N,theta,lamda,d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes a matrix that has all the relative phases for each
% element in a pseudo van atta array with M subarrays each with N elements
% (are you confused yet?) calculated for every arrival angle theta and
% splits it into a cell that has N cells of M x theta(resolution)
%
%   M: number of subarrays
%   N: number of elements in subarray
%   theta: vector of angles (size defines resolution of plot)
%   i_phi: a Mx1 matrix, each row containing size(theta,2) values. Contains
%       incident phases for each element in all subarrays at all angles of
%       theta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_phi_1 = -N*(M-1);                              
phi_compensation_multiplier = linspace(delta_phi_1,-N,M/2);
phi_compensation_multiplier = [phi_compensation_multiplier -fliplr(phi_compensation_multiplier)];

[phi_tx] = RxToTxPhase(N,M,theta,lamda,d);
[i_phi] = IncidentPhases(M,N,theta,lamda,d);    % calculate incident phases for all elements
[i_phi_sub] = elementphases(M,N,theta,i_phi);

compensation_i_phi_sub = [];
for i = 1:M
   sub_array_i = cell2mat(i_phi_sub(i));
   compensation_i_phi_temp = phi_compensation_multiplier(i)*sub_array_i;
   compensation_i_phi_sub = [compensation_i_phi_sub; compensation_i_phi_temp];
end

compensated_phases = compensation_i_phi_sub + phi_tx;
end