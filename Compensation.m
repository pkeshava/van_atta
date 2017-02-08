function [compensation_i_phi_sub] = Compensation(M,N,theta,lamda,d,error,res)
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
%   d: element spacing
%   lamda: wavelength
%   error: error for compensation in percent
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

delta_phi_1 = -N*(M-1);                              
phi_compensation_multiplier = linspace(delta_phi_1,-N,M/2);
phi_compensation_multiplier = [phi_compensation_multiplier -fliplr(phi_compensation_multiplier)];
[i_phi] = IncidentPhases(M,N,theta,lamda,d);
i_phi = i_phi(2,:);

for i = 1:M
    comp_vector(i,:) = phi_compensation_multiplier(i)*i_phi;
end

for i = 1:M
    for j = 1:size(theta,2)
        comp_vect_r(i,j) = round2(comp_vector(i,j),res);
    end
end

%comp_vector = comp_vector*(1+error);

[phi_tx] = RxToTxPhase(N,M,theta,lamda,d);
[tx_phi_sub] = elementphases(M,N,theta,phi_tx);
compensation_i_phi_sub = [];
for i = 1:M
    sub_array_i = cell2mat(tx_phi_sub(i)); %select subarray
    for j = 1:N
       sub_array_i(j,:) = sub_array_i(j,:) + comp_vector(i,:);
    end
    compensation_i_phi_sub = [compensation_i_phi_sub; sub_array_i];
end

%compensated_phases = compensation_i_phi_sub + phi_tx;
end