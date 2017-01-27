function [i_phi] = IncidentPhases(M,N,theta,lamda,d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes a matrix that has all the relative phases for each
% element in a pseudo van atta array with M subarrays each with N elements
% (are you confused yet?) calculated for every arrival angle theta and
% splits it into a cell that has N cells of M x theta(resolution)
%
%   M: number of subarrays
%   N: number of elements in subarray
%   theta: vector of angles (size defines resolution of plot)
%   i_phi: a 4x1 cell, each cell containing size(theta,2) values. Contains
%       incident phases for each element in all subarrays at all angles of
%       theta
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 0:N*M-1
    for i = 1:size(theta,2)
        i_phi(n+1,i) = n*2*pi/lamda*d*cos(theta(i));
    end
end
end