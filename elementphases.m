function [i_phi_sub] = elementphases(M,N,theta,i_phi)
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
A = repmat(N,M);
A = A(:,2)';
i_phi_sub = mat2cell(i_phi,A,size(theta,2));
end
