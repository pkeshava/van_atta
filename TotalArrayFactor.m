function [AF_T] = TotalArrayFactor(N,M,theta,i_phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes the cell created by elementphases.m and calculates
% the array factor specified by the index i
%
%   i: index of the subarray chosen for array factor calculation values
%      from 1-M
%   N: number of elements in subarray
%   theta: vector of angles (size defines resolution of plot)
%   i_phi_sub: cell of incident phase angles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for n = 0:M*N-1 
    for i = 1:size(theta,2)
        AF_T(n+1,i) = exp(-1i*i_phi(n+1,i));
    end
end
AF_T = sum(AF_T(1:M*N,:))/(N*M);
end