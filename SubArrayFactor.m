function [subAF] = SubArrayFactor(i,N,theta,i_phi_sub)
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

i_phi_subi = i_phi_sub(i);
i_phi_subi = cell2mat(i_phi_subi);
for n = 0:N-1
    for i = 1:size(theta,2)
        subAF(n+1,i) = exp(-1i*i_phi_subi(n+1,i));
    end
end
subAF = sum(subAF(1:4,:));
end