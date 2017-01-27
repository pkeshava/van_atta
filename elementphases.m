function [i_phi_sub] = elementphases(M,N,theta,i_phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes a matrix that has all the relative phases for each
% element in a pseudo van atta array with M subarrays each with N elements
% (are you confused yet?) calculated for every arrival angle theta and
% splits it into a cell that has N cells of M x theta(resolution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A = repmat(N,M);
A = A(:,2)';
i_phi_sub = mat2cell(i_phi,A,size(theta,2));
end
