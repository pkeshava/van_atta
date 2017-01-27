function [i_phi_sub] = elementphases(M,N,theta,i_phi)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes a matrix that has all the relative phases for each
% element in a pseudo van atta array with M subarrays each with N elements
% (are you confused yet?) calculated for every arrival angle theta and
% splits it into a cell that has N cells of M x theta(resolution)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = reshape(i_phi,4,[]);
B = repmat(size(theta),N);
B = B(:,2)';
i_phi_sub = mat2cell(X,M,B);
end
