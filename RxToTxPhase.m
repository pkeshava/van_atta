function [phi_tx] = RxToTxPhase(N,M,theta,lamda,d)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function takes phases for a subarray van atta with regards to
% incident phase and passes the corresponding relative phases to the now
% transmitting element (mimicing Van Atta behaviour). This is useful for
% the next step which is to simulate beam parameters withr respect to 
% compensation error
%
%   M: Number of subarrays
%   N: number of elements in subarray
%   theta: vector of angles (size defines resolution of plot)
%   lamda: free space wavelength
%   d: element spacing in wavelengths
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[i_phi] = IncidentPhases(M,N,theta,lamda,d);    % calculate incident phases for all elements
[i_phi_sub] = elementphases(M,N,theta,i_phi);   % segment results into appropriate cells of subarrays
phi_tx = [];
for i = 1:M
    i_phi = cell2mat(i_phi_sub(i));
    phi_tx_i = flipud(i_phi);                   % flips relative phases of elements within array to transmit element i.e. simulating Van Atta behaviour
    phi_tx = [phi_tx; phi_tx_i];    
end

end