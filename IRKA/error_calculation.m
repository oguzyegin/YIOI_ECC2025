function [H_red,H_red_fro, H_E_2, H_E_fro]=error_calculation(A_red,B_red,C_red,E_red,s,H_full)
%calculateFrequencyResponse Calculate the frequency response of a LTI first
%   order system
%
% Syntax
% ======
%   [H_red,H_red_fro, H_E_2, H_E_fro]=error_calculation(A_red,B_red,C_red,s,H_full)
%
%
% Description
% ===========
%   Calculates the frequency response of a reduced LTI system and the error
%   between the frequency response of full system and the reduced system
%
%   The frequency response matrix and the error is calculated at the frequencies s=i*omega
%
%   Reminder: It is considered that the frequency response matrix of the
%   full system is calculated at the same frequencies
%
% Input Arguments
% ===============
%
%   Required input arguments - pass in this order
%   ---------------------------------------------
%
%        -'A': system matrix
%        -'B': input matrix
%        -'C': output matrix
%        -'s': Vector of the imaginary points where the system response
%           matrix is calculated
%        -'E': descriptor matrix
%
%
% Output Arguments
% ================
%       -H_red: System Response Matrix at  s
%       -H_red_fro: Frobenius Norm of the System Response Matrix at
%              s - FroNorm_H(afreq)=norm(H{afreq},'fro');
%       - H_E_2: 2-Norm of the  error
%       - H_E_fro: Frobenius norm of the relative error
%
% See also
%   norm
% History
% =======
%  Written by J. Fehr 11th January 2016
%           after calculateFrequencyResponse
%           backslash operator instead of lu(...)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Calculate transfer function of the reduced system.
npts = (length(s));
% E = eye(size(A_red));
H_red = zeros(npts,size(C_red,1),size(B_red,2));
H_red_fro = zeros(1,npts);
H_E_2 = zeros(1,npts);
H_E_fro = zeros(1,npts);

for k = 1:npts
    sEmA = s(k)*E_red-A_red;
    H_red(k,:,:) = C_red * (sEmA\B_red);
    H_red_fro(k) = norm( H_red(k,:,:),'fro');
end

for k = 1:npts
    % take index (1,1) of H_full matrix sind B and C are now only the
    % first dof of the first node
    H_E_2(k)   = norm(H_full{k}-H_red(k,:,:),2);
    H_E_fro(k) = norm(H_full{k}-H_red(k,:,:),'fro');
end

