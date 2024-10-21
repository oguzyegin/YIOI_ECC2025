function [V,W,A_red,B_red,C_red,E_red,sigma]=IRKA(A,B,C,E,sigma,b,c,maxiter,contol,varargin)
%  function   [V,W,A_red,B_red,C_red,sigma]=IRKA(A,B,C,sigma,n,contol,varargin)
%  IRKA (iterative rational Krylov algorithm) is an interpolation-based model
%  reduction method for SISO linear time invariant systems based on a
%  Rational Krylov Algorithm.
%  The algorithm calculates in an iterative manner a reduced order model
%  that satisfies a first-order necessary conditions for H2 optimality
%  model reduction. The algorithms finds a local minimizer of the H2 MOR
%  problem
%
%  Input:  A,B,C,E:  Matrices of an SISO LTI-System
%          sigma :   initial guess of complex shifts given as a row vector
%          b,c :     initial tangential directions
%          maxiter:  maximum number of iterations taken
%          contol:   convergence tolerance of IRKA
%          Debug :   (optional=1) debug flag, if 1 (default) the transfer
%                    functions of original and iteratively reduced systems
%                    are calculated and visualized. In addition the H2 and
%                    H-infiniy norm of the error systems ||H-H_red||_H2/Hinf
%                    are calculated, be careful error calculation is very
%                    expensive (calls lyap (dense!) for the error system)
%
%  Output: A_red,B_red,C_red,E_red: Matrices of a reduced LTI-System dim nred
%          V,W:           right and left reduction matrices
%          sigma_out:	  final shifts
%  Example:
%  [V,W,A_red,B_red,C_red,s]=IRKA(A,B,C,E,sigma,nred,tol,1)
%
%  original IRKA idea from
%  S. Gugercin; A.C. Antoulas; C. Beattie: "H2 Model Reduction for Large-
%  Scale Linear Dynamical Systems", SIAM. J. Matrix Anal. & Appl., vol.30, no.2,
%  pp.609-638, 2008.
%  Algorithm 4.1 (which only works for SISO systems)
%  Therefore, a tangential Krylov is used, the implementation of this version
%  uses ideas from 
%  H.~K.~F. Panzer: "Model Order Reduction by Krylov Subspace Methods
%  with Global Error Bounds and Automatic Choice of Parameters", Dissertation,
%  Technische Universität München, München, 2014.
%
% This code is published under the BSD3-Clause License.
% Copyright (c) 2015, Joerg.Fehr ITM University of Stuttgart,
% http://www.itm.uni-stuttgart.de
% Copyright (c) 2016, Jens Saak CSC MPI Magdeburg
% http://www.mpi-madgeburg.mpg.de/csc/
%


%%
%  Redistribution and use in source and binary forms, with or without modi-
%  fication, are permitted provided that the following conditions are met:
%  1. Redistributions of source code must retain the above copyright
%   notice, this list of conditions and the following disclaimer.
%  2. Redistributions in binary form must reproduce the above copyright
%   notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution.
%  3. Neither the name of the copyright holder nor the names of its
%   contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
%  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
%  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
%  A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
%  HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
%  iNCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
%  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
%  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
%  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
%  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
%  DAMAGE.
% All rights reserved. (c) 2015
% Joerg.Fehr, University of Stuttgart (c) 2015
% Jens Saak, MPI Magdeburg (c) 2016


% initial guess for sigma = 1 (since the algorithm is instable, a stable
% sigma can vary from system to system)
% sigma = 1;
switch nargin
    case 9
        Debug = 1;
        disp('verbose error computation: on!');
    case 10
        Debug = varargin{1};
        if Debug
            disp('verbose error computation: on!');
        else
            disp('verbose error computation: off!');
        end
    otherwise
        Debug = varargin{1};
        if Debug
            disp('verbose error computation: on!');
        else
            disp('verbose error computation: off!');
        end
        warning('To many inputs 9th and further inputs are not considered')
end

% Octave seems to be slightly less acurate in eigenvalue computation. We
% need to increase the tolerance for eigenvalue pairing a bit there:
if exist('OCTAVE_VERSION', 'builtin')
    cplxpair_tol = 1000.0*eps; % increased tolerance for Octave
else
    cplxpair_tol = 100.0*eps; % the default value used in Matlab.
end
%%  initalize standard parameters
n_org=size(E,1);    sigma_change = inf;    iter=1; 
fprintf('IRKA iteration: %3d\n',iter);

%% To check the quality of the reduced model a comparison with the original
% system is performed if desired, therefore the p-coded non-disclosure function
% calculateFreqeuncyResponse.m derived from the non-public Morembs package
% is used.
if Debug==1
    %%  initalize standard parameters
    npts=10000;  H2_norm = inf(1,maxiter); Hinf_norm = inf(1,maxiter);  
    sigmas= zeros(length(sigma),maxiter);
    H2_norm_approx = inf(1,maxiter);
    
    %% Settings used to calculate the transfer function
    fmin = 0.1; fmax = 1000; afreq_interval = [2.0*pi*fmin, 2.0*pi*fmax];
    afreq = logspace(log10(afreq_interval(1)),log10(afreq_interval(2)),npts);
    s = 1i*afreq;
    
    [H_full,FroNorm_H,s,~] = calculateFrequencyResponse(s,A,B,C,E);
    figure(1),loglog(s/2.0/pi/1i,FroNorm_H),hold on
end
%%
% Use Tangential Krylov algorithm by Heiko Panzer 
% from [Panzer14] to calculate V and W. Only real systems are
% calculated
% Sorting 
sigma = cplxpair(sigma, cplxpair_tol);
[V, ~, ~, W, ~, ~] = TangentialKrylov(A,B,C,E,sigma,b,c);
A_red = full(W'*(A*V));   E_red = full(W'*E*V);    
B_red = full(W'*B);    C_red = full(C*V);
[eigenvecs, eigenvals] = eig(A_red,E_red);  % calculate the eigenvalues
sigma  =  cplxpair(-diag(eigenvals), cplxpair_tol);
b = (eigenvecs\(E_red\B_red)).';
c = (C_red*eigenvecs).';

while (sigma_change > contol && iter <= maxiter)
    [V, ~, ~, W, ~, ~] = TangentialKrylov(A,B,C,E,sigma,b,c);
    sigma_old = sigma;
    A_red = full(W'*(A*V)); E_red=full(W'*E*V);
    B_red = full(W'*B);    C_red = full(C*V);
    % Compute new shifts and tangential directions and continue algorithm.
    [eigenvecs, eigenvals] = eig(A_red,E_red);  % calculate the eigenvalues
    sigma  =  cplxpair(-diag(eigenvals), cplxpair_tol); %pair conjugates
    % update directions:
    b = (eigenvecs\(E_red\B_red)).'; 
    c = (C_red*eigenvecs).';

                                          
    % Compute the relative change for the stopping criterion
    sigma_change = norm(sigma - sigma_old)/norm(sigma_old); % difference
    iter = iter+1;
    fprintf('IRKA iteration: %3d\n',iter);

    if Debug
        % Calculate the Frobenius norm of the reduced system,
        % the 2-norm and the Frobenius norm of the output error
        % e=|H(i2\pi f)-H_red(i2\pi f)|_{2/Fro}
        n_red = size(A_red,1);
        [~, H_IRKA_fro, H_E_2, H_E_fro] = error_calculation(A_red,B_red,C_red,E_red,s,H_full);
        Error_System.E = [E,zeros(n_org,n_red);zeros(n_red,n_org),E_red];
        Error_System.A = [A,zeros(n_org,n_red); zeros(n_red,n_org),A_red];
        Error_System.B = [B; B_red];
        Error_System.C = [C, -C_red];
        P = lyap( Error_System.A, Error_System.B*Error_System.B', ...
                  [], Error_System.E);
        H2_norm_lyap = sqrt(trace(Error_System.C*P*Error_System.C'));
        H2_norm(iter-1) = H2_norm_lyap;
        %  approximated H2-error
        H2_norm_approx(iter-1) = sqrt(1.0/2.0/pi*2.0*trapz(afreq,(H_E_fro).^2));
        Hinf_norm(iter-1) = max(H_E_2);
        % Plot the errors
        color = rand(3,1);
        figure(1)
        loglog(s/2.0/pi/1i,H_IRKA_fro,'Color',color,...
            'DisplayName',['Iter',num2str(size(H2_norm_approx,1))])
        hold on
        xlabel('frequency [1/s]');   ylabel('|| H ||_F [-]');
        drawnow;
        
        figure(2)
        loglog(s/2/pi/1i,H_E_fro,'Color',color,...
            'DisplayName',['Iter',num2str(size(H2_norm_approx,1))]),
        hold on
        xlabel('frequency [1/s]');   ylabel('error in H2 norm [-]');
        drawnow;
        
        % Save the sigmas
        sigmas(:,iter-1) = sigma;
        figure(3)
        plot3((iter-1)*ones(size(sigma,1),size(sigma,2)),...
              real(sigma),imag(sigma),'Color',color, 'LineStyle','none',...
              'Marker','*','DisplayName',['Iter',num2str(size(H2_norm,2))])
        hold on
        title('evolution of the IRKA shifts')
        xlabel('iteration [-]');    
        ylabel('real part [1/s]');   
        zlabel('imag part [1/s]');
        drawnow;
    end % Debug
    
end
B_red = full(W'*B); C_red = full(C*V); 

% Plot the errors of each iteration.
if Debug==1
    iter = iter-1;
    figure(4);
    plot(1:iter, H2_norm(1:iter),'r','DisplayName','lyap');
    hold on
    plot(1:iter, H2_norm_approx(1:iter),'b', 'DisplayName','approx');
    title('IRKA of LTI system');
    xlabel('iteration [-]');   ylabel('error in H2-norm [-]');
    legend('show');
    
    figure(5);
    plot(1:iter, Hinf_norm(1:iter),'b');
    title('IRKA of LTI system');   legend('H_{infty} error'); 
    xlabel('iteration [-]');   ylabel('error [-]');
end

