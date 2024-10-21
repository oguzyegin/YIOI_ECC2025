%%
% This code is published under the BSD3-Clause License.
% Copyright (c) 2015, Joerg.Fehr
%%
%   All rights reserved.
%  Redistribution and use in source and binary forms, with or without modification,
%  are permitted provided that the following conditions are met:
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
% All rights reserved. 
% Joerg.Fehr, University of Stuttgart (c) 2015
% Jens Saak, MPI Magdeburg (c) 2016


% Algorithm tested with 
%    Matlab 8.6.0.267246 (R2015b) Linux/distru 17-Mar-2016
%    Matlab 8.2.0.701 (R2013b) on Win. 7 /distru 17-Mar-2016
%    Matlab 8.4.0.150421 (R2014b) Linux/Ubuntu 12.04LTS 14-Apr-2016
%    Octave 4.0 with control toolbox on Win. 7/distru 17-Mar-2016
%    Octave 4.0.1 with control toolbox 3.0.0 Linux/CentOS 7.1 14-Apr-2016

% The example used to test the IRKA algorithm is the
% FOM example from the SLICOT Library
% http://slicot.org/20-site/126-benchmark-examples-for-model-reduction
% downloaded at 21.12.2015 by Joerg Fehr, ITM Uni Stuttgart
clear all;
FOM_Ex = load('fom.mat');
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%')
disp('% The FOM example from the Slicot collection')
disp('%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%%
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% The single shift case, i.e. reduced order 1')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
close all

b = randn(size(FOM_Ex.B,2),1);
c = randn(1,size(FOM_Ex.C,1));

[V1,W1,A_red1,B_red1,C_red1,E_red1,s1]=IRKA(FOM_Ex.A,... % system Mat. FOM ex.
    FOM_Ex.B,...                                  % input  Mat. FOM ex.
    FOM_Ex.C,...                                  % output Mat. FOM ex.
    speye(size(FOM_Ex.A,1)),...                   % descriptor Matrix
    1,...                                         % initial shift value
    b,...                                         % tangential direction for B
    c,...                                         % tangential direction for C
    25,...                                        % maximum iteration number
    1.0e-3,...                                    % convergence tolerance
    1);                                           % debug flag
disp('Press any key to continue')
pause;
%% Test IRKA with multiple shift
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% Double shift -- reduced order 2')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
close all

nred = 2 ;
b = randn(size(FOM_Ex.B,2),nred);
c = randn(nred,size(FOM_Ex.C,1));

[V2,W2,A_red2,B_red2,C_red2,E_red2,s2]=IRKA(FOM_Ex.A,... % system Mat. FOM ex.
    FOM_Ex.B,...                                  % input  Mat. FOM ex.
    FOM_Ex.C,...                                  % output Mat. FOM ex.
    speye(size(FOM_Ex.A,1)),...                   % descriptor Matrix
    [1 3 ]',...                                   % initial shift values
    b,...                                         % tangential direction for B
    c,...                                         % tangential direction for C
    25,...                                        % maximum iteration number
    1.0e-3,...                                    % convergence tolerance
    1);                                           % debug flag
disp('Press any key to continue')
pause;
%% Test IRKA with odd shift number
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('% Multishift example -- reduced order 5')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
close all
nred = 5;
b = randn(size(FOM_Ex.B,2),nred);
c = randn(nred,size(FOM_Ex.C,1));

[V5,W5,A_red5,B_red5,C_red5,E_red5,s5]=IRKA(FOM_Ex.A,... % system Mat. FOM ex.
    FOM_Ex.B,...                                  % input  Mat. FOM ex.
    FOM_Ex.C,...                                  % output Mat. FOM ex.
    speye(size(FOM_Ex.A,1)),...                   % descriptor Matrix
    [1 3 4 5 243 ]',...                           % initial shift values
    b,...                                         % tangential direction for B
    c,...                                         % tangential direction for C
    25,...                                        % maximum iteration number
    1.0e-3,...                                    % convergence tolerance
    1);                                           % debug output on
disp('Press any key to continue')
pause;

%% Test IRKA with a second Example
% The example used to test the IRKA algorithm is the
% Rail Example from the Oberwolfach Model Reduction Benchmark Collection
% https://portal.uni-freiburg.de/imteksimulation/downloads/benchmark/Steel%20Profiles%20(38881)
% 17-Mar-2016 
% To import the rail example MatrixMarket I/O Functions for Matlab 
% http://math.nist.gov/MatrixMarket/ downloaded at 17-Mar-2016 was used
% and later stored in the file Rail_20209.mat
% [Rail_20209.A] = mmread('rail_20209_c60.A');
% [Rail_20209.B] = mmread('rail_20209_c60.B');
% [Rail_20209.C] = mmread('rail_20209_c60.C');
% [Rail_20209.E] = mmread('rail_20209_c60.E');
% save('Rail_20209.mat','Rail_20209')

load('Rail_20209.mat')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%')
disp('% The Rail example from the Oberwolfach collection')
disp('%')
disp('% Here only input 6 and output 2 are used to have a SISO system')
disp('% (as in [Gugercin/Antoulas/Beattie 08])')
disp('%')
disp('% Multishift example -- reduced order 6')
disp('%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Similar to IRKA example from GugercinAntoulasBeattie08 the SISO example
% with the 6th Input from B and the 2nd output from C is reduced with the
% IRKA algorithm
close all
% Transfer function
% Settings used to calculate the transfer function of original system
fmin = 10^(-8)/2/pi; fmax = 100/2/pi;  afreq_interval = [2*pi*fmin, 2*pi*fmax];
npts=1000;afreq = logspace(log10(afreq_interval(1)),log10(afreq_interval(2)),npts); 
s = 1i*afreq;
%
% Computation of the transfer function of the original system
[H_full_siso,FroNorm_H_siso,s,~]=calculateFrequencyResponse(s,Rail_20209.A,Rail_20209.B(:,6),...
    Rail_20209.C(2,:),Rail_20209.E);
figure(6),loglog(s/1i,FroNorm_H_siso,'r','DisplayName','original'),hold on
%%
%
nred = 6;
V = orth(randn(20209,nred));
sigma0 = -eig(V'*(Rail_20209.A*V),V'*(Rail_20209.E*V));
b = randn(1,nred);
c = randn(nred,1);

[V6s,W6s,A_red,B_red,C_red,E_red,~]=IRKA(Rail_20209.A,... % system Mat. Rail ex.
    Rail_20209.B(:,6),...                             % input Mat. Rail ex.
    Rail_20209.C(2,:),...                             % output Mat. Rail ex.
    Rail_20209.E,...                                  % descriptor Matrix
    sigma0,...                                        % initial shift values
    b,...                                             % tangential direction for B
    c,...                                             % tangential direction for C
    25,...                                            % maximum iteration number
    1.0e-3,...                                        % convergence tolerance
    0);                                               % debug output off
% comparison of the transfer fucntion of the reduced system with the 
% transfer function of the original system
s = 1i*afreq;
[H_fullred_siso,FroNorm_Hred_siso,s,~]=calculateFrequencyResponse(s,A_red,B_red,C_red,E_red);
figure(6),  loglog(s/1i,FroNorm_Hred_siso,'b--'),  hold on
xlabel('frequency \omega [rad/sec]'),  ylabel('||H(i\omega)||_F')
legend show;

disp('Press any key to continue')
pause;

%%
%
close all
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%')
disp('% The Rail example from the Oberwolfach collection')
disp('%')
disp('% The MIMO case using tangential interpolation')
disp('%')
disp('% Multishift example -- reduced order 6')
disp('%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%
% Computation of the transfer function of the original system
[H_full_mimo,FroNorm_H_mimo,s,~]=calculateFrequencyResponse(s,Rail_20209.A,Rail_20209.B,...
    Rail_20209.C,Rail_20209.E);
figure(6),loglog(s/1i,FroNorm_H_mimo,'r','DisplayName','original'),hold on
%%
nred = 50;
V = orth(randn(20209,nred));
sigma0 = -eig(V'*(Rail_20209.A*V),V'*(Rail_20209.E*V));
b = randn(size(Rail_20209.B,2),nred);
c = randn(nred,size(Rail_20209.C,1));

[V,W,A_red,B_red,C_red,E_red,~]=IRKA(Rail_20209.A,... % system Mat. Rail ex.
    Rail_20209.B,...                                  % input Mat. Rail ex.
    Rail_20209.C,...                                  % output Mat. Rail ex.
    Rail_20209.E,...                                  % descriptor Matrix
    sigma0,...                                        % initial shift values
    b,...                                             % tangential direction for B
    c,...                                             % tangential direction for C
    25,...                                            % maximum iteration number
    1.0e-3,...                                        % convergence tolerance
    0);                                               % debug output off
% comparison of the transfer fucntion of the reduced system with the 
% transfer function of the original system
s = 1i*afreq;
[H_fullred_mimo,FroNorm_Hred_mimo,s,~]=calculateFrequencyResponse(s,A_red,B_red,C_red,E_red);
figure(6),  loglog(s/1i,FroNorm_Hred_mimo,'b--','DisplayName','reduced'),  hold on
xlabel('frequency \omega [rad/sec]'),  ylabel('||H(i\omega)||_F')
legend show;

