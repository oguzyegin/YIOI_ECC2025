%    THREAD_ROOTLOCI Comptutes the root locus of a SISO delay system given
%    by its transfer function G :
%  
%           t(s) + SUM(t_i(s)*exp(-i*tau*s))
%    G(s) = ----------------------------------    i=1:N' , j=1:N
%           p(s) + SUM(q_j(s)*exp(-j*tau*s))
%  
%    The delay system must be of :
%    - retarded type
%    - neutral type with a finite number of poles in { Re s > -a, a>0}.
% 
%    Syntax:
%    T = thread_RootLoci(iPolyMatrix,iDelayVect,iAlpha,iTau,ibPlotOption,iDeltaTau)
%    T = thread_RootLoci(iPolyMatrix,iDelayVect,iAlpha,iTau,ibPlotOption)
%    T = thread_RootLoci(iPolyMatrix,iDelayVect,iAlpha,iTau)
%  
%    Inputs: 
%    - iPolyMatrix: the quasi polynomial
%           p(s) + q_1(s)*exp(-tau*s) + ... + q_N(s)*exp(-N*tau*s)
%    - iDelayVect: the delay vector, which is a vector of values of k for
%      which q_(k) is not null of length N (p(s) is assumed non zero).
%    - iAlpha: a real between 0 and 1 describing the fractionnal power alpha
%      for the fractionnary equation in s^(alpha).
%    - iTau: the value of the nominal delay tau.
%    - ibPlotOption: root locus plot option (1 if plot, 0 otherwise), root
%      locus for delays from 0 to iTau, default value = 1
%    - iDeltaTau: precision of integration procedure for computing root
%      locus, default value = 1e-4
%  
%    Output: Structure giving information on the:
%    - roots of the system without delay
%  
%    And for retarded systems or neutral systems with a finite number of
%    unstable poles in { Re s > -a, a>0}:
%    - crossing table
%    - set of all imaginary roots for a delay between 0 and the nominal delay
%      tau.
%    - position of unstable poles where the delay is equal to tau
%    - error on the unstable poles
%    - root locus
%