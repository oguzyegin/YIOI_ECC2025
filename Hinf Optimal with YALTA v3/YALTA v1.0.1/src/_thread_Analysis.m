%    THREAD_ANALYSIS Launches the analysis of a SISO delay system given by
%    its transfer function G :
%  
%           t(s) + SUM_i(t_i(s)*exp(-i*tau*s))
%    G(s) = ----------------------------------    i=1:N' , j=1:N
%           p(s) + SUM_j(q_j(s)*exp(-j*tau*s))
%  
%    The delay system must be of :
%    - retarded type.
%    - neutral type (but some functions will be available only for neutral systems 
%      with a finite number of poles in { Re s > -a, a>0}).
%  
%    Syntax:
%    T = thread_Analysis(iPolyMatrix, iDelayVect, iAlpha, iTau)
%  
%    Inputs : 
%    - iPolyMatrix: the quasi polynomial
%           p(s) + q_1(s)*exp(-tau*s) + ... + q_N(s)*exp(-N*tau*s)
%    - iDelayVect: the delay vector, which is a vector of values of k for
%      which q_(j) is not null of length N (p(s) is assumed non zero).
%    - iAlpha: a real between 0 and 1 describing the fractionnal power alpha
%      for the fractionnary equation in s^(alpha).
%    - iTau: the value of the nominal delay tau.
%  
%    Output: Structure giving information on the:
%    - asymptotic stability of the system
%    - type of the system: retarded, neutral 
%      or advanced (in case the user was wrong in defining his/her system). 
%    - roots of the system without delay
%    - position of the asymptotic axes of chains of poles (in the case of neutral systems)
%      and if applicable the information that the system has a infinite number of unstable poles 
%    - the position of the chain of poles relative to the asymptotic axis 
%      in the case the asymptotic axis is the imaginary axis.
%  
%    And for retarded systems or neutral systems with a finite number of
%    unstable poles in { Re s > -a, a>0}:
%    - crossing table
%    - set of all imaginary roots for a delay between 0 and the delay tau.
%