%    THREAD_STABILITYWINDOWS Comptutes the stability windows of the system
%    given by its transfer function G :
%  
%           t(s) + SUM(t_i(s)*exp(-i*tau*s))
%    G(s) = ----------------------------------    i=1:N' , j=1:N
%           p(s) + SUM(q_j(s)*exp(-j*tau*s))
% 
%    The delay system must be of :
%    - retarded type.
%    - neutral type with a finite number of poles in { Re s > -a, a>0}.
%  
%    Syntax:
%    T = thread_StabilityWindows(iPolyMatrix,iDelayVect,iAlpha,iTau,tminsw,tmaxsw,plot)
%    T = thread_StabilityWindows(iPolyMatrix,iDelayVect,iAlpha,iTau,tminsw,tmaxsw)
%    T = thread_StabilityWindows(iPolyMatrix,iDelayVect,iAlpha,iTau,tmaxsw)
%    T = thread_StabilityWindows(iPolyMatrix,iDelayVect,iAlpha,iTau)
%  
%    Inputs: 
%    - iPolyMatrix: the quasi polynomial
%           p(s) + q_1(s)*exp(-tau*s) + ... + q_N(s)*exp(-N*tau*s)
%    - iDelayVect: the delay vector, which is a vector of values of k for
%      which q_(k) is not null of length N (p(s) is assumed non zero).
%    - iAlpha: a real between 0 and 1 describing the fractionnal power alpha
%      for the fractionnary equation in s^(alpha).
%    - iTau: the value of the nominal delay tau.
%    - tminsw: minimum delay of stability window, default value = 0
%    - tmaxsw: maximum delay of stability window, default value = iTau
%    - plot: option if we want to plot the stabitlity windows.
%  
%    Output : Structure giving information on the :
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
%    - stability window for a delay between the minimum delay tminsw and the 
%      maximum delay tmaxsw
%    - value of the delays for which the number of unstable poles changes
%