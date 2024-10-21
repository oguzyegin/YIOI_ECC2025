%    DELAYFREQUENCYANALYSISMIN Main function of the toolbox. It launches
%    system analysis, computes and displays stability windows and root locus
%    of a SISO delay system with transfer function G of the type :
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
%    T = delayFrequencyAnalysisMin(iPolyMatrix,iDelayVect,iAlpha,iTau,ibPlotOption,iDeltaTau,tminsw,tmaxsw)
%    T = delayFrequencyAnalysisMin(iPolyMatrix,iDelayVect,iAlpha,iTau,ibPlotOption,iDeltaTau)
%    T = delayFrequencyAnalysisMin(iPolyMatrix,iDelayVect,iAlpha,iTau,ibPlotOption)
%    T = delayFrequencyAnalysisMin(iPolyMatrix,iDelayVect,iAlpha,iTau)
%  
%    Inputs: 
%    - iPolyMatrix: the quasi polynomial
%           p(s) + q_1(s)*exp(-tau*s) + ... + q_N(s)*exp(-N*tau*s)
%    - iDelayVect: the delay vector, which is a vector of values of k for
%    which q_(k) is not null of length N (p(s) is assumed non zero).
%    - iAlpha: a real between 0 and 1 describing the fractionnal power alpha
%    for the fractionnary equation in s^(alpha).
%    - iTau: the value of the nominal delay tau.
%    - ibPlotOption: root locus plot option (1 if plot, 0 otherwise), root
%    locus for delays from 0 to iTau, default value = 1
%    - iDeltaTau: precision of intergration procedure for computing root
%    locus, default value = 1e-4
%    - tminsw: minimum delay of stability window, default value = 0
%    - tmaxsw: maximum delay of stability window, default value = iTau
% 
%    Example: Please refer to the documentation UserDocYalta.pdf to see a
%    example for function delayFrequencyAnalysisMin.
%  
%    Output: Structure giving information on the: 
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
%    - set of all imaginary roots for a delay between 0 and the delay tau
%    - stability window for a delay between 0 and the maximum delay tmaxsw
%    - value of the delays for which the number of unstable poles changes
%    - number of unstable poles  when the delay is equal to tau
%    - position of unstable poles when the delay is equal to tau
%    - error on the unstable poles
%    - root locus.
%