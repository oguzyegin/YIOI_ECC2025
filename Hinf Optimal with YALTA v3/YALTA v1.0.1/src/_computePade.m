%    COMPUTEPADE Computes the Pade-2 approximation of a non fractionnal
%    SISO delay system with transfer function G of the type :
%  
%              p(s) + SUM_j(q_j(s)*exp(-j*tau*s))
%    C(s) =  -----------------------------------   j=1:N
%                    (s + 1) ^ iDelta
%  
%    with iDelta >= (n + 1) where n is the degree of p(s). 
% 
%    There are two modes for computing the Pade approximation (iMode argument).
%    iMode is either 'ORDER' and thus iModArg is the order of the
%    approximation, or iMode is 'NORM' and iModArg represents the maximum
%    limit of the difference (H-infinity-norm) between C(s) and its
%    approximation.
% 
%    Output oPadeStruct is a compact structure with following informations:
%    - NumApprox : vector of coefficients of the numerator of the transfer
%      function of the approximation.
%    - DenApprox : vector of coefficients of the denominator of the transfer
%      function of the approximation.
%    - ErrorNorm : the H-infinity-norm of the difference between C(s) and
%      its approximation.
%    - PadeOrder : the order of the approximation.
%    - Roots : the roots of the approximation.
%    - RootsError : an array of differences between the computed roots of
%      C(s) (by the function thread_RootLoci) and the roots of the
%      approximation.
% 
%    Syntax:
%    oPadeStruct = COMPUTEPADE(iPolyMatrix, iDelta, iTau, iDelayVector, ...
%                              iModArg, iMode)
%  
%    Example:
%    If iPolyMatrix = [6, -66, 180;
%                      0,  -2, 12;
%                      0,  6,  30;
%                      0,  0,  -2];
%  
%    and iDelayVector = [1, 2, 3];
%  
%    then oPadeStruct = computePade(iPolyMatrix, 3, 1, iDelayVector, 4, 'ORDER')
%    ans =
%  
%    NumApprox: [1x37 double]
%    DenApprox: [1x37 double]
%    ErrorNorm: 8.7100e-04
%    PadeOrder: 4
%    Roots: [5.9996 5.0032]
%    RootsError: [2.9745e-04 4.9066e-04]
%