function [nump,denp] = spectral_factorization(nums,dens) 

numK=nums(1)/dens(1)*(-1)^((length(nums)-length(dens))/2);
% numK=nums(1)/dens(1)*(-1)^((length(nums)-length(dens)));
if numK<0,
   error='Cannot do spec-fact.'
   nump=[];
   denp=[];
%end
else
    numK=sqrt(numK);

    nump=numK*pspec(nums);
    denp=pspec(dens);

    nump=real(nump);
    denp=real(denp);
end;

