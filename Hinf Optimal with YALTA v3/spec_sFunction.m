function[nump,denp]=spec_sFunction(numFunc,denFunc)
nums=sym2poly(numFunc);
dens=sym2poly(denFunc);
numK=nums(1)/dens(1)*(-1)^((length(nums)-length(dens))/2);

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
