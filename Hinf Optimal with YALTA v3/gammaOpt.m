% % Computes optimal performance level by using iterative algorithm 
function [gamma_min,xvalue,yvalue,newFunc,beta_fin]=gammaOpt(nw1,dw1,nw2,dw2,alpha,Npoints,gmin,gmax,EPS,Mn)
rootsEPS = 1e-4;
minValCoeff=0.1;
iterNu = 0;
ind = 1;
while 1
    xval=[];
    yval=[];
    iterNu = iterNu+1;
    minValCoeff=minValCoeff*0.05;
    for gamma=gmin:(gmax-gmin)/Npoints:gmax,
    a1=polyadd(conv(nw1,star(nw1)),-gamma^2*conv(dw1,star(dw1)));
    a1 = a1(find(a1,1):end);
    
    b1=gamma^2*conv(dw1,star(dw1)); 
    b1 = b1(find(b1,1):end);

    a2=polyadd(conv(nw2,star(nw2)),-gamma^2*conv(dw2,star(dw2)));
    a2 = a2(find(a2,1):end);
    b2=gamma^2*conv(dw2,star(dw2));
    b2 = b2(find(b2,1):end);
    A1=conv(a1,a2);
    B1=conv(b1,b2);
    A=polyadd(B1,-A1);
    B=B1;
    [dG,nG]=spec(A,B); 
    eta=roots(dw1);
    n1=length(eta);
    nF=conv(nG,poly(-eta));
    dF=conv(dG,poly(eta));

    beta=roots(a1);
    beta=beta(find((real(beta)<rootsEPS) | ((abs(real(beta))<=rootsEPS) & (imag(beta)>=0))));
    l=length(alpha);
    if abs(gamma-(gmin+gmax)/2)<1e-6
        beta_fin = beta;
    end
    betas = beta;
       M=[];
       for k=1:l,
           s=alpha(k);
           M=[M; alpha(k).^[0:n1+l-1] eval(Mn)*polyval(nF,alpha(k)) ...
              /polyval(dF,alpha(k))*alpha(k).^[0:n1+l-1]];
       end;
       for k=1:n1,
           s=beta(k);
           M=[M; beta(k).^[0:n1+l-1] eval(Mn)*polyval(nF,beta(k)) ...
              /polyval(dF,beta(k))*beta(k).^[0:n1+l-1]];
       end;
       for k=1:l,
           s=alpha(k);
           M=[M;(-alpha(k)).^[0:n1+l-1]*eval(Mn)*polyval(nF,alpha(k)) ...
              /polyval(dF,alpha(k)) (-alpha(k)).^[0:n1+l-1]];
       end;

       for k=1:n1,
           s=beta(k);
           M=[M; (-beta(k)).^[0:n1+l-1]*eval(Mn)*polyval(nF,beta(k)) ...
              /polyval(dF,beta(k)) (-beta(k)).^[0:n1+l-1]];
       end;
       assignin('base','M',M);
       xval=[xval gamma];
       try
            yval=[yval min(svd(M))];
       catch error
           error
           break;
       end

    end;
    
    yval=fliplr(yval);
    xval=fliplr(xval);
    if iterNu==1
        xvalue = xval;
        yvalue=abs(yval);
        plot(xvalue,yvalue)
        newFunc=-yval;
        figure;
        plot(xval,newFunc)
        [~,b]=findpeaks(newFunc);
        gamma_opt=xval(b);
        try
            minValCoeff=0.5*(gamma_opt(ind:ind+1)*[1;-1]);
        catch
            if gamma_opt(ind)>1e-2
                minValCoeff = 1e-2;
            else
                minValCoeff = 3*gamma_pos/4;
            end
        end
        gamma_pos = gamma_opt(ind);
        temp = yval(b);
    else
        [~,b]=min(yval);
        gamma_pos=xval(b);
    end
    if yval(b)<EPS
        Jmatrix = zeros(1,n1);
        for ind=1:n1 
            Jmatrix(ind) = (-1)^(ind+1);
        end
        Jmatrix = diag(Jmatrix); s = beta(1:n1);
        Dn1 = diag(eval(Mn).*polyval(nF,beta(1:n1))./polyval(nF,beta(1:n1)));
%         if min(abs(Dn1-Jmatrix))<1e-5 || min(abs(Dn1+Jmatrix))<1e-5
%             iterNu = 1;
%             ind = ind+1;
%             gamma_pos = gamma_opt(ind);
%             try
%                 minValCoeff=0.5*(gamma_opt(ind:ind+1)*[1;-1]);
%             catch
%                 if gamma_opt(ind)>1e-2
%                     minValCoeff = 1e-2;
%                 else
%                     minValCoeff = 3*gamma_pos/4;
%                 end
%             end
%             temp = yval(b);min(temp)
%             gmin=gamma_pos-minValCoeff
%             gmax=gamma_pos+minValCoeff
%         else
            break;
%         end
    else
        if iterNu>5 && yval(b)>(4*min(temp)/5)
            iterNu = 1;
            ind = ind+1;
            gamma_pos = gamma_opt(ind);
            try
                minValCoeff=0.5*(gamma_opt(ind:ind+1)*[1;-1]);
            catch
                if gamma_opt(ind)>1e-2
                    minValCoeff = 1e-2;
                else
                    minValCoeff = 3*gamma_pos/4;
                end
            end
        end
        temp = yval(b);min(temp)
        gmin=gamma_pos-minValCoeff
        gmax=gamma_pos+minValCoeff
    end
end
gamma_min = gamma_pos;