function W = conj_poly(W)
    for k=1:length(W)
        W(length(W)+1-k) = (-1)^(mod(k+1,2))*W(length(W)+1-k);
    end
end