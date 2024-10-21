function [V,S_V,Crt,W,S_W,Brt] = TangentialKrylov(A,B,C,E,s0,t_B,t_C)
% Tangential Input and Output Rational Krylov Subspaces
%
% Input:
% A,B,C,E: HFM matrices;
% s0: vector of shifts;
% t_B,t_C: matrices of tangential directions as column/row vectors
%
% Output: 
%    A*V - E*V*S_V - B*Crt = 0,
%    W.'*A - S_W*W.'*E - Brt*C = 0
%
% All rights reserved. (c)2014 Heiko K. F. Panzer,Techn. Univ. Muenchen.
% This file is published under the BSD3-Clause License.

% initialization and preallocation
N=size(A,1); n = length(s0); m = size(B,2); p = size(C,1); i=1;
V=zeros(N,n); S_V = zeros(n,n); Crt = zeros(m,n); W = V; S_W=S_V; Brt = zeros(n,p);
while i<=n
    s = s0(i);
    % compute new basis vectors
    [L,U,P,Q,R] = lu(sparse(A-s*E)); % ==> P*(R\A)*Q = L*U
    tempV = Q*(U\(L\(P*(R\(B*t_B(:,i))))));
    tempW = (t_C(i,:)*C*Q/U/L*P/R).';
    if ~isreal(s)
        % complex conjugated pair of shifts -> two columns
        V(:,i:(i+1)) = [real(tempV), imag(tempV)];
        Crt(:,i:(i+1)) = [real(t_B(:,i)), imag(t_B(:,i))];
        S_V(i:(i+1),i:(i+1)) = [real(s), imag(s); -imag(s), real(s)];
        W(:,i:(i+1)) = [real(tempW), imag(tempW)];
        Brt(i:(i+1),:) = [real(t_C(i,:)); imag(t_C(i,:))];
        S_W(i:(i+1),i:(i+1)) = [real(s), -imag(s); imag(s), real(s)];
        i = i+2;
    else
        % real shift -> one column
        V(:,i) = real(tempV);
        Crt(:,i) = real(t_B(:,i)); S_V(i,i) = s;
        W(:,i) = real(tempW);
        Brt(i,:) = real(t_C(i,:)); S_W(i,i) = s;
        i = i+1;
    end
end
% orthogonalization
[V,S_V,Crt] = GramSchmidt(V,S_V,Crt,[1 n]);
[W,S_W,Brt] = GramSchmidt(W,S_W.',Brt.',[1 n]); S_W = S_W.'; Brt = Brt.';

