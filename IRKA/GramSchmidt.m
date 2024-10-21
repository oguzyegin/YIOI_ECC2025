function [X,Y,Z] = GramSchmidt(X,Y,Z,cols)
%  function [X,Y,Z] = GramSchmidt(X,Y,Z,cols)
% Gram-Schmidt orthonormalization
% Input:    X,[Y,[Z]]:  matrices in Sylvester eq.: V,S_V,Crt or W.',S_W.',Brt.'
%           cols:       2-dim. vector: number of first and last column to be treated
% Output:   X,[Y,[Z]]: solution of Sylvester eq. with X.'*X = I
%
% This code is derived from the file GramSchmidt.m [Panzer14] published
% under the BSD3-Clause License.
% All rights reserved. (c) 2014
% Code used for the lecture: Modelreduction of Mechanical Systems
% Institut fuer Technische und Numerische Mechanik, Uni Stuttgart
% Profs Eberhard / Fehr / Hanss
% SS 2015
% Jun. Prof. J. Fehr    Dipl.-Math. D. Grunert
if nargin<4, cols=[1 size(X,2)]; end
for k=cols(1):cols(2)
    for j=1:(k-1)
        % orthogonalization
        T = eye(size(X,2)); T(j,k)=-X(:,k)'*X(:,j);
        X = X*T;
        if nargout>=2, Y = T\Y*T; end
        if nargout>=3, Z = Z*T; end
    end
    h = norm(X(:,k));
    X(:,k) = X(:,k)/h; % normalization
    if nargout>=2, Y(:,k) = Y(:,k)/h; Y(k,:) = Y(k,:)*h; end
    if nargout==3, Z(:,k) = Z(:,k)/h; end
end

