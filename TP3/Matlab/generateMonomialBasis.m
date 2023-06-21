function [mono, dimPk] = generateMonomialBasis(k,ip,mesh)
%
%%% dimension of Pk
%
dimPk = (k+1)*(k+2)/2;
%
% centroid of the polygon
%
xb = mesh.polygon(ip).xb;
yb = mesh.polygon(ip).yb;
%
% diameter of the polygon
%
h = mesh.polygon(ip).h;
%
%
X = @(x) (x-xb)/h;
Y = @(y) (y-yb)/h;
%
% iniziatization
%
mono = struct;
%
% construction of local basis
%
for ialpha=1:dimPk
    % alpha
    alpha = o2m(ialpha);
    % mono
    mono(ialpha).p = @(x,y) X(x).^alpha(1) .* Y(y).^alpha(2);
    mono(ialpha).indexX=alpha(1);
    mono(ialpha).indexY=alpha(2);
    mono(ialpha).coeff=1.;
end

