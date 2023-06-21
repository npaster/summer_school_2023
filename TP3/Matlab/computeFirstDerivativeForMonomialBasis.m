function [mono_x, mono_y, dimPk] = computeFirstDerivativeForMonomialBasis(k,ip,mesh)
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
mono_x =  struct;
mono_y = struct;
%
% construction of local basis
%
for ialpha=1:dimPk
    % alpha
    alpha = o2m(ialpha);
    % mono_x
    if alpha(1)==0
        mono_x(ialpha).p=@(x,y) 0;
        mono_x(ialpha).indexX =alpha(1)-1;
        mono_x(ialpha).indexY =alpha(2);
        mono_x(ialpha).coeff = 0;
    else
        mono_x(ialpha).p=@(x,y) alpha(1)*X(x).^(alpha(1)-1)*(1/h) .* Y(y).^alpha(2);
        mono_x(ialpha).indexX =alpha(1)-1;
        mono_x(ialpha).indexY =alpha(2);
        mono_x(ialpha).coeff = alpha(1)/h;
    end
    % mono_y
    if alpha(2)==0
        mono_y(ialpha).p=@(x,y) 0;
        mono_y(ialpha).indexX =alpha(1);
        mono_y(ialpha).indexY =alpha(2)-1;
        mono_y(ialpha).coeff = 0;
    else
        mono_y(ialpha).p=@(x,y) X(x).^alpha(1) .* alpha(2).*Y(y).^(alpha(2)-1)*(1/h); 
        mono_y(ialpha).indexX = alpha(1);
        mono_y(ialpha).indexY = alpha(2)-1;
        mono_y(ialpha).coeff = alpha(2)/h;
    end
end

end



