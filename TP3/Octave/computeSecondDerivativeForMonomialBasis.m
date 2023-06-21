function [mono_xx, mono_yy, mono_xy, dimPk] = computeSecondDerivativeForMonomialBasis(k,ip,mesh)
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
mono_xx = struct;
mono_xy = struct;
mono_yy = struct;
%
% construction of local basis
%
for ialpha=1:dimPk
    % alpha
    alpha = o2m(ialpha);
    % mono_xx
    if alpha(1)<=1
        mono_xx(ialpha).p=@(x,y) 0;
        mono_xx(ialpha).indexX=alpha(1)-2;
        mono_xx(ialpha).indexY=alpha(2);
        mono_xx(ialpha).coeff = 0.;
    else
        mono_xx(ialpha).p = @(x,y) alpha(1)*(alpha(1)-1)*X(x).^(alpha(1)-2)*(1/h)*(1/h).* Y(y).^alpha(2);
        mono_xx(ialpha).indexX=alpha(1)-2;
        mono_xx(ialpha).indexY=alpha(2);
        mono_xx(ialpha).coeff = alpha(1)*(alpha(1)-1)*(1/h)*(1/h);
    end
    % mono_yy
    if alpha(2)<=1
        mono_yy(ialpha).p=@(x,y) 0;
        mono_yy(ialpha).indexX=alpha(1);
        mono_yy(ialpha).indexY=alpha(2)-2;
        mono_yy(ialpha).coeff = 0;
    else
        mono_yy(ialpha).p=@(x,y) X(x).^alpha(1).*alpha(2)*(alpha(2)-1).*Y(y).^(alpha(2)-2)*(1/h)*(1/h); 
        mono_yy(ialpha).indexX=alpha(1);
        mono_yy(ialpha).indexY=alpha(2)-2;
        mono_yy(ialpha).coeff = alpha(2)*(alpha(2)-1)*(1/h)*(1/h);
    end
    % mono_xy
    if (alpha(1)==0) || (alpha(2)==0)
        mono_xy(ialpha).p=@(x,y) 0;
        mono_xy(ialpha).indexX=alpha(1)-1;
        mono_xy(ialpha).indexY=alpha(2)-1;
        mono_xy(ialpha).coeff = 0;
    else
        mono_xy(ialpha).p=@(x,y) alpha(1)*X(x).^(alpha(1)-1).*alpha(2).*Y(y).^(alpha(2)-1)*(1/h)*(1/h);
        mono_xy(ialpha).indexX=alpha(1)-1;
        mono_xy(ialpha).indexY=alpha(2)-1;
        mono_xy(ialpha).coeff = alpha(1)*alpha(2)*(1/h)*(1/h);
    end
end

end



