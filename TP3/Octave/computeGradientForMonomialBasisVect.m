function [grad, dim] = computeGradientForMonomialBasisVect(k,ip,mesh)
%
%%% compute Derivative 
%
[mono_x,mono_y,dimPk]=computeFirstDerivativeForMonomialBasis(k,ip,mesh);
%
%%% dimension
%
dim=2*dimPk;
%
% iniziatization
%
grad = struct;
%
% null function
%
nullFun=@(x,y) 0*x.*y;
%
% construction of local basis
%
for ialpha=1:dimPk
    grad(ialpha).p = @(x,y) [mono_x(ialpha).p(x,y), mono_y(ialpha).p(x,y); nullFun(x,y) , nullFun(x,y)];
    grad(ialpha+dimPk).p = @(x,y) [nullFun(x,y) , nullFun(x,y) ;mono_x(ialpha).p(x,y), mono_y(ialpha).p(x,y)];
end
   
end

