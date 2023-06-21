function [grad_eps, dim] = computeSymmetricGradientForMonomialBasisVect(k,ip,mesh)
%
%%% compute Derivative 
%
[grad,dimPk]=computeGradientForMonomialBasisVect(k,ip,mesh);
%
%%% dimension
%
dim=dimPk;
%
% iniziatization
%
grad_eps = struct;
%
% construction of local basis
%
for ialpha=1:dimPk
    grad_eps(ialpha).p = @(x,y) (grad(ialpha).p(x,y)+grad(ialpha).p(x,y)')/2;
end
   
end

