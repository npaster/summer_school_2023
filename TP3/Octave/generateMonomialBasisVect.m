function [monoVect, dim] = generateMonomialBasisVect(k,ip,mesh)
%
%%% function to generate monomial basis
%
[mono, dimPk] = generateMonomialBasis(k,ip,mesh);
%
%%% dimension of vector Pk
%
dim = 2*dimPk;
%
% iniziatization
%
monoVect = struct;
%
% null function
%
nullFun=@(x,y) 0*x.*y;
%
% construction of local basis
%
for ialpha=1:dimPk
    monoVect(ialpha).p = @(x,y) [mono(ialpha).p(x,y); nullFun(x,y)];
    monoVect(ialpha+dimPk).p = @(x,y) [nullFun(x,y); mono(ialpha).p(x,y)]; 
end
   
end

