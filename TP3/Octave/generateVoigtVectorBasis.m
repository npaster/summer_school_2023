function [basis, dim] = generateVoigtVectorBasis(k,ip,mesh)
%
%%% function to generate monomial basis
%
[mono, dimPk] = generateMonomialBasis(k,ip,mesh);
%
% dimension for Voigt vectors
%
dim = 3 * dimPk;
%
% iniziatization
%
basis=struct;
%
% null function
%
nullFun=@(x,y) 0*x.*y;
%
% construction of local basis
%
for ialpha=1:dimPk
    basis(ialpha).p =@(x,y) [mono(ialpha).p(x,y); nullFun(x,y); nullFun(x,y)];
    basis(ialpha+dimPk).p =@(x,y) [nullFun(x,y); mono(ialpha).p(x,y); nullFun(x,y)];
    basis(ialpha+2*dimPk).p =@(x,y) [nullFun(x,y); nullFun(x,y); mono(ialpha).p(x,y)];
end

end