function tensorC=computeElasticityTensorC(lambda, mu, isVoigtNotation)
if nargin==2
    isVoigtNotation=false;
end
%
%   |2*mu+lambda   lambda      0 |
% C=| lambda     2*mu+lambda   0 |
%   |   0            0       2*mu| 
tensorC=zeros(3,3);
% block A
% diagonal elements
tensorC(1,1)=2*mu+lambda;
tensorC(2,2)=2*mu+lambda;
%
% If we use voigt notation the term C(2,2) is 4*mu and not 2*mu
%
if(isVoigtNotation) 
    tensorC(3,3)=4*mu; 
else
    tensorC(3,3)=2*mu;
end
% Upper diagonal part
tensorC(1,2)=lambda;
% Lower diagonal part
tensorC(2,1)=lambda;
end

