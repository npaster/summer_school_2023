function ialpha=m2o(alpha)
%
%IALPHA=M2O(ALPHA)
%
% Multi to One
%
% input:  ALPHA = multiindex
% output: IALPHA = linear index
%
%
alpha1 = alpha(1);
alpha2 = alpha(2);
%
kmax=100;
%
ialpha=1;
%
for beta1=0:kmax
    %
    for beta2=beta1:-1:0
        %
        if alpha1==beta2 && alpha2==beta1-beta2
            return
        end
        %
        ialpha = ialpha+1;
        %
    end
    %
end
%
fprintf('%s: error in m2o',mfilename);
%
end



