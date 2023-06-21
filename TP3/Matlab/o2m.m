function alpha=o2m(ialpha)
%
%ALPHA=02M(IALPHA)
%
% One to Multi
%
% input:  IALPHA = linear index
% output: ALPHA = multiindex
%

%
kmax=100;
%
ibeta = 1;
%
for beta1=0:kmax
    %
    for beta2=beta1:-1:0
        %
        if ibeta==ialpha
            %
            alpha1 = beta2;
            alpha2 = beta1-beta2;
            %
            alpha = [alpha1,alpha2];
            %
            return
        end
        %
        ibeta = ibeta+1;
        %
    end
    %
end
%
fprintf('%s: error in o2m',mfilename);
%
end



