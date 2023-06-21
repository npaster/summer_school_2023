function numver(mesh,ivrange)
%
% NUMVER(MESH,IVRANGE)
%
% Prints number of vertices on the mesh MESH.
%
% If IVRANGE is missing, numbers all vertices;
% otherwise, numbers vertices in range IVRANGE
%
if nargin==1
    ivrange = 1:mesh.NV;
end
%
%
if isscalar(ivrange)
    if ivrange==0
        ivrange = 1:mesh.NV;
    end
end
%
in_range = find(ivrange<=mesh.NV & ivrange>=1);
%
if length(in_range)<length(ivrange)
    fprintf('%s: ivrange out of vertices range; cropping\n',mfilename);
    ivrange = ivrange(in_range);
end
%
str= '';
for iv=ivrange
    %
    str = strvcat(str,int2str(iv));
    %
end
%
xv = [mesh.vertex(ivrange).x];
yv = [mesh.vertex(ivrange).y];
%    
%
text(xv,yv,1e-3*ones(size(ivrange)),str,'Color','b','BackgroundColor','w','HorizontalAlignment','center')
%
end
