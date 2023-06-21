function numedge(mesh,ierange)
%
% NUMEDGE(MESH,IERANGE)
%
% Prints number of edges on the mesh MESH.
%
% If IERANGE is missing, numbers all edges;
% otherwise, numbers edges in range IERANGE
%
if nargin==1
    ierange = 1:mesh.NE;
end
%
%
if isscalar(ierange)
    if ierange==0
        ierange = 1:mesh.NE;
    end
end
%
in_range = find(ierange<=mesh.NE & ierange>=1);
%
if length(in_range)<length(ierange)
    fprintf('%s: ierange out of edges range; cropping\n',mfilename);
    ierange = ierange(in_range);
end
%

str= '';
for ie=ierange
    %
    str = strvcat(str,int2str(ie));
    %
end
%
x1 = [mesh.vertex([mesh.edge(ierange).v1]).x];
y1 = [mesh.vertex([mesh.edge(ierange).v1]).y];
%    
x2 = [mesh.vertex([mesh.edge(ierange).v2]).x];
y2 = [mesh.vertex([mesh.edge(ierange).v2]).y];
%
xe = (x1+x2)/2;
ye = (y1+y2)/2;
%
text(xe,ye,1e-3*ones(size(ierange)),str,'Color','g','BackgroundColor','w','HorizontalAlignment','center')
%
end
