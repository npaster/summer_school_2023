function numele(mesh,iprange)
%
% NUMELE(MESH,IPRANGE)
%
% Prints number of polygons on the mesh MESH.
%
% If IPRANGE is missing, numbers all polygons;
% otherwise, numbers polygons in range IPRANGE
%
if nargin==1
    iprange = 1:mesh.NP;
end
%
if isscalar(iprange)
    if iprange==0
        iprange = 1:mesh.NP;
    end
end
%
in_range = find(iprange<=mesh.NP & iprange>=1);
%
if length(in_range)<length(iprange)
    fprintf('%s: iprange out of polygons range; cropping\n',mfilename);
    iprange = iprange(in_range);
end
%
str= '';
for ip=iprange
    %
    str = strvcat(str,int2str(ip));
    %
end
%
text([mesh.polygon(iprange).xb],[mesh.polygon(iprange).yb],1e-3*ones(size(iprange)),...
    str,'Color','r','BackgroundColor','w','HorizontalAlignment','center')
%
end
