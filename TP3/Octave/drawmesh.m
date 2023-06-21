function drawmesh(mesh,iprange,colmesh)
%
% DRAWMESH(MESH,IPRANGE,COLMESH)
%
% Draws the mesh MESH.
%
% If IPRANGE and COLMESH are missing, draws all polygons transparent;
% if COLMESH is missing, draws polygons in range IPRANGE transparent;
% otherwise, draws polygons in range IPRANGE of color COLMESH.
% IF IPRANGE==0, draws all polygons.
%
switch nargin
    case 1
        iprange = 1:mesh.NP;
        colmesh = 'none';
    case 2
        colmesh = 'none';
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
Faces = NaN*ones(mesh.NP,max([mesh.polygon.NV]));
%
for ip=iprange
    Faces(ip,1:mesh.polygon(ip).NV) = [mesh.polygon(ip).vertices]';    
end
%
Xv = [mesh.vertex.x]';
Yv = [mesh.vertex.y]';
Zv = zeros(size(Xv));
%
patch('Vertices',[Xv Yv],'Faces',Faces,'FaceColor',colmesh)
%
hmax = max([mesh.polygon.h]);
hmean = sum([mesh.polygon.h])/mesh.NP;
%
% print mesh parameters
%
stitle = sprintf('%s\n%d polygons, %d edges, %d vertices',mesh.name,mesh.NP,mesh.NE,mesh.NV);
stitle = [stitle sprintf('\nhmax = %0.2e, hmean = %0.2e',hmax,hmean)];
%
% if k is defined, print also k and number of dofs
%
if isfield(mesh,'EZVEM')
    if isfield(mesh.EZVEM,'K')
        %
        k = mesh.EZVEM.K;
        %
        Ndof_global = mesh.NV + (k-1)*mesh.NE + dimP(k-2)*mesh.NP;
        %
        stitle = [stitle sprintf(', k=%d, Ndof=%d',k,Ndof_global)];
        %
    end
end
if isfield(mesh,'ELASTICITYVEM')
    if isfield(mesh.ELASTICITYVEM,'K')
        %
        k = mesh.ELASTICITYVEM.K;
        %
        Ndof_global = 2*(mesh.NV + (k-1)*mesh.NE + dimP(k-2)*mesh.NP);
        %
        stitle = [stitle sprintf(', k=%d, Ndof=%d',k,Ndof_global)];
        %
    end
end
%
h = title(stitle);
%
axis tight
axis equal
zlim([-0.1 0.1])
%axis vis3d
%
end
