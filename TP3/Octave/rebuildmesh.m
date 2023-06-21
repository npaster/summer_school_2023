function meshout=rebuildmesh(mesh)
%
% rebuilds a mesh starting from minimal information:
% mesh.vertex.x, mesh.vertex.y
% mesh.polygon.vertices
%

%
% find maximum number of vertices in polygons
%
NP = length(mesh.polygon);
%
NVmax = 0;
for ip=1:NP
    NV = length(mesh.polygon(ip).vertices);
    if NV>NVmax
        NVmax = NV;
    end
end
%
nodi = zeros(NP,NVmax);
%
for ip=1:NP
    NV = length(mesh.polygon(ip).vertices);
    if NV>0
        nodi(ip,1:NV)=mesh.polygon(ip).vertices;
    end
end
%
meshout = conn2mesh([mesh.vertex.x],[mesh.vertex.y],nodi);
%
if isfield(mesh,'name')
    meshout.name = mesh.name;
else
    meshout.name = '';
end

return
