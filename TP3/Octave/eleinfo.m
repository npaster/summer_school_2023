function eleinfo(ip,mesh)
% 
% ELEINFO(IP,MESH)
%
% draws all information about polygon IP of mesh MESH
%
if ip<1 || ip>mesh.NP
    fprintf('%s: ip out of range\n',mfilename);
    return
end
%
% drawmesh(mesh,ip,'yellow');
drawmesh(mesh);
%
numele(mesh,ip);
%
numedge(mesh,mesh.polygon(ip).edges);
%
numver(mesh,mesh.polygon(ip).vertices);
%
title(sprintf('polygon %d, h=%0.2e, area=%0.4e',ip,mesh.polygon(ip).h,mesh.polygon(ip).area));
%
axis normal
axis tight
axis equal
zlim([-0.1*mesh.polygon(ip).h, 0.1*mesh.polygon(ip).h])
%axis vis3d
%
end
