function drawcontour(mesh,NL)
%
if nargin==1
    NL = 20;
end
%
N = 100;
[X,Y] = meshgrid(linspace(0,1,N+1),linspace(0,1,N));
%
k = mesh.EZVEM.K;
%
P0kuh = zeros(size(X));
%
for ip=1:mesh.NP
    %
    % global dofs of the polygon
    %
    global_dofs = get_global_dofs(k,ip,mesh);
    %
    % coefficients of m_alpha of the two projections of uh
    %
    sP0k = mesh.polygon(ip).P0kstar*mesh.uh(global_dofs);
    sPNk = mesh.polygon(ip).PNkstar*mesh.uh(global_dofs);
    %
    xv = [mesh.vertex([mesh.polygon(ip).vertices]).x];
    yv = [mesh.vertex([mesh.polygon(ip).vertices]).y];
    %
    [rin,cin] = find(inpolygon(X,Y,xv,yv));
    %
    Xin = X(rin,cin);
    Yin = Y(rin,cin);
    %
    Zin = 0;
    %
    for ialpha=1:dimP(k)
        Zin = Zin + sPNk(ialpha)*m(ialpha,Xin,Yin,ip,mesh);
    end
    %
    P0kuh(rin,cin) = Zin;
end
%contourf(X,Y,P0kuh,NL,'LineWidth',2)
contour(X,Y,P0kuh,NL,'LineWidth',2)
%
