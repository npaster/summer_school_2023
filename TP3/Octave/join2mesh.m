function mesh=join2mesh(mesh1,mesh2,x1,y1,x2,y2)
%
% MESH=JOIN2MESH(MESH1,MESH2,X1,Y1,X2,Y2)
%
%   join meshes MESH1 and MESH2 along vertical and horizontal lines
%
%   vertical lines are specified by x=X1 and x=X2
%   horizontal lines are specified by y=Y1 and y=Y2
%
%   argument x1 is mandatory
%
%   to join along horizontal line y=Y1 only, set X1=NaN
%
VERBOSE = false;
%
mesh1c = mesh1;
mesh2c = mesh2;
%
switch nargin
    case 3
        y1 = NaN;
        x2 = NaN;
        y2 = NaN;
    case 4
        x2 = NaN;
        y2 = NaN;
    case 5
        y2 = NaN;
end
%
for iv1=1:mesh1.NV
    if ~isnan(x1)
        if iseq(mesh1.vertex(iv1).x,x1)
            xa = mesh1.vertex(iv1).x;
            ya = mesh1.vertex(iv1).y;
            mesh2c = addvertextoedge(xa,ya,mesh2c);
        end
    end
    if ~isnan(x2)
        if iseq(mesh1.vertex(iv1).x,x2)
            xa = mesh1.vertex(iv1).x;
            ya = mesh1.vertex(iv1).y;
            mesh2c = addvertextoedge(xa,ya,mesh2c);
        end
    end
    if ~isnan(y1)
        if iseq(mesh1.vertex(iv1).y,y1)
            xa = mesh1.vertex(iv1).x;
            ya = mesh1.vertex(iv1).y;
            mesh2c = addvertextoedge(xa,ya,mesh2c);
        end
    end
    if ~isnan(y2)
        if iseq(mesh1.vertex(iv1).y,y2)
            xa = mesh1.vertex(iv1).x;
            ya = mesh1.vertex(iv1).y;
            mesh2c = addvertextoedge(xa,ya,mesh2c);
        end
    end
end
%
for iv2=1:mesh2.NV
    if ~isnan(x1)
        if iseq(mesh2.vertex(iv2).x,x1)
            xa = mesh2.vertex(iv2).x;
            ya = mesh2.vertex(iv2).y;
            mesh1c = addvertextoedge(xa,ya,mesh1c);
        end
    end
    if ~isnan(x2)
        if iseq(mesh2.vertex(iv2).x,x2)
            xa = mesh2.vertex(iv2).x;
            ya = mesh2.vertex(iv2).y;
            mesh1c = addvertextoedge(xa,ya,mesh1c);
        end
    end
    if ~isnan(y1)
        if iseq(mesh2.vertex(iv2).y,y1)
            xa = mesh2.vertex(iv2).x;
            ya = mesh2.vertex(iv2).y;
            mesh1c = addvertextoedge(xa,ya,mesh1c);
        end
    end
    if ~isnan(y2)
        if iseq(mesh2.vertex(iv2).y,y2)
            xa = mesh2.vertex(iv2).x;
            ya = mesh2.vertex(iv2).y;
            mesh1c = addvertextoedge(xa,ya,mesh1c);
        end
    end
end

% join the two meshes mesh1c and mesh2c

mesh = mesh1c;
                
for ip2=1:mesh2c.NP
    mesh.polygon(mesh.NP+ip2) = mesh2c.polygon(ip2);
    mesh.polygon(mesh.NP+ip2).vertices = [mesh.polygon(mesh.NP+ip2).vertices] + mesh.NV;
end

for iv2=1:mesh2c.NV
    mesh.vertex(mesh.NV+iv2) = mesh2c.vertex(iv2);
end
                
mesh.NP = mesh.NP + mesh2c.NP;
mesh.NV = mesh.NV + mesh2c.NV;
                    
% find duplicate vertices

TOL = 1e-12;

found = true;
ivstart = 1;

while found
    
    found = false;
    
    for iv=ivstart:mesh.NV
        for jv=1:mesh.NV
            if iv~=jv
                xiv = mesh.vertex(iv).x;
                yiv = mesh.vertex(iv).y;
                xjv = mesh.vertex(jv).x;
                yjv = mesh.vertex(jv).y;
                %
                if ((xiv-xjv)^2+(yiv-yjv)^2)<= TOL^2
                    found = true;
                    break
                end
                if found
                    break
                end
            end
            if found
                break
            end
        end
        if found
            break
        end
    end
    if found
        %
        ivs = iv;
        jvs = jv;
        %
        % elimino il piu' grande
        %
        % se ivs e' piu' piccolo di jvs li scambio
        %
        if ivs<jvs
            tmp = ivs;
            ivs = jvs;
            jvs = tmp;
        end
        %
        % ivs ora e' il piu' grande
        %
        if VERBOSE
            fprintf('%d and %d are the same vertex\n',ivs,jvs)
        end
        %
        % poligoni
        for ip=1:mesh.NP
            vv = [mesh.polygon(ip).vertices];
            ivsind = find(vv==ivs);
            if ~isempty(ivsind)
                vv(ivsind) = jvs;
            end
            ivsind = find(vv>ivs);
            vv(ivsind)=vv(ivsind)-1;
            mesh.polygon(ip).vertices = vv;
        end
        %
        % elimino ivs
        %
        for kv=ivs:mesh.NV-1
            mesh.vertex(kv) = mesh.vertex(kv+1);
        end
        mesh.vertex(mesh.NV) = [];
        mesh.NV = mesh.NV-1;
        
        ivstart = jvs;
        
    end
end
    
mesh.name = [mesh1.name '+' mesh2.name];

mesh = rebuildmesh(mesh);

%fprintf('saving mesh %s\n\n',mesh.name)

%save(mesh.name,'mesh')

return


function meshout=addvertextoedge(xa,ya,mesh)

meshout = mesh;

VERBOSE = false;

% vertex ccordinates are (xa,ya)

% exclude exisiting vertex

for iv=1:mesh.NV
    %
    xiv = mesh.vertex(iv).x;
    yiv = mesh.vertex(iv).y;
    %
    if iseq(xa,xiv) && iseq(ya,yiv)
        if VERBOSE
            fprintf('%s: vertex is already present!\n',mfilename);
        end
        return
    end
end

% look for edge

for ie=1:mesh.NE
    %
    x1 = mesh.vertex(mesh.edge(ie).v1).x;
    y1 = mesh.vertex(mesh.edge(ie).v1).y;
    %
    x2 = mesh.vertex(mesh.edge(ie).v2).x;
    y2 = mesh.vertex(mesh.edge(ie).v2).y;
    %
    sx = (x2-xa)/(x2-x1);
    sy = (y2-ya)/(y2-y1);
    %
        
    if (iseq(x1,x2) && iseq(x1,xa) && 0<= sy && sy <= 1) || ...
       (iseq(y1,y2) && iseq(y1,ya) && 0<= sx && sx <= 1) || ...
       (iseq(sx,sy) && 0<= sx && sx <= 1)
        
        if VERBOSE
            fprintf('%s: found edge %d\n',mfilename,ie);
        end
        
        % now we find the polygon
        for ip=1:mesh.NP
            iel = find([mesh.polygon(ip).edges]==ie);
            %
            if length(iel)==1
                %
                % add new vertex
                %
                vv = mesh.polygon(ip).vertices;
                
                meshout.polygon(ip).vertices = [vv(1:iel) meshout.NV+1 vv(iel+1:end)];

                meshout.NV = meshout.NV+1;
                meshout.vertex(meshout.NV).x = xa;
                meshout.vertex(meshout.NV).y = ya;
                
                meshout = rebuildmesh(meshout);

                return
                
            end
        end
    end
end

function iseq=iseq(x,y)

TOL = 1e-12;

if abs(x-y)<TOL
    iseq = true;
else
    iseq = false;
end

return
