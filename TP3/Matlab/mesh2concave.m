function mesh=mesh2concave(mesh)
%
%MESH=MESH2CONCAVE(MESH)
%
%  transforms a convex mesh into a concave mesh
%  by "carving" existing polygons
%

for ip=1:mesh.NP
    
    % modify polygon
    
    v = [mesh.polygon(ip).vertices];
    e = [mesh.polygon(ip).edges];
    
    el = [mesh.edge(e).length];
    
    [emax, iemax] = max(el);    
    [emin, iemin] = min(el);    
    
    emean = (emax+emin)/2;
    
    tmp = find(el>=emean);
    iepicked = tmp(1);
    
    mesh.polygon(ip).vertices = [v(1:iepicked) mesh.NV+ip v((iepicked+1):end)];
    
    % add polygon
    
    iepickedp = iepicked+1;
    if iepickedp==length(e)+1
        iepickedp = 1;
    end
    
    mesh.polygon(mesh.NP+ip).vertices = [v(iepickedp) mesh.NV+ip v(iepicked)];
    
    % add a vertex    
    
    iv1 = mesh.edge(e(iepicked)).v1;
    iv2 = mesh.edge(e(iepicked)).v2;

    xm = (mesh.vertex(iv1).x+mesh.vertex(iv2).x)/2;
    ym = (mesh.vertex(iv1).y+mesh.vertex(iv2).y)/2;
    
    % t = 1;  % centroid
    % t = 0;  % midpoint of the segment
    t = 1;
    
    mesh.vertex(mesh.NV+ip).x = t*mesh.polygon(ip).xb+(1-t)*xm;
    mesh.vertex(mesh.NV+ip).y = t*mesh.polygon(ip).yb+(1-t)*ym;
    
end

mesh = rebuildmesh(mesh);

keg = [];

for ip=1:mesh.NP
    
    if mesh.polygon(ip).NV==3
        
        for ie=1:3
            for je=ie:3
                
                if not(ie==je)                    
            
                    ieg = mesh.polygon(ip).edges(ie);
                    jeg = mesh.polygon(ip).edges(je);
                    
                    pieg = sort([mesh.edge(ieg).pleft mesh.edge(ieg).pright]);
                    pjeg = sort([mesh.edge(jeg).pleft mesh.edge(jeg).pright]);
                    
                    if isequal(pieg,pjeg)
                        
                        exclude = [ieg, jeg];
                        break
                        
                    end
                end
            end
        end
           
        keg = [keg setdiff([mesh.polygon(ip).edges],exclude)];
        
    end

end

mesh = deledge2mesh(keg,mesh);

end

function mesh=deledge2mesh(xedges,oldmesh)

[mesh, ipx] = deledge(xedges,oldmesh);

mesh = rebuildmesh(mesh);

end

function [mesh,ipx]=deledge(xedges,oldmesh)

ipx = 0;

xedges = sort(unique(xedges),'descend');

mesh = oldmesh;

for xedge=xedges
    
    if xedge<1 || xedge>mesh.NE || oldmesh.edge(xedge).marker==1
        
        fprintf('%s: edge %d is on the boundary or out of range\n\n',mfilename,xedge)
        
    else
        
        iv1 = mesh.edge(xedge).v1;
        iv2 = mesh.edge(xedge).v2;
        
        % determino i 2 poligoni che hanno quello come edge
        
        ip1 = 0;
        ip2 = 0;
        
        for ip=1:mesh.NP
            if find(mesh.polygon(ip).edges==xedge)
                if ip1==0
                    ip1 = ip;
                else
                    ip2 = ip;
                end
            end
        end
        
        if ip1==0 || ip2==0
            if ip1>0
                ipx = ip1;
            elseif ip2>0
                ipx = ip2;
            else
                fprintf('%s: edge %d does not belong to any polygon\n\n',mfilename,xedge)
                return
            end
            fprintf('%s: edge %d belongs only to polygon %d\n\n',mfilename,xedge,ipx)
            return
        end
        
        % i poligoni sono ip1 e ip2
        
        % i vertici sono iv1 e iv2
        
        % parto dal vertice iv1 del poligono ip1
        
        Vip1 = [mesh.polygon(ip1).vertices mesh.polygon(ip1).vertices];
        Vip2 = [mesh.polygon(ip2).vertices mesh.polygon(ip2).vertices];
        
        Eip1 = [mesh.polygon(ip1).edges mesh.polygon(ip1).edges];
        Eip2 = [mesh.polygon(ip2).edges mesh.polygon(ip2).edges];
        
        % cerco [iv1 iv2] oppure [iv2 iv1]
        
        for i=1:length(Vip1)-1
            if isequal(Vip1(i:i+1),[iv1 iv2]) || isequal(Vip1(i:i+1),[iv2 iv1])
                % prendo quelli dopo
                Vok1 = Vip1(i+2:i+2+mesh.polygon(ip1).NV-2);
                break
            end
        end
        
        for i=1:length(Vip2)-1
            if isequal(Vip2(i:i+1),[iv1 iv2]) || isequal(Vip2(i:i+1),[iv2 iv1])
                % prendo quelli dopo
                Vok2 = Vip2(i+2:i+2+mesh.polygon(ip2).NV-2);
                break
            end
        end
        
        % cerco xedge
        
        for i=1:length(Eip1)-1
            if isequal(Eip1(i),xedge)
                % prendo quelli dopo
                Eok1 = Eip1(i+1:i+1-1+mesh.polygon(ip1).NE-1);
                break
            end
        end
        
        for i=1:length(Eip2)-1
            if isequal(Eip2(i),xedge)
                % prendo quelli dopo
                Eok2 = Eip2(i+1:i+1-1+mesh.polygon(ip2).NE-1);
                break
            end
        end
        
        % aggiorniamo il poligono ip1
        
        mesh.polygon(ip1).vertices = [Vok1 Vok2];
        mesh.polygon(ip1).NV = length([Vok1 Vok2]);
        mesh.polygon(ip1).edges = [Eok1 Eok2];
        mesh.polygon(ip1).NE = length([Eok1 Eok2]);

        % sposto i vertici avanti di uno
        
        vtmp = mesh.polygon(ip1).vertices;
        mesh.polygon(ip1).vertices = [vtmp(end) vtmp(1:end-1)];
        
        % eliminiamo il poligono ip2
        %
        mesh.polygon(ip2) = [];
        mesh.NP = mesh.NP-1;
        
        % eliminiamo l'edge xedge
        %
        mesh.edge(xedge) = [];
        mesh.NE = mesh.NE-1;
        %
        for ip=1:mesh.NP
            for ie=1:mesh.polygon(ip).NE
                if mesh.polygon(ip).edges(ie)>xedge
                    mesh.polygon(ip).edges(ie) = mesh.polygon(ip).edges(ie)-1;
                end
            end
        end
        
    end
    
end

end