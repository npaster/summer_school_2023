function mesh=conn2mesh(xv,yv,vertices)
%
% CONN2MESH(xv,yv,vertices)
%
% builds a VEM mesh starting from minimal structure
%
% (xv,yv) vertex coordinates
% vertices(ip,:) = vertices of polygon ip (padded with zeros)
%
nver = length(xv);
nele = size(vertices,1);
%
mesh.NV = nver;
%
% creiamo le strutture che conterranno le coord dei vertici e il marker
%
mesh.vertex(mesh.NV).x = [];
mesh.vertex(mesh.NV).y = [];
mesh.vertex(mesh.NV).marker = [];
%
mesh.NE = 0;
mesh.edge = [];
%
% read vertex coordinates
%
for ivg=1:mesh.NV
    %
    mesh.vertex(ivg).x = xv(ivg);
    mesh.vertex(ivg).y = yv(ivg);
    %
end
%
% numbers of polygons
%
mesh.NP = nele;
%
if mesh.NP==1
    %
    % in case there is only 1 polygon global=local!!!
    %
    % set vertices and edges
    %
    mesh.polygon(1).NV = nver;
    mesh.polygon(1).vertices = 1:nver;
    %
    mesh.polygon(1).NE = nver;
    mesh.polygon(1).edges = 1:nver;
    %
    mesh.NE = nver;
    %
    for iv=1:nver
        %
        mesh.vertex(iv).marker = 1;
        %
        ivp = iv+1;
        if ivp==nver+1
            ivp = 1;
        end
        %
        mesh.edge(iv).v1 = iv;
        mesh.edge(iv).v2 = ivp;
        %
        x1 = xv(iv);
        y1 = yv(iv);
        %
        x2 = xv(ivp);
        y2 = yv(ivp);
        %
        mesh.edge(iv).length = norm([x2-x1, y2-y1]);
        %
        mesh.edge(iv).type = 'segment';
        %
        mesh.edge(iv).pleft = 1;
        mesh.edge(iv).pright = -1;
        %
        mesh.edge(iv).marker = 1;
        %
    end
    %
else
    %
    % prepariamo la struttura dei polygons
    %
    mesh.polygon(mesh.NP).NV = 0;
    mesh.polygon(mesh.NP).vertices = [];
    %
    mesh.polygon(mesh.NP).NE = 0;
    mesh.polygon(mesh.NP).edges = [];
    %
    mesh.polygon(mesh.NP).xb = [];
    mesh.polygon(mesh.NP).yb = [];
    %
    mesh.polygon(mesh.NP).h = [];
    mesh.polygon(mesh.NP).area = [];
    %
    for ip=1:mesh.NP
        %
        vv = vertices(ip,:);
        %
        % prendo solo i non zeri
        %
        mesh.polygon(ip).vertices = vv(vv~=0);
        %
        mesh.polygon(ip).NV = length(mesh.polygon(ip).vertices);
        mesh.polygon(ip).NE = mesh.polygon(ip).NV;
    end
    %
    % generiamo la struttura degli edge
    %
    iedge = 0;
    edgematrix = zeros(1,2);
    %
    for ip=1:mesh.NP
        for iv=1:mesh.polygon(ip).NV
            ivp = iv+1;
            if ivp==mesh.polygon(ip).NV+1
                ivp = 1;
            end
            iedge = iedge+1;
            edgematrix(iedge,1) = mesh.polygon(ip).vertices(iv);
            edgematrix(iedge,2) = mesh.polygon(ip).vertices(ivp);
        end
    end
    %
    % ordiniamo gli estremi dal minore al maggiore
    %
    edgematrix = sort(edgematrix,2);
    %
    % eliminiamo gli uguali
    %
    edgematrix = unique(edgematrix,'rows');
    %
    % allochiamo la struttura edge
    %
    mesh.NE = size(edgematrix,1);
    mesh.edge(mesh.NE).v1 = [];
    mesh.edge(mesh.NE).v2 = [];
    %
    mesh.edge(mesh.NE).marker = [];
    mesh.edge(mesh.NE).length = [];
    %
    for ieg=1:mesh.NE
        %
        v1 = edgematrix(ieg,1);
        v2 = edgematrix(ieg,2);
        %
        mesh.edge(ieg).v1 = v1;
        mesh.edge(ieg).v2 = v2;
        %
        x1 = mesh.vertex(v1).x;
        y1 = mesh.vertex(v1).y;
        %
        x2 = mesh.vertex(v2).x;
        y2 = mesh.vertex(v2).y;
        %
        mesh.edge(ieg).length = norm([x2-x1, y2-y1]);
        %
        mesh.edge(ieg).type = 'segment';
        %
    end
    %
    % associamo ad ogni poligono i suoi edges
    %
    edges = [mesh.edge.v1]'+sqrt(-1)*[mesh.edge.v2]';
    %
    for ip=1:mesh.NP
        for ie=1:mesh.polygon(ip).NE
            iep = ie+1;
            if iep==mesh.polygon(ip).NV+1
                iep = 1;
            end
            v1 = mesh.polygon(ip).vertices(ie);
            v2 = mesh.polygon(ip).vertices(iep);
            %
            % mesh.polygon(ip).edges(ie) deve contenere l'edge e12 cui estremi sono v1 e v2
            %
            if v1<v2
                e12 = find(edges==v1+sqrt(-1)*v2);
                %e12 = find(ismembc(edges,v1+sqrt(-1)*v2));
            else
                e12 = find(edges==v2+sqrt(-1)*v1);
                %e12 = find(ismembc(edges,v2+sqrt(-1)*v1));
            end
            %
            if isempty(e12)
                fprintf('%s: edge not found\n\n',mfilename)
            else
                mesh.polygon(ip).edges(ie) = e12;
            end
        end
    end
    %
    % assegniamo a ogni edge i suoi poligoni
    %
    for ip=1:mesh.NP
        for ie=1:mesh.polygon(ip).NE
            iep = ie+1;
            if iep==mesh.polygon(ip).NV+1
                iep = 1;
            end
            v1 = mesh.polygon(ip).vertices(ie);
            v2 = mesh.polygon(ip).vertices(iep);
            %
            ieg = mesh.polygon(ip).edges(ie);
            %
            if mesh.edge(ieg).v1==v1 && mesh.edge(ieg).v2==v2
                % same orientation
                mesh.edge(ieg).pleft = ip;
            elseif mesh.edge(ieg).v1==v2 && mesh.edge(ieg).v2==v1
                % opposite orientation
                mesh.edge(ieg).pright = ip;
            else
                fprintf('%s: internal error!\n',mfilename);
            end
        end
    end
    
    for ieg=1:mesh.NE
        if isempty(mesh.edge(ieg).pright)
            mesh.edge(ieg).pright = 0;
        end
        if isempty(mesh.edge(ieg).pleft)
            mesh.edge(ieg).pleft = 0;
        end
    end
    %
    % marker dei vertici e degli edge
    %
    for ivg=1:mesh.NV
        mesh.vertex(ivg).marker = 0;
    end
    %
    for ieg=1:mesh.NE
        if mesh.edge(ieg).pright==0 || mesh.edge(ieg).pleft==0
            %
            % set edge marker = 1
            %
            mesh.edge(ieg).marker = 1;
            %
            % set vertex marker of the edge = 1
            %
            v1 = mesh.edge(ieg).v1;
            v2 = mesh.edge(ieg).v2;
            %
            mesh.vertex(v1).marker = 1;
            mesh.vertex(v2).marker = 1;
            %
        else
            %
            % set edge marker = 0
            %
            mesh.edge(ieg).marker = 0;
            %
        end
    end
end
%
% diametro, area e baricentro
%
for ip=1:mesh.NP
    
    % estraggo le coordinate dei vertici del poligono
    %
    x = [mesh.vertex([mesh.polygon(ip).vertices]).x];
    y = [mesh.vertex([mesh.polygon(ip).vertices]).y];
    
    % aggiungo in fondo il primo vertice
    %
    x = [x x(1)];
    y = [y y(1)];
    
    NV = mesh.polygon(ip).NV;
    
    % diametro (massimo delle distanze tra tutti i vertici)
    %
    h = 0;
    %
    for iv=1:NV
        for jv=1:NV
            d = norm([x(iv)-x(jv) y(iv)-y(jv)]);
            if d>h
                h = d;
            end
        end
    end
    %
    mesh.polygon(ip).h = h;
    
    % area of polygon
    %
    tmp = zeros(NV,1);
    %
    for iv=1:NV
        tmp(iv) = x(iv)*y(iv+1)-x(iv+1)*y(iv);
    end
    %
    area = 0.5*sum(tmp);
    %
    mesh.polygon(ip).area = area;
    
    % centroid of polygon
    %
    xb = 0;
    yb = 0;
    %
    for iv=1:NV
        xb = xb + (x(iv)+x(iv+1))*tmp(iv);
        yb = yb + (y(iv)+y(iv+1))*tmp(iv);
    end
    %
    mesh.polygon(ip).xb = xb/(6*area);
    mesh.polygon(ip).yb = yb/(6*area);
    
end
%
% mesh name
%
mesh.name = ['Generated by conn2mesh on ' datestr(now)];
%
return
