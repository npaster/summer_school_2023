function VEM=vembasis(p,ip,meshp,type,n1,n2,enhanced)
%
% VEM=VEMBASIS(K,IP,MESHP,TYPE,N1,N2,ENAHNCED)
%
% computes true VEM basis functions.
%
% K = VEM degree
% IP = polygon
% MESH = mesh
%
% TYPE can be 'vertex','edge','internal'
% switch TYPE
%    case {'vertex','internal'}
%       N1 = number of dof, N2 = dummy
%    case {'edge'}
%       N1 = edge number, N2 = dof on the edge
% end
%
% ENHANCED can be true or false
%

%
% other parameters:
%
%VEM.CHECK = true;
VEM.CHECK = false;
%
VEM.DRAW = true;
%VEM.DRAW = false;
%
%VEM.COMPUTE_ALL = true;
VEM.COMPUTE_ALL = false;
%
% polynomial degree for the FEM approximation
%
VEM.k = p;
%

%
%close all
%
VEM.enhanced = enhanced;
%
VEM.NV = meshp.polygon(ip).NV;
%
VEM.xv = [meshp.vertex([meshp.polygon(ip).vertices]).x];
VEM.yv = [meshp.vertex([meshp.polygon(ip).vertices]).y];
%
[xq,yq,wq] = quadrature_vianello(2*p,ip,meshp);
%
% storing quadrature nodes and weights in the mesh structure
%
meshp.polygon(ip).xq = xq;
meshp.polygon(ip).yq = yq;
meshp.polygon(ip).wq = wq;
%
q = vemprojectors(p,ip,meshp);
%
Id = eye(size(q.PNk));
%
% polynomial degree for the FEM approximation
%
k = VEM.k;
%
% definition of polygon and mesh
%
VEM.p = p;
%
% code starts
%
xv = VEM.xv;
yv = VEM.yv;
%
VEM.Ndof = p * VEM.NV + dimP(p-2);
%
% diametro (massimo delle distanze tra tutti i vertici)
%
VEM.h = meshp.polygon(ip).h;
%
VEM.area = meshp.polygon(ip).area;
%
% centroid of polygon
%
VEM.xb = meshp.polygon(ip).xb;
VEM.yb = meshp.polygon(ip).yb;
%
VEM.NE = VEM.NV;
%
%
%%% GAUSS-LOBATTO STUFF
%
% (tgl,wgl) = Gauss-Lobatto quadrature nodes and weigths on [-1 +1] 
%
[VEM.tgl,VEM.wgl] = mylglnodes(VEM.p);
%
% write poly file
%
omega = './triangle/tmp/polygon-automatic';
%
fid = fopen([omega '.poly'],'w');
%
fprintf(fid,'# numero vertici, 2, 0, marker imposto (1 si, 0 no)\n');
%
fprintf(fid,'%d 2 0 1\n',p*VEM.NV);
%
fprintf(fid,'# vertice, x, y, marker\n');
%
ic = 0;
for iv=1:VEM.NV
    ic = ic+1;
    fprintf(fid,'%d %e %e %d\n',ic,xv(iv),yv(iv),100*iv);
    %
    ivp = iv+1;
    if ivp==VEM.NV+1
        ivp = 1;
    end
    %
    x1 = xv(iv);
    y1 = yv(iv);
    x2 = xv(ivp);
    y2 = yv(ivp);
    for kp=1:p-1
        xkp = (x2-x1)/2*VEM.tgl(kp+1)+(x2+x1)/2;
        ykp = (y2-y1)/2*VEM.tgl(kp+1)+(y2+y1)/2;
        ic = ic+1;
        fprintf(fid,'%d %e %e %d\n',ic,xkp,ykp,iv);
    end
end
%
fprintf(fid,'# number of segments, marker imposto (1 si, 0 no)\n');
%
fprintf(fid,'%d 1\n',p*VEM.NV);
%
fprintf(fid,'# lato del bordo, estremo 1, estremo 2, marker\n');
%
ic = 0;
for iv=1:VEM.NV
    for kp=0:p-1
        ic = ic+1;
        icp = ic+1;
        if icp==p*VEM.NV+1
            icp = 1;
        end
        fprintf(fid,'%d %d %d %d\n',ic,ic,icp,iv);
    end
end
%
fprintf(fid,'# numero di buchi\n');
%
fprintf(fid,'0\n');
%
fprintf(fid,'# numero di vincoli di area\n');
%
fprintf(fid,'0\n');
%
%fprintf(fid,'# regione, punto nella regione, vincolo di area\n');
%
%fprintf(fid,'1 0.0 0.0 0.001\n');
%
fclose(fid);
%
% define mesh refinement
%
%sarea = '0.001';
sarea = sprintf('%f',meshp.polygon(ip).area/1000);
%
switch computer
    %
    %64-Bit Platforms
    %      PCWIN64  - Microsoft Windows on x64       1     0     0   win64
    %      GLNXA64  - Linux on x86_64                0     1     0   glnxa64
    %      MACI64   - Apple Mac OS X on x86_64       0     1     1   maci64
    %
    case 'PCWIN64'
        %
        % Microsoft Windows on x64
        %
        system(['.\triangle\PCWIN64\triangle.exe -pqnIa' sarea ' ' omega]);
        %
    case 'GLNXA64'
        %
        % Linux on x86_64
        %
        system(['./triangle/GLNXA64/triangle -pqnIa' sarea ' ' omega]);
        %
    case 'MACI64'
        %
        % Apple Mac OS X on x86_64
        %
        system(['./triangle/MACI64/triangle -pqnIa' sarea ' ' omega]);
        %
end
%
mesh = readmesh(omega);
%
%drawmesh(mesh)
%drawnow
%
if mesh.NV-mesh.NE+mesh.NP == 1
    fprintf('%s: Euler formula OK for mesh %s\n',mfilename,mesh.name)
end
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FEMk for triangles starts
%
% gradi di liberta' totali
%
Ndof = mesh.NV + (k-1)*mesh.NE + dimP(k-3)*mesh.NP;
%
% select quadrature to integrate exactly all polynomials
%
% 2k-2: grad grad term
% basic: (p-2)+k: right hand side (right hand side is a monomial up to degree p-2)
% enhanced: p+k: right hand side (right hand side is a monomial up to degree p)
%
if VEM.enhanced
    fdq = ['degree=' num2str(max([2*k-2,p+k]))];
else
    fdq = ['degree=' num2str(max([2*k-2,(p-2)+k]))];
end
%
%fdq = 'degree=13';
%
[Nq,xhq,yhq,whq] = quadratura(fdq);
%
% funzioni di base dell'elemento di riferimento
% calcolate nei nodi di quadratura
%
phihq = zeros(dimP(k),Nq);
%
for i=1:dimP(k)
    phihq(i,:) = phih(i,xhq,yhq,k);
end
%
% gradienti delle funzioni di base calcolati nei nodi di quadratura
%
gphihq = zeros(2,Nq,dimP(k));
%
% gphihq(<componente>,<nodo di quadratura>,<nodo di base>)
%
for i=1:dimP(k)
    [ghxq, ghyq] = gradhphih(i,xhq,yhq,k);
    %
    gphihq(1,:,i) = ghxq;
    gphihq(2,:,i) = ghyq;
end
%
if VEM.enhanced
    eNpm2 = dimP(p);
else
    eNpm2 = dimP(p-2);
end
%
% various matrices
%
A = sparse(Ndof,Ndof);
B = zeros(Ndof,eNpm2);
%
fprintf('%s: info: assembling stiffness matrix\n',mfilename)
tic
%
% compute matrix A
%
for ip=1:mesh.NP
    %
    if mod(ip,100)==0
        fprintf('%s: triangle %d out of %d\n',mesh.name,ip,mesh.NP);
    end
    %
    % get global dofs
    %
    gdof = getgdof(k,ip,mesh);
    %
    % coordinate dei nodi del triangolo ip
    %
    xv = [mesh.vertex([mesh.polygon(ip).vertices]).x];
    yv = [mesh.vertex([mesh.polygon(ip).vertices]).y];
    %
    % vettori dei lati del triangolo iele
    % e(i,:) = vettore del lato i
    %
    e(1,1) = xv(3)-xv(2); e(1,2) = yv(3)-yv(2);
    e(2,1) = xv(1)-xv(3); e(2,2) = yv(1)-yv(3);
    e(3,1) = xv(2)-xv(1); e(3,2) = yv(2)-yv(1);
    %
    % trasformazione F : T^ ---> T
    %
    % Jacobiana
    %
    J = [+e(3,1) -e(2,1)
        +e(3,2) -e(2,2)];
    %
    % nodi di quadratura sull'elemento T
    %
    xyq = [xv(1)*ones(1,Nq); yv(1)*ones(1,Nq)] + J * [xhq yhq]';
    %
    xq = xyq(1,:);
    yq = xyq(2,:);
    %
    % jacobiana dell'inversa trasposta
    %
    Jit = 1/(2*mesh.polygon(ip).area) * [
        -e(2,2) -e(3,2)
        +e(2,1) +e(3,1)
        ];
    %
    % diffusion
    %
    ATD = zeros(dimP(k),dimP(k));
    %
    for i=1:dimP(k)
        %
        B1 = Jit(1,1)*gphihq(1,:,i)+Jit(1,2)*gphihq(2,:,i);
        B2 = Jit(2,1)*gphihq(1,:,i)+Jit(2,2)*gphihq(2,:,i);
        %
        for j=1:i-1
            ATD(i,j) = ATD(j,i);
        end
        %
        for j=i:dimP(k)
            %
            A1 = Jit(1,1)*gphihq(1,:,j)+Jit(1,2)*gphihq(2,:,j);
            A2 = Jit(2,1)*gphihq(1,:,j)+Jit(2,2)*gphihq(2,:,j);
            %
            ATD(i,j) = 2*mesh.polygon(ip).area*sum((A1.*B1+A2.*B2).*whq');
            %
        end
    end
    %
    A(gdof,gdof) = A(gdof,gdof) + ATD;
    %
    % calcoliamo la matrice per colonne
    %
    %
    %BBT = zeros(length(gdof),length(eNpm2));
    %
    for igamma=1:eNpm2
        %
        gamma = o2m(igamma);
        %
        % monomio scalato valutato nei nodi di quadratura dell'elem. corrente
        %
        mq = VEM_m(gamma,xq,yq,VEM);
        %
        BT = zeros(dimP(k),1);
        %
        for i=1:dimP(k)
            %
            BT(i) = 2*mesh.polygon(ip).area * sum(phihq(i,:).*mq.*whq');
            %
        end
        %
        %BBT(:,igamma) = BT;
        %
        %B(gdof,igamma) = B(gdof,igamma) + BT;
        B(gdof,igamma) = BT;
        %
    end
    %B(gdof,1:eNpm2) = B(gdof,1:eNpm2) + BBT;
    %
end
%
Z = zeros(eNpm2,eNpm2);
%
% full mixed matrix without the bc
%
C = [A  B
     B' Z];
%
% qui inizia il codice per le phih VEM
%
VEMphi = zeros(Ndof,VEM.Ndof);
myVEMphi = zeros(Ndof,VEM.Ndof);
%
% coefficienti del polinomio
%
VEMc = zeros(eNpm2,VEM.Ndof);
myVEMc = zeros(eNpm2,VEM.Ndof);
%
for dof=1:VEM.Ndof
    %
    gdofL = [];
    %
    VEM = splitdof(VEM,dof);
    %
    uh = zeros(Ndof+eNpm2,1);
    b = zeros(Ndof+eNpm2,1);
    %
    if (strcmp(type,'vertex') && n1==VEM.iv) ||...
            (strcmp(type,'edge') && n1==VEM.ie && n2==VEM.ip) ||...
            (strcmp(type,'internal') && n1==VEM.ialpha) || VEM.COMPUTE_ALL    
        toc
        fprintf('%s: info: enforcing boundary conditions\n',mfilename)
        %
        %%% BOUNDARY CONDITIONS
        %
        % Dirichlet bc on vertices - update load term
        %
        for iv=1:mesh.NV
            ivmarker = mesh.vertex(iv).marker;
            if ivmarker > 0
                uh(iv) = g(VEM,ivmarker,mesh.vertex(iv).x,mesh.vertex(iv).y);
                b = b - uh(iv)*C(:,iv);
            else
                gdofL = [gdofL iv];
            end
        end
        %
        for ieg=1:mesh.NE
            %
            % edge marker
            %
            iegmarker = mesh.edge(ieg).marker;
            %
            % extrema of the edge globally ordered
            %
            x1 = mesh.vertex(mesh.edge(ieg).v1).x;
            y1 = mesh.vertex(mesh.edge(ieg).v1).y;
            %
            x2 = mesh.vertex(mesh.edge(ieg).v2).x;
            y2 = mesh.vertex(mesh.edge(ieg).v2).y;
            %
            for ik=1:k-1
                %
                % indice globale del nodo interno all'edge
                %
                igik = mesh.NV + (ieg-1)*(k-1)+ik;
                %
                if iegmarker > 0
                    %
                    % equispaced points
                    %
                    xik = (1-ik/k)*x1+(ik/k)*x2;
                    yik = (1-ik/k)*y1+(ik/k)*y2;
                    %
                    % Dirichlet boundary condition on the edge - update load term
                    %
                    uh(igik) = g(VEM,iegmarker,xik,yik);
                    b = b - uh(igik)*C(:,igik);
                    %
                else
                    %
                    % add ig to free nodes
                    %
                    gdofL = [gdofL igik];
                    %
                end
            end
        end
        %
        % adding all internal dofs to the free dofs
        %
        for ip=1:mesh.NP
            idofs = mesh.NV+mesh.NE*(k-1) + (((ip-1)*dimP(k-3)+1):((ip-1)*dimP(k-3)+dimP(k-3)));
            gdofL = [gdofL idofs];
        end
        %
        % adding polynomials dof
        %
        NL = [gdofL Ndof+(1:eNpm2)];
        %
        % update right-hand-side for polynomial dofs
        %
        if VEM.ialpha>0
            b(Ndof+VEM.ialpha) = b(Ndof+VEM.ialpha) + VEM.area;
        end
        %
        if VEM.enhanced
            %
            % questi ci sono sempre!!!!!
            %
            for igamma=(dimP(p-2)+1):dimP(p)
                b(Ndof+igamma) = b(Ndof+igamma) + q.C(igamma,dof);
            end
            %
        end
        %
        % extract submatrices corresponding to the free dofs
        %
        Ch = C(NL,NL);
        %
        bh = b(NL);
        %
        uh(NL) = Ch\bh;
        %
        VEMphi(:,dof) = uh(1:Ndof);
        VEMc(:,dof) = uh(Ndof+1:end);
        %
        % Schur stuff
        %
        % gdofL are the fem dof
        % Npm2L are the others
        %
        NdofL = (1:length(gdofL));
        Npm2L = length(gdofL)+(1:eNpm2);
        %
        myA = Ch(NdofL,NdofL);
        myB = Ch(NdofL,Npm2L);
        %
        myBt = Ch(Npm2L,NdofL);
        %
        mya = bh(NdofL);
        myb = bh(Npm2L);
        %
        myVEMc(:,dof) = (myB'*(myA\myB))\(myB'*(myA\mya)-myb);
        %
        myVEMphi(:,dof) = uh(1:Ndof);
        myVEMphi(gdofL,dof) = (myA\(mya-myB*myVEMc(:,dof)));
        %
    end
    %
    %
    if not(strcmp(type,'edge'))
        n2 = 0;
    end
    %
    if VEM.DRAW
        if (strcmp(type,'vertex') && n1==VEM.iv) ||...
            (strcmp(type,'edge') && n1==VEM.ie && n2==VEM.ip) ||...
            (strcmp(type,'internal') && n1==VEM.ialpha)
            %
            %ffigure(dof)
            clf
            if VEM.enhanced
                stitle = ['k=' num2str(p) ' - enhanced - '];
            else
                stitle = ['k=' num2str(p) ' - basic - '];
            end
            %
            title([stitle 'vertex:' num2str(VEM.iv) ' edge:' num2str(VEM.ie) ',' num2str(VEM.ip) ' internal:' num2str(VEM.ialpha)]);
            drawuh(VEMphi(:,dof),k,mesh)
            view(3)
            drawnow
            %
        end
    end
    %
end
%
% check if necessary
%
if VEM.CHECK
    %
    Mcheck = zeros(dimP(p-2),VEM.Ndof);
    %
    for igamma=1:dimP(p-2)
        for dof=1:VEM.Ndof
            gamma = o2m(igamma);
            %
            for ip=1:mesh.NP
                %
                if mod(ip,100)==0
                    fprintf('%s: triangle %d out of %d\n',mesh.name,ip,mesh.NP);
                end
                %
                % get global dofs of the element
                %
                gdof = getgdof(k,ip,mesh);
                %
                % coordinate dei nodi del triangolo ip
                %
                xv = [mesh.vertex([mesh.polygon(ip).vertices]).x];
                yv = [mesh.vertex([mesh.polygon(ip).vertices]).y];
                %
                % vettori dei lati del triangolo iele
                % e(i,:) = vettore del lato i
                %
                e(1,1) = xv(3)-xv(2); e(1,2) = yv(3)-yv(2);
                e(2,1) = xv(1)-xv(3); e(2,2) = yv(1)-yv(3);
                e(3,1) = xv(2)-xv(1); e(3,2) = yv(2)-yv(1);
                %
                % Area del triangolo
                %
                area = mesh.polygon(ip).area;
                %
                % trasformazione F : T^ ---> T
                %
                % trasformiamo i nodi di quadratura
                %
                % Jacobiana
                %
                J = [+e(3,1) -e(2,1)
                    +e(3,2) -e(2,2)];
                %
                % jacobiana dell'inversa trasposta
                %
                Jit = 1/(2*area) * [
                    -e(2,2) -e(3,2)
                    +e(2,1) +e(3,1)
                    ];
                %
                % punti di quadratura sull'elemento T
                %
                xyq = [xv(1)*ones(1,Nq); yv(1)*ones(1,Nq)] + J * [xhq yhq]';
                xq = xyq(1,:);
                yq = xyq(2,:);
                %
                % uh valutato nei nodi di quadratura dell'elemento di riferimento
                %
                VEMphiq = 0;
                %
                for i=1:dimP(k)
                    %
                    VEMphiq = VEMphiq + VEMphi(gdof(i),dof)*phihq(i,:);
                    %
                end
                %
                % monomio scalato valutato nei nodi di quadratura dell'elem. corrente
                %
                mq = VEM_m(gamma,xq,yq,VEM);
                %
                % momento
                %
                Mcheck(igamma,dof) = Mcheck(igamma,dof) + 2*area * sum(VEMphiq.*mq.*whq');
                %
            end
        end
    end
    %
    VEM.Mcheck = 1/VEM.area * Mcheck;
    %
end
%
VEM.phi= VEMphi;
%
% true stiffness matrix
%
VEM.KP = VEMphi'*A*VEMphi;
%
VEM.Kc = q.PNkstar'*q.Gt*q.PNkstar;
VEM.Ks = (Id-q.PNk)'*Id*(Id-q.PNk);
%
% vem stiffness matrix
%
VEM.KPh = q.PNkstar'*q.Gt*q.PNkstar + VEM.Ks;

end

function mesh=readmesh(omega)
    
% STRUTTURA DATI

% iv: generico vertice                iv=1:    nver
% iele: generico triangolo (elemento) iele=1:  nele
% iedge: generico edge (lato)         iedge=1: nedge

% (xv(iv),yv(iv)) = coordinate del vertice iv=1:nver
% vertexmarker(iv) = flag iv=1:nver

% vertices(iele,:) = vertici del triangolo iele [5 6 7], iele=1:nele

% edges(iele,:) = edges del triangolo iele [67* 73* 9*], iele=1:nele

% endpoints(iedge,:) = estremi dell'edge iedge [89 2], iedge=1:nedge
% edgemarker(iedge) = flag, iedge=1:nedge

% legge la mesh contenuta in omega

fid = fopen([omega '.node']);

tmp = fscanf(fid,'%d',4); % tmp = vettore
nver = tmp(1);

xv = zeros(nver,1); % vettore colonna di nver zeri
yv = zeros(nver,1);
vertexmarker = zeros(nver,1);

% leggiamo le coord e i flag dei nodi

for iv=1:nver
    %
    tmp = fscanf(fid,'%f',4);
    %
    xv(iv) = tmp(2);         % coord x
    yv(iv) = tmp(3);         % coord y
    %
    vertexmarker(iv) = tmp(4); % flag
    %
end

fclose(fid); % chiudo il file dei nodi

fid = fopen([omega '.ele']);
tmp = fscanf(fid,'%d',3);

nele = tmp(1);
vertices = zeros(nele,3);

for iele=1:nele
    tmp = fscanf(fid,'%d',4);
    vertices(iele,1) = tmp(2);
    vertices(iele,2) = tmp(3);
    vertices(iele,3) = tmp(4);
end

fclose(fid);

fid = fopen([omega '.neigh']);
tmp = fscanf(fid,'%d',2);

nele = tmp(1);
neigh = zeros(nele,3);

for iele=1:nele
    tmp = fscanf(fid,'%d',4);
    neigh(iele,1) = tmp(2);
    neigh(iele,2) = tmp(3);
    neigh(iele,3) = tmp(4);
end

fclose(fid);

% crea struttura edge

S = sparse(nele,nele);
T = sparse(nele,3);

nedge = 0;

endpoints = zeros(1,2);
edgemarker = zeros(2,1);

for iele=1:nele
    for i=1:3
        ielen = neigh(iele,i);
        if ielen==-1
            nedge = nedge + 1;
            T(iele,i) = nedge;
            endpoints(nedge,1:2) = setdiff(vertices(iele,:),vertices(iele,i));
            edgemarker(nedge) = 1;
        else
            if ~S(iele,ielen)
                nedge = nedge + 1;
                S(iele,ielen) = nedge;
                S(ielen,iele) = nedge;
                endpoints(nedge,1:2) = intersect(vertices(iele,:),vertices(ielen,:));
                edgemarker(nedge) = 0;
            end
        end
    end
end

edges1 = zeros(nele,3);
edges = zeros(nele,3);

for iele=1:nele
    fS = full(S(iele,:));
    fT = full(T(iele,:));
    edges1(iele,:) = [fS(find(fS)) fT(find(fT))];
end

for iele=1:nele
    v1 = vertices(iele,1);
    v2 = vertices(iele,2);
    v3 = vertices(iele,3);
    for i=1:3
        e = edges1(iele,i);
        v = endpoints(e,:);
        if length(intersect([v1 v2],v))==2
            edges(iele,3) = e;
        elseif length(intersect([v2 v3],v))==2
            edges(iele,1) = e;
        elseif length(intersect([v3 v1],v))==2
            edges(iele,2) = e;
        end
    end
end

mesh.NP = nele;

for iele=1:nele
    mesh.polygon(iele).NV = 3;
    mesh.polygon(iele).vertices = vertices(iele,:);
    mesh.polygon(iele).NE = 3;
    mesh.polygon(iele).edges = [edges(iele,3) edges(iele,1) edges(iele,2)];
end

mesh.NV = nver;

for iv=1:nver
    mesh.vertex(iv).x = xv(iv);
    mesh.vertex(iv).y = yv(iv);
    mesh.vertex(iv).marker = vertexmarker(iv);
end

mesh.NE = nedge;

for ie=1:nedge
    mesh.edge(ie).v1 = endpoints(ie,1);
    mesh.edge(ie).v2 = endpoints(ie,2);
    if edgemarker(ie)>0
        mesh.edge(ie).marker = min([mesh.vertex(endpoints(ie,:)).marker]);
    else
        mesh.edge(ie).marker = 0;
    end
    x1 = mesh.vertex(mesh.edge(ie).v1).x;
    y1 = mesh.vertex(mesh.edge(ie).v1).y;
    %
    x2 = mesh.vertex(mesh.edge(ie).v2).x;
    y2 = mesh.vertex(mesh.edge(ie).v2).y;
    mesh.edge(ie).length = norm([x2-x1,y2-y1]);
    %
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
            if d>h;
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
% nome della mesh
%

mesh.name = ['polygon - ' datestr(now)];

end

function z=phih(i,x,y,k)

switch k
    case 1
        switch i
            case 1
                z = 1 - y - x;
            case 2
                z = x;
            case 3
                z = y;
        end
    case 2
        switch i
            case 1
                z = 4.*x.*y - 3.*y - 3.*x + 2.*x.^2 + 2.*y.^2 + 1;
            case 2
                z = 2.*x.^2 - x;
            case 3
                z = 2.*y.^2 - y;
            case 4
                z = 4.*x.*y;
            case 5
                z = 4.*y - 4.*x.*y - 4.*y.^2;
            case 6
                z = 4.*x - 4.*x.*y - 4.*x.^2;
        end
    case 3
        switch i
            case 1
                z = 18.*x.*y - (11.*y)/2 - (11.*x)/2 - (27.*x.*y.^2)/2 - (27.*x.^2.*y)/2 + 9.*x.^2 - (9.*x.^3)/2 + 9.*y.^2 - (9.*y.^3)/2 + 1;
            case 2
                z = x - (9.*x.^2)/2 + (9.*x.^3)/2;
            case 3
                z = y - (9.*y.^2)/2 + (9.*y.^3)/2;
            case 4
                z = (27.*x.^2.*y)/2 - (9.*x.*y)/2;
            case 5
                z = (27.*x.*y.^2)/2 - (9.*x.*y)/2;
            case 6
                z = (9.*x.*y)/2 - (9.*y)/2 - (27.*x.*y.^2)/2 + 18.*y.^2 - (27.*y.^3)/2;
            case 7
                z = 9.*y - (45.*x.*y)/2 + 27.*x.*y.^2 + (27.*x.^2.*y)/2 - (45.*y.^2)/2 + (27.*y.^3)/2;
            case 8
                z = 9.*x - (45.*x.*y)/2 + (27.*x.*y.^2)/2 + 27.*x.^2.*y - (45.*x.^2)/2 + (27.*x.^3)/2;
            case 9
                z = (9.*x.*y)/2 - (9.*x)/2 - (27.*x.^2.*y)/2 + 18.*x.^2 - (27.*x.^3)/2;
            case 10
                z = 27.*x.*y - 27.*x.*y.^2 - 27.*x.^2.*y;
        end
    case 4
        switch i
            case 1
                z = 64.*x.^2.*y.^2 - (25.*y)/3 - (25.*x)/3 + (140.*x.*y)/3 - 80.*x.*y.^2 - 80.*x.^2.*y + (128.*x.*y.^3)/3 + (128.*x.^3.*y)/3 + (70.*x.^2)/3 - (80.*x.^3)/3 + (32.*x.^4)/3 + (70.*y.^2)/3 - (80.*y.^3)/3 + (32.*y.^4)/3 + 1;
            case 2
                z = (22.*x.^2)/3 - x - 16.*x.^3 + (32.*x.^4)/3;
            case 3
                z = (22.*y.^2)/3 - y - 16.*y.^3 + (32.*y.^4)/3;
            case 4
                z = (16.*x.*y)/3 - 32.*x.^2.*y + (128.*x.^3.*y)/3;
            case 5
                z = 64.*x.^2.*y.^2 + 4.*x.*y - 16.*x.*y.^2 - 16.*x.^2.*y;
            case 6
                z = (16.*x.*y)/3 - 32.*x.*y.^2 + (128.*x.*y.^3)/3;
            case 7
                z = (16.*y)/3 - (16.*x.*y)/3 + 32.*x.*y.^2 - (128.*x.*y.^3)/3 - (112.*y.^2)/3 + (224.*y.^3)/3 - (128.*y.^4)/3;
            case 8
                z = 64.*x.^2.*y.^2 - 12.*y + 28.*x.*y - 144.*x.*y.^2 - 16.*x.^2.*y + 128.*x.*y.^3 + 76.*y.^2 - 128.*y.^3 + 64.*y.^4;
            case 9
                z = 16.*y - 128.*x.^2.*y.^2 - (208.*x.*y)/3 + 192.*x.*y.^2 + 96.*x.^2.*y - 128.*x.*y.^3 - (128.*x.^3.*y)/3 - (208.*y.^2)/3 + 96.*y.^3 - (128.*y.^4)/3;
            case 10
                z = 16.*x - 128.*x.^2.*y.^2 - (208.*x.*y)/3 + 96.*x.*y.^2 + 192.*x.^2.*y - (128.*x.*y.^3)/3 - 128.*x.^3.*y - (208.*x.^2)/3 + 96.*x.^3 - (128.*x.^4)/3;
            case 11
                z = 64.*x.^2.*y.^2 - 12.*x + 28.*x.*y - 16.*x.*y.^2 - 144.*x.^2.*y + 128.*x.^3.*y + 76.*x.^2 - 128.*x.^3 + 64.*x.^4;
            case 12
                z = (16.*x)/3 - (16.*x.*y)/3 + 32.*x.^2.*y - (128.*x.^3.*y)/3 - (112.*x.^2)/3 + (224.*x.^3)/3 - (128.*x.^4)/3;
            case 13
                z = 256.*x.^2.*y.^2 + 96.*x.*y - 224.*x.*y.^2 - 224.*x.^2.*y + 128.*x.*y.^3 + 128.*x.^3.*y;
            case 14
                z = 32.*x.*y.^2 - 32.*x.*y - 128.*x.^2.*y.^2 + 160.*x.^2.*y - 128.*x.^3.*y;
            case 15
                z = 160.*x.*y.^2 - 32.*x.*y - 128.*x.^2.*y.^2 + 32.*x.^2.*y - 128.*x.*y.^3;
        end
    case 5
        switch i
            case 1
                z = (1875.*x.^2.*y.^2)/4 - (137.*y)/12 - (137.*x)/12 - (3125.*x.^2.*y.^3)/12 - (3125.*x.^3.*y.^2)/12 + (375.*x.*y)/4 - (2125.*x.*y.^2)/8 - (2125.*x.^2.*y)/8 + (625.*x.*y.^3)/2 + (625.*x.^3.*y)/2 - (3125.*x.*y.^4)/24 - (3125.*x.^4.*y)/24 + (375.*x.^2)/8 - (2125.*x.^3)/24 + (625.*x.^4)/8 - (625.*x.^5)/24 + (375.*y.^2)/8 - (2125.*y.^3)/24 + (625.*y.^4)/8 - (625.*y.^5)/24 + 1;
            case 2
                z = (x.*(875.*x.^2 - 250.*x - 1250.*x.^3 + 625.*x.^4 + 24))/24;
            case 3
                z = (y.*(875.*y.^2 - 250.*y - 1250.*y.^3 + 625.*y.^4 + 24))/24;
            case 4
                z = (25.*x.*y.*(55.*x - 150.*x.^2 + 125.*x.^3 - 6))/24;
            case 5
                z = (25.*x.*y.*(5.*y - 1).*(25.*x.^2 - 15.*x + 2))/12;
            case 6
                z = (25.*x.*y.*(5.*x - 1).*(25.*y.^2 - 15.*y + 2))/12;
            case 7
                z = (25.*x.*y.*(55.*y - 150.*y.^2 + 125.*y.^3 - 6))/24;
            case 8
                z = -(25.*y.*(x + y - 1).*(55.*y - 150.*y.^2 + 125.*y.^3 - 6))/24;
            case 9
                z = (25.*y.*(25.*y.^2 - 15.*y + 2).*(10.*x.*y - 9.*y - 9.*x + 5.*x.^2 + 5.*y.^2 + 4))/12;
            case 10
                z = -(25.*y.*(5.*y - 1).*(47.*x + 47.*y - 120.*x.*y + 75.*x.*y.^2 + 75.*x.^2.*y - 60.*x.^2 + 25.*x.^3 - 60.*y.^2 + 25.*y.^3 - 12))/12;
            case 11
                z = (25.*y.*(750.*x.^2.*y.^2 - 154.*y - 154.*x + 710.*x.*y - 1050.*x.*y.^2 - 1050.*x.^2.*y + 500.*x.*y.^3 + 500.*x.^3.*y + 355.*x.^2 - 350.*x.^3 + 125.*x.^4 + 355.*y.^2 - 350.*y.^3 + 125.*y.^4 + 24))/24;
            case 12
                z = (25.*x.*(750.*x.^2.*y.^2 - 154.*y - 154.*x + 710.*x.*y - 1050.*x.*y.^2 - 1050.*x.^2.*y + 500.*x.*y.^3 + 500.*x.^3.*y + 355.*x.^2 - 350.*x.^3 + 125.*x.^4 + 355.*y.^2 - 350.*y.^3 + 125.*y.^4 + 24))/24;
            case 13
                z = -(25.*x.*(5.*x - 1).*(47.*x + 47.*y - 120.*x.*y + 75.*x.*y.^2 + 75.*x.^2.*y - 60.*x.^2 + 25.*x.^3 - 60.*y.^2 + 25.*y.^3 - 12))/12;
            case 14
                z = (25.*x.*(25.*x.^2 - 15.*x + 2).*(10.*x.*y - 9.*y - 9.*x + 5.*x.^2 + 5.*y.^2 + 4))/12;
            case 15
                z = -(25.*x.*(x + y - 1).*(55.*x - 150.*x.^2 + 125.*x.^3 - 6))/24;
            case 16
                z = -(125.*x.*y.*(47.*x + 47.*y - 120.*x.*y + 75.*x.*y.^2 + 75.*x.^2.*y - 60.*x.^2 + 25.*x.^3 - 60.*y.^2 + 25.*y.^3 - 12))/6;
            case 17
                z = -(125.*x.*y.*(x + y - 1).*(25.*x.^2 - 15.*x + 2))/6;
            case 18
                z = -(125.*x.*y.*(x + y - 1).*(25.*y.^2 - 15.*y + 2))/6;
            case 19
                z = -(125.*x.*y.*(5.*x - 1).*(5.*y - 1).*(x + y - 1))/4;
            case 20
                z = (125.*x.*y.*(5.*y - 1).*(10.*x.*y - 9.*y - 9.*x + 5.*x.^2 + 5.*y.^2 + 4))/4;
            case 21
                z = (125.*x.*y.*(5.*x - 1).*(10.*x.*y - 9.*y - 9.*x + 5.*x.^2 + 5.*y.^2 + 4))/4;
        end
    case 6
        switch i
            case 1
                z = 1890.*x.^2.*y.^2 - (147.*y)/10 - (147.*x)/10 - 2268.*x.^2.*y.^3 - 2268.*x.^3.*y.^2 + 972.*x.^2.*y.^4 + 1296.*x.^3.*y.^3 + 972.*x.^4.*y.^2 + (812.*x.*y)/5 - (1323.*x.*y.^2)/2 - (1323.*x.^2.*y)/2 + 1260.*x.*y.^3 + 1260.*x.^3.*y - 1134.*x.*y.^4 - 1134.*x.^4.*y + (1944.*x.*y.^5)/5 + (1944.*x.^5.*y)/5 + (406.*x.^2)/5 - (441.*x.^3)/2 + 315.*x.^4 - (1134.*x.^5)/5 + (324.*x.^6)/5 + (406.*y.^2)/5 - (441.*y.^3)/2 + 315.*y.^4 - (1134.*y.^5)/5 + (324.*y.^6)/5 + 1;
            case 2
                z = (x.*(137.*x - 675.*x.^2 + 1530.*x.^3 - 1620.*x.^4 + 648.*x.^5 - 10))/10;
            case 3
                z = (y.*(137.*y - 675.*y.^2 + 1530.*y.^3 - 1620.*y.^4 + 648.*y.^5 - 10))/10;
            case 4
                z = (18.*x.*y.*(105.*x.^2 - 25.*x - 180.*x.^3 + 108.*x.^4 + 2))/5;
            case 5
                z = (9.*x.*y.*(6.*y - 1).*(11.*x - 36.*x.^2 + 36.*x.^3 - 1))/2;
            case 6
                z = 4.*x.*y.*(18.*x.^2 - 9.*x + 1).*(18.*y.^2 - 9.*y + 1);
            case 7
                z = (9.*x.*y.*(6.*x - 1).*(11.*y - 36.*y.^2 + 36.*y.^3 - 1))/2;
            case 8
                z = (18.*x.*y.*(105.*y.^2 - 25.*y - 180.*y.^3 + 108.*y.^4 + 2))/5;
            case 9
                z = -(18.*y.*(x + y - 1).*(105.*y.^2 - 25.*y - 180.*y.^3 + 108.*y.^4 + 2))/5;
            case 10
                z = (9.*y.*(11.*y - 36.*y.^2 + 36.*y.^3 - 1).*(12.*x.*y - 11.*y - 11.*x + 6.*x.^2 + 6.*y.^2 + 5))/2;
            case 11
                z = -4.*y.*(18.*y.^2 - 9.*y + 1).*(37.*x + 37.*y - 90.*x.*y + 54.*x.*y.^2 + 54.*x.^2.*y - 45.*x.^2 + 18.*x.^3 - 45.*y.^2 + 18.*y.^3 - 10);
            case 12
                z = (9.*y.*(6.*y - 1).*(216.*x.^2.*y.^2 - 57.*y - 57.*x + 238.*x.*y - 324.*x.*y.^2 - 324.*x.^2.*y + 144.*x.*y.^3 + 144.*x.^3.*y + 119.*x.^2 - 108.*x.^3 + 36.*x.^4 + 119.*y.^2 - 108.*y.^3 + 36.*y.^4 + 10))/2;
            case 13
                z = -(18.*y.*(87.*x + 87.*y - 2160.*x.^2.*y.^2 + 1080.*x.^2.*y.^3 + 1080.*x.^3.*y.^2 - 580.*x.*y + 1395.*x.*y.^2 + 1395.*x.^2.*y - 1440.*x.*y.^3 - 1440.*x.^3.*y + 540.*x.*y.^4 + 540.*x.^4.*y - 290.*x.^2 + 465.*x.^3 - 360.*x.^4 + 108.*x.^5 - 290.*y.^2 + 465.*y.^3 - 360.*y.^4 + 108.*y.^5 - 10))/5;
            case 14
                z = -(18.*x.*(87.*x + 87.*y - 2160.*x.^2.*y.^2 + 1080.*x.^2.*y.^3 + 1080.*x.^3.*y.^2 - 580.*x.*y + 1395.*x.*y.^2 + 1395.*x.^2.*y - 1440.*x.*y.^3 - 1440.*x.^3.*y + 540.*x.*y.^4 + 540.*x.^4.*y - 290.*x.^2 + 465.*x.^3 - 360.*x.^4 + 108.*x.^5 - 290.*y.^2 + 465.*y.^3 - 360.*y.^4 + 108.*y.^5 - 10))/5;
            case 15
                z = (9.*x.*(6.*x - 1).*(216.*x.^2.*y.^2 - 57.*y - 57.*x + 238.*x.*y - 324.*x.*y.^2 - 324.*x.^2.*y + 144.*x.*y.^3 + 144.*x.^3.*y + 119.*x.^2 - 108.*x.^3 + 36.*x.^4 + 119.*y.^2 - 108.*y.^3 + 36.*y.^4 + 10))/2;
            case 16
                z = -4.*x.*(18.*x.^2 - 9.*x + 1).*(37.*x + 37.*y - 90.*x.*y + 54.*x.*y.^2 + 54.*x.^2.*y - 45.*x.^2 + 18.*x.^3 - 45.*y.^2 + 18.*y.^3 - 10);
            case 17
                z = (9.*x.*(11.*x - 36.*x.^2 + 36.*x.^3 - 1).*(12.*x.*y - 11.*y - 11.*x + 6.*x.^2 + 6.*y.^2 + 5))/2;
            case 18
                z = -(18.*x.*(x + y - 1).*(105.*x.^2 - 25.*x - 180.*x.^3 + 108.*x.^4 + 2))/5;
            case 19
                z = 54.*x.*y.*(216.*x.^2.*y.^2 - 57.*y - 57.*x + 238.*x.*y - 324.*x.*y.^2 - 324.*x.^2.*y + 144.*x.*y.^3 + 144.*x.^3.*y + 119.*x.^2 - 108.*x.^3 + 36.*x.^4 + 119.*y.^2 - 108.*y.^3 + 36.*y.^4 + 10);
            case 20
                z = -54.*x.*y.*(x + y - 1).*(11.*x - 36.*x.^2 + 36.*x.^3 - 1);
            case 21
                z = -54.*x.*y.*(x + y - 1).*(11.*y - 36.*y.^2 + 36.*y.^3 - 1);
            case 22
                z = -36.*x.*y.*(6.*y - 1).*(x + y - 1).*(18.*x.^2 - 9.*x + 1);
            case 23
                z = -36.*x.*y.*(6.*x - 1).*(x + y - 1).*(18.*y.^2 - 9.*y + 1);
            case 24
                z = 36.*x.*y.*(18.*y.^2 - 9.*y + 1).*(12.*x.*y - 11.*y - 11.*x + 6.*x.^2 + 6.*y.^2 + 5);
            case 25
                z = -36.*x.*y.*(6.*y - 1).*(37.*x + 37.*y - 90.*x.*y + 54.*x.*y.^2 + 54.*x.^2.*y - 45.*x.^2 + 18.*x.^3 - 45.*y.^2 + 18.*y.^3 - 10);
            case 26
                z = -36.*x.*y.*(6.*x - 1).*(37.*x + 37.*y - 90.*x.*y + 54.*x.*y.^2 + 54.*x.^2.*y - 45.*x.^2 + 18.*x.^3 - 45.*y.^2 + 18.*y.^3 - 10);
            case 27
                z = 36.*x.*y.*(18.*x.^2 - 9.*x + 1).*(12.*x.*y - 11.*y - 11.*x + 6.*x.^2 + 6.*y.^2 + 5);
            case 28
                z = 27.*x.*y.*(6.*x - 1).*(6.*y - 1).*(12.*x.*y - 11.*y - 11.*x + 6.*x.^2 + 6.*y.^2 + 5);
        end
    case 7
        switch i
            case 1
                z = (16807.*x.^2.*y.^2)/3 - (363.*y)/20 - (363.*x)/20 - (386561.*x.^2.*y.^3)/36 - (386561.*x.^3.*y.^2)/36 + (117649.*x.^2.*y.^4)/12 + (117649.*x.^3.*y.^3)/9 + (117649.*x.^4.*y.^2)/12 - (823543.*x.^2.*y.^5)/240 - (823543.*x.^3.*y.^4)/144 - (823543.*x.^4.*y.^3)/144 - (823543.*x.^5.*y.^2)/240 + (22981.*x.*y)/90 - (331681.*x.*y.^2)/240 - (331681.*x.^2.*y)/240 + (33614.*x.*y.^3)/9 + (33614.*x.^3.*y)/9 - (386561.*x.*y.^4)/72 - (386561.*x.^4.*y)/72 + (117649.*x.*y.^5)/30 + (117649.*x.^5.*y)/30 - (823543.*x.*y.^6)/720 - (823543.*x.^6.*y)/720 + (22981.*x.^2)/180 - (331681.*x.^3)/720 + (16807.*x.^4)/18 - (386561.*x.^5)/360 + (117649.*x.^6)/180 - (117649.*x.^7)/720 + (22981.*y.^2)/180 - (331681.*y.^3)/720 + (16807.*y.^4)/18 - (386561.*y.^5)/360 + (117649.*y.^6)/180 - (117649.*y.^7)/720 + 1;
            case 2
                z = (x.*(79576.*x.^2 - 12348.*x - 252105.*x.^3 + 420175.*x.^4 - 352947.*x.^5 + 117649.*x.^6 + 720))/720;
            case 3
                z = (y.*(79576.*y.^2 - 12348.*y - 252105.*y.^3 + 420175.*y.^4 - 352947.*y.^5 + 117649.*y.^6 + 720))/720;
            case 4
                z = (49.*x.*y.*(1918.*x - 11025.*x.^2 + 29155.*x.^3 - 36015.*x.^4 + 16807.*x.^5 - 120))/720;
            case 5
                z = (49.*x.*y.*(7.*y - 1).*(1715.*x.^2 - 350.*x - 3430.*x.^3 + 2401.*x.^4 + 24))/240;
            case 6
                z = (49.*x.*y.*(49.*y.^2 - 21.*y + 2).*(77.*x - 294.*x.^2 + 343.*x.^3 - 6))/144;
            case 7
                z = (49.*x.*y.*(49.*x.^2 - 21.*x + 2).*(77.*y - 294.*y.^2 + 343.*y.^3 - 6))/144;
            case 8
                z = (49.*x.*y.*(7.*x - 1).*(1715.*y.^2 - 350.*y - 3430.*y.^3 + 2401.*y.^4 + 24))/240;
            case 9
                z = (49.*x.*y.*(1918.*y - 11025.*y.^2 + 29155.*y.^3 - 36015.*y.^4 + 16807.*y.^5 - 120))/720;
            case 10
                z = -(49.*y.*(x + y - 1).*(1918.*y - 11025.*y.^2 + 29155.*y.^3 - 36015.*y.^4 + 16807.*y.^5 - 120))/720;
            case 11
                z = (49.*y.*(1715.*y.^2 - 350.*y - 3430.*y.^3 + 2401.*y.^4 + 24).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/240;
            case 12
                z = -(49.*y.*(77.*y - 294.*y.^2 + 343.*y.^3 - 6).*(107.*x + 107.*y - 252.*x.*y + 147.*x.*y.^2 + 147.*x.^2.*y - 126.*x.^2 + 49.*x.^3 - 126.*y.^2 + 49.*y.^3 - 30))/144;
            case 13
                z = (49.*y.*(49.*y.^2 - 21.*y + 2).*(2058.*x.^2.*y.^2 - 638.*y - 638.*x + 2506.*x.*y - 3234.*x.*y.^2 - 3234.*x.^2.*y + 1372.*x.*y.^3 + 1372.*x.^3.*y + 1253.*x.^2 - 1078.*x.^3 + 343.*x.^4 + 1253.*y.^2 - 1078.*y.^3 + 343.*y.^4 + 120))/144;
            case 14
                z = -(49.*y.*(7.*y - 1).*(2754.*x + 2754.*y - 51450.*x.^2.*y.^2 + 24010.*x.^2.*y.^3 + 24010.*x.^3.*y.^2 - 16450.*x.*y + 36015.*x.*y.^2 + 36015.*x.^2.*y - 34300.*x.*y.^3 - 34300.*x.^3.*y + 12005.*x.*y.^4 + 12005.*x.^4.*y - 8225.*x.^2 + 12005.*x.^3 - 8575.*x.^4 + 2401.*x.^5 - 8225.*y.^2 + 12005.*y.^3 - 8575.*y.^4 + 2401.*y.^5 - 360))/240;
            case 15
                z = (49.*y.*(607110.*x.^2.*y.^2 - 8028.*y - 8028.*x - 648270.*x.^2.*y.^3 - 648270.*x.^3.*y.^2 + 252105.*x.^2.*y.^4 + 336140.*x.^3.*y.^3 + 252105.*x.^4.*y.^2 + 71456.*x.*y - 244755.*x.*y.^2 - 244755.*x.^2.*y + 404740.*x.*y.^3 + 404740.*x.^3.*y - 324135.*x.*y.^4 - 324135.*x.^4.*y + 100842.*x.*y.^5 + 100842.*x.^5.*y + 35728.*x.^2 - 81585.*x.^3 + 101185.*x.^4 - 64827.*x.^5 + 16807.*x.^6 + 35728.*y.^2 - 81585.*y.^3 + 101185.*y.^4 - 64827.*y.^5 + 16807.*y.^6 + 720))/720;
            case 16
                z = (49.*x.*(607110.*x.^2.*y.^2 - 8028.*y - 8028.*x - 648270.*x.^2.*y.^3 - 648270.*x.^3.*y.^2 + 252105.*x.^2.*y.^4 + 336140.*x.^3.*y.^3 + 252105.*x.^4.*y.^2 + 71456.*x.*y - 244755.*x.*y.^2 - 244755.*x.^2.*y + 404740.*x.*y.^3 + 404740.*x.^3.*y - 324135.*x.*y.^4 - 324135.*x.^4.*y + 100842.*x.*y.^5 + 100842.*x.^5.*y + 35728.*x.^2 - 81585.*x.^3 + 101185.*x.^4 - 64827.*x.^5 + 16807.*x.^6 + 35728.*y.^2 - 81585.*y.^3 + 101185.*y.^4 - 64827.*y.^5 + 16807.*y.^6 + 720))/720;
            case 17
                z = -(49.*x.*(7.*x - 1).*(2754.*x + 2754.*y - 51450.*x.^2.*y.^2 + 24010.*x.^2.*y.^3 + 24010.*x.^3.*y.^2 - 16450.*x.*y + 36015.*x.*y.^2 + 36015.*x.^2.*y - 34300.*x.*y.^3 - 34300.*x.^3.*y + 12005.*x.*y.^4 + 12005.*x.^4.*y - 8225.*x.^2 + 12005.*x.^3 - 8575.*x.^4 + 2401.*x.^5 - 8225.*y.^2 + 12005.*y.^3 - 8575.*y.^4 + 2401.*y.^5 - 360))/240;
            case 18
                z = (49.*x.*(49.*x.^2 - 21.*x + 2).*(2058.*x.^2.*y.^2 - 638.*y - 638.*x + 2506.*x.*y - 3234.*x.*y.^2 - 3234.*x.^2.*y + 1372.*x.*y.^3 + 1372.*x.^3.*y + 1253.*x.^2 - 1078.*x.^3 + 343.*x.^4 + 1253.*y.^2 - 1078.*y.^3 + 343.*y.^4 + 120))/144;
            case 19
                z = -(49.*x.*(77.*x - 294.*x.^2 + 343.*x.^3 - 6).*(107.*x + 107.*y - 252.*x.*y + 147.*x.*y.^2 + 147.*x.^2.*y - 126.*x.^2 + 49.*x.^3 - 126.*y.^2 + 49.*y.^3 - 30))/144;
            case 20
                z = (49.*x.*(1715.*x.^2 - 350.*x - 3430.*x.^3 + 2401.*x.^4 + 24).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/240;
            case 21
                z = -(49.*x.*(x + y - 1).*(1918.*x - 11025.*x.^2 + 29155.*x.^3 - 36015.*x.^4 + 16807.*x.^5 - 120))/720;
            case 22
                z = -(343.*x.*y.*(2754.*x + 2754.*y - 51450.*x.^2.*y.^2 + 24010.*x.^2.*y.^3 + 24010.*x.^3.*y.^2 - 16450.*x.*y + 36015.*x.*y.^2 + 36015.*x.^2.*y - 34300.*x.*y.^3 - 34300.*x.^3.*y + 12005.*x.*y.^4 + 12005.*x.^4.*y - 8225.*x.^2 + 12005.*x.^3 - 8575.*x.^4 + 2401.*x.^5 - 8225.*y.^2 + 12005.*y.^3 - 8575.*y.^4 + 2401.*y.^5 - 360))/120;
            case 23
                z = -(343.*x.*y.*(x + y - 1).*(1715.*x.^2 - 350.*x - 3430.*x.^3 + 2401.*x.^4 + 24))/120;
            case 24
                z = -(343.*x.*y.*(x + y - 1).*(1715.*y.^2 - 350.*y - 3430.*y.^3 + 2401.*y.^4 + 24))/120;
            case 25
                z = -(343.*x.*y.*(7.*y - 1).*(x + y - 1).*(77.*x - 294.*x.^2 + 343.*x.^3 - 6))/48;
            case 26
                z = -(343.*x.*y.*(x + y - 1).*(49.*x.^2 - 21.*x + 2).*(49.*y.^2 - 21.*y + 2))/36;
            case 27
                z = -(343.*x.*y.*(7.*x - 1).*(x + y - 1).*(77.*y - 294.*y.^2 + 343.*y.^3 - 6))/48;
            case 28
                z = (343.*x.*y.*(77.*y - 294.*y.^2 + 343.*y.^3 - 6).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/48;
            case 29
                z = -(343.*x.*y.*(49.*y.^2 - 21.*y + 2).*(107.*x + 107.*y - 252.*x.*y + 147.*x.*y.^2 + 147.*x.^2.*y - 126.*x.^2 + 49.*x.^3 - 126.*y.^2 + 49.*y.^3 - 30))/36;
            case 30
                z = (343.*x.*y.*(7.*y - 1).*(2058.*x.^2.*y.^2 - 638.*y - 638.*x + 2506.*x.*y - 3234.*x.*y.^2 - 3234.*x.^2.*y + 1372.*x.*y.^3 + 1372.*x.^3.*y + 1253.*x.^2 - 1078.*x.^3 + 343.*x.^4 + 1253.*y.^2 - 1078.*y.^3 + 343.*y.^4 + 120))/48;
            case 31
                z = (343.*x.*y.*(7.*x - 1).*(2058.*x.^2.*y.^2 - 638.*y - 638.*x + 2506.*x.*y - 3234.*x.*y.^2 - 3234.*x.^2.*y + 1372.*x.*y.^3 + 1372.*x.^3.*y + 1253.*x.^2 - 1078.*x.^3 + 343.*x.^4 + 1253.*y.^2 - 1078.*y.^3 + 343.*y.^4 + 120))/48;
            case 32
                z = -(343.*x.*y.*(49.*x.^2 - 21.*x + 2).*(107.*x + 107.*y - 252.*x.*y + 147.*x.*y.^2 + 147.*x.^2.*y - 126.*x.^2 + 49.*x.^3 - 126.*y.^2 + 49.*y.^3 - 30))/36;
            case 33
                z = (343.*x.*y.*(77.*x - 294.*x.^2 + 343.*x.^3 - 6).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/48;
            case 34
                z = -(343.*x.*y.*(7.*x - 1).*(7.*y - 1).*(107.*x + 107.*y - 252.*x.*y + 147.*x.*y.^2 + 147.*x.^2.*y - 126.*x.^2 + 49.*x.^3 - 126.*y.^2 + 49.*y.^3 - 30))/24;
            case 35
                z = (343.*x.*y.*(7.*y - 1).*(49.*x.^2 - 21.*x + 2).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/24;
            case 36
                z = (343.*x.*y.*(7.*x - 1).*(49.*y.^2 - 21.*y + 2).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/24;
        end
    case 8
        switch i
            case 1
                z = (68416.*x.^2.*y.^2)/5 - (761.*y)/35 - (761.*x)/35 - 36864.*x.^2.*y.^3 - 36864.*x.^3.*y.^2 + 53248.*x.^2.*y.^4 + (212992.*x.^3.*y.^3)/3 + 53248.*x.^4.*y.^2 - (196608.*x.^2.*y.^5)/5 - 65536.*x.^3.*y.^4 - 65536.*x.^4.*y.^3 - (196608.*x.^5.*y.^2)/5 + (524288.*x.^2.*y.^6)/45 + (1048576.*x.^3.*y.^5)/45 + (262144.*x.^4.*y.^4)/9 + (1048576.*x.^5.*y.^3)/45 + (524288.*x.^6.*y.^2)/45 + (118124.*x.*y)/315 - (12816.*x.*y.^2)/5 - (12816.*x.^2.*y)/5 + (136832.*x.*y.^3)/15 + (136832.*x.^3.*y)/15 - 18432.*x.*y.^4 - 18432.*x.^4.*y + (106496.*x.*y.^5)/5 + (106496.*x.^5.*y)/5 - (65536.*x.*y.^6)/5 - (65536.*x.^6.*y)/5 + (1048576.*x.*y.^7)/315 + (1048576.*x.^7.*y)/315 + (59062.*x.^2)/315 - (4272.*x.^3)/5 + (34208.*x.^4)/15 - (18432.*x.^5)/5 + (53248.*x.^6)/15 - (65536.*x.^7)/35 + (131072.*x.^8)/315 + (59062.*y.^2)/315 - (4272.*y.^3)/5 + (34208.*y.^4)/15 - (18432.*y.^5)/5 + (53248.*y.^6)/15 - (65536.*y.^7)/35 + (131072.*y.^8)/315 + 1;
            case 2
                z = (x.*(6534.*x - 52528.*x.^2 + 216608.*x.^3 - 501760.*x.^4 + 659456.*x.^5 - 458752.*x.^6 + 131072.*x.^7 - 315))/315;
            case 3
                z = (y.*(6534.*y - 52528.*y.^2 + 216608.*y.^3 - 501760.*y.^4 + 659456.*y.^5 - 458752.*y.^6 + 131072.*y.^7 - 315))/315;
            case 4
                z = (64.*x.*y.*(6496.*x.^2 - 882.*x - 23520.*x.^3 + 44800.*x.^4 - 43008.*x.^5 + 16384.*x.^6 + 45))/315;
            case 5
                z = (16.*x.*y.*(8.*y - 1).*(274.*x - 1800.*x.^2 + 5440.*x.^3 - 7680.*x.^4 + 4096.*x.^5 - 15))/45;
            case 6
                z = (64.*x.*y.*(32.*y.^2 - 12.*y + 1).*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3))/45;
            case 7
                z = (4.*x.*y.*(44.*x - 192.*x.^2 + 256.*x.^3 - 3).*(44.*y - 192.*y.^2 + 256.*y.^3 - 3))/9;
            case 8
                z = (64.*x.*y.*(32.*x.^2 - 12.*x + 1).*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3))/45;
            case 9
                z = (16.*x.*y.*(8.*x - 1).*(274.*y - 1800.*y.^2 + 5440.*y.^3 - 7680.*y.^4 + 4096.*y.^5 - 15))/45;
            case 10
                z = (64.*x.*y.*(6496.*y.^2 - 882.*y - 23520.*y.^3 + 44800.*y.^4 - 43008.*y.^5 + 16384.*y.^6 + 45))/315;
            case 11
                z = -(64.*y.*(x + y - 1).*(6496.*y.^2 - 882.*y - 23520.*y.^3 + 44800.*y.^4 - 43008.*y.^5 + 16384.*y.^6 + 45))/315;
            case 12
                z = (16.*y.*(274.*y - 1800.*y.^2 + 5440.*y.^3 - 7680.*y.^4 + 4096.*y.^5 - 15).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/45;
            case 13
                z = -(64.*y.*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3).*(73.*x + 73.*y - 168.*x.*y + 96.*x.*y.^2 + 96.*x.^2.*y - 84.*x.^2 + 32.*x.^3 - 84.*y.^2 + 32.*y.^3 - 21))/45;
            case 14
                z = (4.*y.*(44.*y - 192.*y.^2 + 256.*y.^3 - 3).*(1536.*x.^2.*y.^2 - 533.*y - 533.*x + 2008.*x.*y - 2496.*x.*y.^2 - 2496.*x.^2.*y + 1024.*x.*y.^3 + 1024.*x.^3.*y + 1004.*x.^2 - 832.*x.^3 + 256.*x.^4 + 1004.*y.^2 - 832.*y.^3 + 256.*y.^4 + 105))/9;
            case 15
                z = -(64.*y.*(32.*y.^2 - 12.*y + 1).*(743.*x + 743.*y - 11520.*x.^2.*y.^2 + 5120.*x.^2.*y.^3 + 5120.*x.^3.*y.^2 - 4140.*x.*y + 8520.*x.*y.^2 + 8520.*x.^2.*y - 7680.*x.*y.^3 - 7680.*x.^3.*y + 2560.*x.*y.^4 + 2560.*x.^4.*y - 2070.*x.^2 + 2840.*x.^3 - 1920.*x.^4 + 512.*x.^5 - 2070.*y.^2 + 2840.*y.^3 - 1920.*y.^4 + 512.*y.^5 - 105))/45;
            case 16
                z = (16.*y.*(8.*y - 1).*(170880.*x.^2.*y.^2 - 3069.*y - 3069.*x - 168960.*x.^2.*y.^3 - 168960.*x.^3.*y.^2 + 61440.*x.^2.*y.^4 + 81920.*x.^3.*y.^3 + 61440.*x.^4.*y.^2 + 24308.*x.*y - 75240.*x.*y.^2 - 75240.*x.^2.*y + 113920.*x.*y.^3 + 113920.*x.^3.*y - 84480.*x.*y.^4 - 84480.*x.^4.*y + 24576.*x.*y.^5 + 24576.*x.^5.*y + 12154.*x.^2 - 25080.*x.^3 + 28480.*x.^4 - 16896.*x.^5 + 4096.*x.^6 + 12154.*y.^2 - 25080.*y.^3 + 28480.*y.^4 - 16896.*y.^5 + 4096.*y.^6 + 315))/45;
            case 17
                z = -(64.*y.*(4329.*x + 4329.*y - 772800.*x.^2.*y.^2 + 1308160.*x.^2.*y.^3 + 1308160.*x.^3.*y.^2 - 1075200.*x.^2.*y.^4 - 1433600.*x.^3.*y.^3 - 1075200.*x.^4.*y.^2 + 344064.*x.^2.*y.^5 + 573440.*x.^3.*y.^4 + 573440.*x.^4.*y.^3 + 344064.*x.^5.*y.^2 - 48860.*x.*y + 221088.*x.*y.^2 + 221088.*x.^2.*y - 515200.*x.*y.^3 - 515200.*x.^3.*y + 654080.*x.*y.^4 + 654080.*x.^4.*y - 430080.*x.*y.^5 - 430080.*x.^5.*y + 114688.*x.*y.^6 + 114688.*x.^6.*y - 24430.*x.^2 + 73696.*x.^3 - 128800.*x.^4 + 130816.*x.^5 - 71680.*x.^6 + 16384.*x.^7 - 24430.*y.^2 + 73696.*y.^3 - 128800.*y.^4 + 130816.*y.^5 - 71680.*y.^6 + 16384.*y.^7 - 315))/315;
            case 18
                z = -(64.*x.*(4329.*x + 4329.*y - 772800.*x.^2.*y.^2 + 1308160.*x.^2.*y.^3 + 1308160.*x.^3.*y.^2 - 1075200.*x.^2.*y.^4 - 1433600.*x.^3.*y.^3 - 1075200.*x.^4.*y.^2 + 344064.*x.^2.*y.^5 + 573440.*x.^3.*y.^4 + 573440.*x.^4.*y.^3 + 344064.*x.^5.*y.^2 - 48860.*x.*y + 221088.*x.*y.^2 + 221088.*x.^2.*y - 515200.*x.*y.^3 - 515200.*x.^3.*y + 654080.*x.*y.^4 + 654080.*x.^4.*y - 430080.*x.*y.^5 - 430080.*x.^5.*y + 114688.*x.*y.^6 + 114688.*x.^6.*y - 24430.*x.^2 + 73696.*x.^3 - 128800.*x.^4 + 130816.*x.^5 - 71680.*x.^6 + 16384.*x.^7 - 24430.*y.^2 + 73696.*y.^3 - 128800.*y.^4 + 130816.*y.^5 - 71680.*y.^6 + 16384.*y.^7 - 315))/315;
            case 19
                z = (16.*x.*(8.*x - 1).*(170880.*x.^2.*y.^2 - 3069.*y - 3069.*x - 168960.*x.^2.*y.^3 - 168960.*x.^3.*y.^2 + 61440.*x.^2.*y.^4 + 81920.*x.^3.*y.^3 + 61440.*x.^4.*y.^2 + 24308.*x.*y - 75240.*x.*y.^2 - 75240.*x.^2.*y + 113920.*x.*y.^3 + 113920.*x.^3.*y - 84480.*x.*y.^4 - 84480.*x.^4.*y + 24576.*x.*y.^5 + 24576.*x.^5.*y + 12154.*x.^2 - 25080.*x.^3 + 28480.*x.^4 - 16896.*x.^5 + 4096.*x.^6 + 12154.*y.^2 - 25080.*y.^3 + 28480.*y.^4 - 16896.*y.^5 + 4096.*y.^6 + 315))/45;
            case 20
                z = -(64.*x.*(32.*x.^2 - 12.*x + 1).*(743.*x + 743.*y - 11520.*x.^2.*y.^2 + 5120.*x.^2.*y.^3 + 5120.*x.^3.*y.^2 - 4140.*x.*y + 8520.*x.*y.^2 + 8520.*x.^2.*y - 7680.*x.*y.^3 - 7680.*x.^3.*y + 2560.*x.*y.^4 + 2560.*x.^4.*y - 2070.*x.^2 + 2840.*x.^3 - 1920.*x.^4 + 512.*x.^5 - 2070.*y.^2 + 2840.*y.^3 - 1920.*y.^4 + 512.*y.^5 - 105))/45;
            case 21
                z = (4.*x.*(44.*x - 192.*x.^2 + 256.*x.^3 - 3).*(1536.*x.^2.*y.^2 - 533.*y - 533.*x + 2008.*x.*y - 2496.*x.*y.^2 - 2496.*x.^2.*y + 1024.*x.*y.^3 + 1024.*x.^3.*y + 1004.*x.^2 - 832.*x.^3 + 256.*x.^4 + 1004.*y.^2 - 832.*y.^3 + 256.*y.^4 + 105))/9;
            case 22
                z = -(64.*x.*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3).*(73.*x + 73.*y - 168.*x.*y + 96.*x.*y.^2 + 96.*x.^2.*y - 84.*x.^2 + 32.*x.^3 - 84.*y.^2 + 32.*y.^3 - 21))/45;
            case 23
                z = (16.*x.*(274.*x - 1800.*x.^2 + 5440.*x.^3 - 7680.*x.^4 + 4096.*x.^5 - 15).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/45;
            case 24
                z = -(64.*x.*(x + y - 1).*(6496.*x.^2 - 882.*x - 23520.*x.^3 + 44800.*x.^4 - 43008.*x.^5 + 16384.*x.^6 + 45))/315;
            case 25
                z = (256.*x.*y.*(170880.*x.^2.*y.^2 - 3069.*y - 3069.*x - 168960.*x.^2.*y.^3 - 168960.*x.^3.*y.^2 + 61440.*x.^2.*y.^4 + 81920.*x.^3.*y.^3 + 61440.*x.^4.*y.^2 + 24308.*x.*y - 75240.*x.*y.^2 - 75240.*x.^2.*y + 113920.*x.*y.^3 + 113920.*x.^3.*y - 84480.*x.*y.^4 - 84480.*x.^4.*y + 24576.*x.*y.^5 + 24576.*x.^5.*y + 12154.*x.^2 - 25080.*x.^3 + 28480.*x.^4 - 16896.*x.^5 + 4096.*x.^6 + 12154.*y.^2 - 25080.*y.^3 + 28480.*y.^4 - 16896.*y.^5 + 4096.*y.^6 + 315))/45;
            case 26
                z = -(256.*x.*y.*(x + y - 1).*(274.*x - 1800.*x.^2 + 5440.*x.^3 - 7680.*x.^4 + 4096.*x.^5 - 15))/45;
            case 27
                z = -(256.*x.*y.*(x + y - 1).*(274.*y - 1800.*y.^2 + 5440.*y.^3 - 7680.*y.^4 + 4096.*y.^5 - 15))/45;
            case 28
                z = -(256.*x.*y.*(8.*y - 1).*(x + y - 1).*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3))/15;
            case 29
                z = -(128.*x.*y.*(x + y - 1).*(32.*y.^2 - 12.*y + 1).*(44.*x - 192.*x.^2 + 256.*x.^3 - 3))/9;
            case 30
                z = -(128.*x.*y.*(x + y - 1).*(32.*x.^2 - 12.*x + 1).*(44.*y - 192.*y.^2 + 256.*y.^3 - 3))/9;
            case 31
                z = -(256.*x.*y.*(8.*x - 1).*(x + y - 1).*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3))/15;
            case 32
                z = (256.*x.*y.*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/15;
            case 33
                z = -(128.*x.*y.*(44.*y - 192.*y.^2 + 256.*y.^3 - 3).*(73.*x + 73.*y - 168.*x.*y + 96.*x.*y.^2 + 96.*x.^2.*y - 84.*x.^2 + 32.*x.^3 - 84.*y.^2 + 32.*y.^3 - 21))/9;
            case 34
                z = (128.*x.*y.*(32.*y.^2 - 12.*y + 1).*(1536.*x.^2.*y.^2 - 533.*y - 533.*x + 2008.*x.*y - 2496.*x.*y.^2 - 2496.*x.^2.*y + 1024.*x.*y.^3 + 1024.*x.^3.*y + 1004.*x.^2 - 832.*x.^3 + 256.*x.^4 + 1004.*y.^2 - 832.*y.^3 + 256.*y.^4 + 105))/9;
            case 35
                z = -(256.*x.*y.*(8.*y - 1).*(743.*x + 743.*y - 11520.*x.^2.*y.^2 + 5120.*x.^2.*y.^3 + 5120.*x.^3.*y.^2 - 4140.*x.*y + 8520.*x.*y.^2 + 8520.*x.^2.*y - 7680.*x.*y.^3 - 7680.*x.^3.*y + 2560.*x.*y.^4 + 2560.*x.^4.*y - 2070.*x.^2 + 2840.*x.^3 - 1920.*x.^4 + 512.*x.^5 - 2070.*y.^2 + 2840.*y.^3 - 1920.*y.^4 + 512.*y.^5 - 105))/15;
            case 36
                z = -(256.*x.*y.*(8.*x - 1).*(743.*x + 743.*y - 11520.*x.^2.*y.^2 + 5120.*x.^2.*y.^3 + 5120.*x.^3.*y.^2 - 4140.*x.*y + 8520.*x.*y.^2 + 8520.*x.^2.*y - 7680.*x.*y.^3 - 7680.*x.^3.*y + 2560.*x.*y.^4 + 2560.*x.^4.*y - 2070.*x.^2 + 2840.*x.^3 - 1920.*x.^4 + 512.*x.^5 - 2070.*y.^2 + 2840.*y.^3 - 1920.*y.^4 + 512.*y.^5 - 105))/15;
            case 37
                z = (128.*x.*y.*(32.*x.^2 - 12.*x + 1).*(1536.*x.^2.*y.^2 - 533.*y - 533.*x + 2008.*x.*y - 2496.*x.*y.^2 - 2496.*x.^2.*y + 1024.*x.*y.^3 + 1024.*x.^3.*y + 1004.*x.^2 - 832.*x.^3 + 256.*x.^4 + 1004.*y.^2 - 832.*y.^3 + 256.*y.^4 + 105))/9;
            case 38
                z = -(128.*x.*y.*(44.*x - 192.*x.^2 + 256.*x.^3 - 3).*(73.*x + 73.*y - 168.*x.*y + 96.*x.*y.^2 + 96.*x.^2.*y - 84.*x.^2 + 32.*x.^3 - 84.*y.^2 + 32.*y.^3 - 21))/9;
            case 39
                z = (256.*x.*y.*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/15;
            case 40
                z = (32.*x.*y.*(8.*x - 1).*(8.*y - 1).*(1536.*x.^2.*y.^2 - 533.*y - 533.*x + 2008.*x.*y - 2496.*x.*y.^2 - 2496.*x.^2.*y + 1024.*x.*y.^3 + 1024.*x.^3.*y + 1004.*x.^2 - 832.*x.^3 + 256.*x.^4 + 1004.*y.^2 - 832.*y.^3 + 256.*y.^4 + 105))/3;
            case 41
                z = (32.*x.*y.*(8.*y - 1).*(44.*x - 192.*x.^2 + 256.*x.^3 - 3).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/3;
            case 42
                z = (32.*x.*y.*(8.*x - 1).*(44.*y - 192.*y.^2 + 256.*y.^3 - 3).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/3;
            case 43
                z = (256.*x.*y.*(32.*x.^2 - 12.*x + 1).*(32.*y.^2 - 12.*y + 1).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/9;
            case 44
                z = -(256.*x.*y.*(8.*x - 1).*(32.*y.^2 - 12.*y + 1).*(73.*x + 73.*y - 168.*x.*y + 96.*x.*y.^2 + 96.*x.^2.*y - 84.*x.^2 + 32.*x.^3 - 84.*y.^2 + 32.*y.^3 - 21))/9;
            case 45
                z = -(256.*x.*y.*(8.*y - 1).*(32.*x.^2 - 12.*x + 1).*(73.*x + 73.*y - 168.*x.*y + 96.*x.*y.^2 + 96.*x.^2.*y - 84.*x.^2 + 32.*x.^3 - 84.*y.^2 + 32.*y.^3 - 21))/9;
        end
    case 9
        switch i
            case 1
                z = (1869885.*x.^2.*y.^2)/64 - (7129.*y)/280 - (7129.*x)/280 - (6589431.*x.^2.*y.^3)/64 - (6589431.*x.^3.*y.^2)/64 + (13286025.*x.^2.*y.^4)/64 + (4428675.*x.^3.*y.^3)/16 + (13286025.*x.^4.*y.^2)/64 - (15411789.*x.^2.*y.^5)/64 - (25686315.*x.^3.*y.^4)/64 - (25686315.*x.^4.*y.^3)/64 - (15411789.*x.^5.*y.^2)/64 + (4782969.*x.^2.*y.^6)/32 + (4782969.*x.^3.*y.^5)/16 + (23914845.*x.^4.*y.^4)/64 + (4782969.*x.^5.*y.^3)/16 + (4782969.*x.^6.*y.^2)/32 - (43046721.*x.^2.*y.^7)/1120 - (14348907.*x.^3.*y.^6)/160 - (43046721.*x.^4.*y.^5)/320 - (43046721.*x.^5.*y.^4)/320 - (14348907.*x.^6.*y.^3)/160 - (43046721.*x.^7.*y.^2)/1120 + (58635.*x.*y)/112 - (122121.*x.*y.^2)/28 - (122121.*x.^2.*y)/28 + (623295.*x.*y.^3)/32 + (623295.*x.^3.*y)/32 - (6589431.*x.*y.^4)/128 - (6589431.*x.^4.*y)/128 + (2657205.*x.*y.^5)/32 + (2657205.*x.^5.*y)/32 - (5137263.*x.*y.^6)/64 - (5137263.*x.^6.*y)/64 + (4782969.*x.*y.^7)/112 + (4782969.*x.^7.*y)/112 - (43046721.*x.*y.^8)/4480 - (43046721.*x.^8.*y)/4480 + (58635.*x.^2)/224 - (40707.*x.^3)/28 + (623295.*x.^4)/128 - (6589431.*x.^5)/640 + (885735.*x.^6)/64 - (5137263.*x.^7)/448 + (4782969.*x.^8)/896 - (4782969.*x.^9)/4480 + (58635.*y.^2)/224 - (40707.*y.^3)/28 + (623295.*y.^4)/128 - (6589431.*y.^5)/640 + (885735.*y.^6)/64 - (5137263.*y.^7)/448 + (4782969.*y.^8)/896 - (4782969.*y.^9)/4480 + 1;
            case 2
                z = (x.*(1063116.*x.^2 - 109584.*x - 5450004.*x.^3 + 16365321.*x.^4 - 29760696.*x.^5 + 32240754.*x.^6 - 19131876.*x.^7 + 4782969.*x.^8 + 4480))/4480;
            case 3
                z = (y.*(1063116.*y.^2 - 109584.*y - 5450004.*y.^3 + 16365321.*y.^4 - 29760696.*y.^5 + 32240754.*y.^6 - 19131876.*y.^7 + 4782969.*y.^8 + 4480))/4480;
            case 4
                z = (81.*x.*y.*(13068.*x - 118188.*x.^2 + 548289.*x.^3 - 1428840.*x.^4 + 2112642.*x.^5 - 1653372.*x.^6 + 531441.*x.^7 - 560))/4480;
            case 5
                z = (81.*x.*y.*(9.*y - 1).*(14616.*x.^2 - 1764.*x - 59535.*x.^3 + 127575.*x.^4 - 137781.*x.^5 + 59049.*x.^6 + 80))/1120;
            case 6
                z = (9.*x.*y.*(81.*y.^2 - 27.*y + 2).*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40))/160;
            case 7
                z = (81.*x.*y.*(33.*y - 162.*y.^2 + 243.*y.^3 - 2).*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8))/320;
            case 8
                z = (81.*x.*y.*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8))/320;
            case 9
                z = (9.*x.*y.*(81.*x.^2 - 27.*x + 2).*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40))/160;
            case 10
                z = (81.*x.*y.*(9.*x - 1).*(14616.*y.^2 - 1764.*y - 59535.*y.^3 + 127575.*y.^4 - 137781.*y.^5 + 59049.*y.^6 + 80))/1120;
            case 11
                z = (81.*x.*y.*(13068.*y - 118188.*y.^2 + 548289.*y.^3 - 1428840.*y.^4 + 2112642.*y.^5 - 1653372.*y.^6 + 531441.*y.^7 - 560))/4480;
            case 12
                z = -(81.*y.*(x + y - 1).*(13068.*y - 118188.*y.^2 + 548289.*y.^3 - 1428840.*y.^4 + 2112642.*y.^5 - 1653372.*y.^6 + 531441.*y.^7 - 560))/4480;
            case 13
                z = (81.*y.*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8).*(14616.*y.^2 - 1764.*y - 59535.*y.^3 + 127575.*y.^4 - 137781.*y.^5 + 59049.*y.^6 + 80))/1120;
            case 14
                z = -(9.*y.*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/160;
            case 15
                z = (81.*y.*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8).*(1458.*x.^2.*y.^2 - 550.*y - 550.*x + 2010.*x.*y - 2430.*x.*y.^2 - 2430.*x.^2.*y + 972.*x.*y.^3 + 972.*x.^3.*y + 1005.*x.^2 - 810.*x.^3 + 243.*x.^4 + 1005.*y.^2 - 810.*y.^3 + 243.*y.^4 + 112))/320;
            case 16
                z = -(81.*y.*(33.*y - 162.*y.^2 + 243.*y.^3 - 2).*(3758.*x + 3758.*y - 51030.*x.^2.*y.^2 + 21870.*x.^2.*y.^3 + 21870.*x.^3.*y.^2 - 19950.*x.*y + 39285.*x.*y.^2 + 39285.*x.^2.*y - 34020.*x.*y.^3 - 34020.*x.^3.*y + 10935.*x.*y.^4 + 10935.*x.^4.*y - 9975.*x.^2 + 13095.*x.^3 - 8505.*x.^4 + 2187.*x.^5 - 9975.*y.^2 + 13095.*y.^3 - 8505.*y.^4 + 2187.*y.^5 - 560))/320;
            case 17
                z = (9.*y.*(81.*y.^2 - 27.*y + 2).*(911250.*x.^2.*y.^2 - 20072.*y - 20072.*x - 852930.*x.^2.*y.^3 - 852930.*x.^3.*y.^2 + 295245.*x.^2.*y.^4 + 393660.*x.^3.*y.^3 + 295245.*x.^4.*y.^2 + 147444.*x.*y - 426465.*x.*y.^2 - 426465.*x.^2.*y + 607500.*x.*y.^3 + 607500.*x.^3.*y - 426465.*x.*y.^4 - 426465.*x.^4.*y + 118098.*x.*y.^5 + 118098.*x.^5.*y + 73722.*x.^2 - 142155.*x.^3 + 151875.*x.^4 - 85293.*x.^5 + 19683.*x.^6 + 73722.*y.^2 - 142155.*y.^3 + 151875.*y.^4 - 85293.*y.^5 + 19683.*y.^6 + 2240))/160;
            case 18
                z = -(81.*y.*(9.*y - 1).*(26792.*x + 26792.*y - 3470040.*x.^2.*y.^2 + 5409180.*x.^2.*y.^3 + 5409180.*x.^3.*y.^2 - 4133430.*x.^2.*y.^4 - 5511240.*x.^3.*y.^3 - 4133430.*x.^4.*y.^2 + 1240029.*x.^2.*y.^5 + 2066715.*x.^3.*y.^4 + 2066715.*x.^4.*y.^3 + 1240029.*x.^5.*y.^2 - 267876.*x.*y + 1089963.*x.*y.^2 + 1089963.*x.^2.*y - 2313360.*x.*y.^3 - 2313360.*x.^3.*y + 2704590.*x.*y.^4 + 2704590.*x.^4.*y - 1653372.*x.*y.^5 - 1653372.*x.^5.*y + 413343.*x.*y.^6 + 413343.*x.^6.*y - 133938.*x.^2 + 363321.*x.^3 - 578340.*x.^4 + 540918.*x.^5 - 275562.*x.^6 + 59049.*x.^7 - 133938.*y.^2 + 363321.*y.^3 - 578340.*y.^4 + 540918.*y.^5 - 275562.*y.^6 + 59049.*y.^7 - 2240))/1120;
            case 19
                z = (81.*y.*(26559414.*x.^2.*y.^2 - 73744.*y - 73744.*x - 62868960.*x.^2.*y.^3 - 62868960.*x.^3.*y.^2 + 81290790.*x.^2.*y.^4 + 108387720.*x.^3.*y.^3 + 81290790.*x.^4.*y.^2 - 54561276.*x.^2.*y.^5 - 90935460.*x.^3.*y.^4 - 90935460.*x.^4.*y.^3 - 54561276.*x.^5.*y.^2 + 14880348.*x.^2.*y.^6 + 29760696.*x.^3.*y.^5 + 37200870.*x.^4.*y.^4 + 29760696.*x.^5.*y.^3 + 14880348.*x.^6.*y.^2 + 1018008.*x.*y - 5796252.*x.*y.^2 - 5796252.*x.^2.*y + 17706276.*x.*y.^3 + 17706276.*x.^3.*y - 31434480.*x.*y.^4 - 31434480.*x.^4.*y + 32516316.*x.*y.^5 + 32516316.*x.^5.*y - 18187092.*x.*y.^6 - 18187092.*x.^6.*y + 4251528.*x.*y.^7 + 4251528.*x.^7.*y + 509004.*x.^2 - 1932084.*x.^3 + 4426569.*x.^4 - 6286896.*x.^5 + 5419386.*x.^6 - 2598156.*x.^7 + 531441.*x.^8 + 509004.*y.^2 - 1932084.*y.^3 + 4426569.*y.^4 - 6286896.*y.^5 + 5419386.*y.^6 - 2598156.*y.^7 + 531441.*y.^8 + 4480))/4480;
            case 20
                z = (81.*x.*(26559414.*x.^2.*y.^2 - 73744.*y - 73744.*x - 62868960.*x.^2.*y.^3 - 62868960.*x.^3.*y.^2 + 81290790.*x.^2.*y.^4 + 108387720.*x.^3.*y.^3 + 81290790.*x.^4.*y.^2 - 54561276.*x.^2.*y.^5 - 90935460.*x.^3.*y.^4 - 90935460.*x.^4.*y.^3 - 54561276.*x.^5.*y.^2 + 14880348.*x.^2.*y.^6 + 29760696.*x.^3.*y.^5 + 37200870.*x.^4.*y.^4 + 29760696.*x.^5.*y.^3 + 14880348.*x.^6.*y.^2 + 1018008.*x.*y - 5796252.*x.*y.^2 - 5796252.*x.^2.*y + 17706276.*x.*y.^3 + 17706276.*x.^3.*y - 31434480.*x.*y.^4 - 31434480.*x.^4.*y + 32516316.*x.*y.^5 + 32516316.*x.^5.*y - 18187092.*x.*y.^6 - 18187092.*x.^6.*y + 4251528.*x.*y.^7 + 4251528.*x.^7.*y + 509004.*x.^2 - 1932084.*x.^3 + 4426569.*x.^4 - 6286896.*x.^5 + 5419386.*x.^6 - 2598156.*x.^7 + 531441.*x.^8 + 509004.*y.^2 - 1932084.*y.^3 + 4426569.*y.^4 - 6286896.*y.^5 + 5419386.*y.^6 - 2598156.*y.^7 + 531441.*y.^8 + 4480))/4480;
            case 21
                z = -(81.*x.*(9.*x - 1).*(26792.*x + 26792.*y - 3470040.*x.^2.*y.^2 + 5409180.*x.^2.*y.^3 + 5409180.*x.^3.*y.^2 - 4133430.*x.^2.*y.^4 - 5511240.*x.^3.*y.^3 - 4133430.*x.^4.*y.^2 + 1240029.*x.^2.*y.^5 + 2066715.*x.^3.*y.^4 + 2066715.*x.^4.*y.^3 + 1240029.*x.^5.*y.^2 - 267876.*x.*y + 1089963.*x.*y.^2 + 1089963.*x.^2.*y - 2313360.*x.*y.^3 - 2313360.*x.^3.*y + 2704590.*x.*y.^4 + 2704590.*x.^4.*y - 1653372.*x.*y.^5 - 1653372.*x.^5.*y + 413343.*x.*y.^6 + 413343.*x.^6.*y - 133938.*x.^2 + 363321.*x.^3 - 578340.*x.^4 + 540918.*x.^5 - 275562.*x.^6 + 59049.*x.^7 - 133938.*y.^2 + 363321.*y.^3 - 578340.*y.^4 + 540918.*y.^5 - 275562.*y.^6 + 59049.*y.^7 - 2240))/1120;
            case 22
                z = (9.*x.*(81.*x.^2 - 27.*x + 2).*(911250.*x.^2.*y.^2 - 20072.*y - 20072.*x - 852930.*x.^2.*y.^3 - 852930.*x.^3.*y.^2 + 295245.*x.^2.*y.^4 + 393660.*x.^3.*y.^3 + 295245.*x.^4.*y.^2 + 147444.*x.*y - 426465.*x.*y.^2 - 426465.*x.^2.*y + 607500.*x.*y.^3 + 607500.*x.^3.*y - 426465.*x.*y.^4 - 426465.*x.^4.*y + 118098.*x.*y.^5 + 118098.*x.^5.*y + 73722.*x.^2 - 142155.*x.^3 + 151875.*x.^4 - 85293.*x.^5 + 19683.*x.^6 + 73722.*y.^2 - 142155.*y.^3 + 151875.*y.^4 - 85293.*y.^5 + 19683.*y.^6 + 2240))/160;
            case 23
                z = -(81.*x.*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(3758.*x + 3758.*y - 51030.*x.^2.*y.^2 + 21870.*x.^2.*y.^3 + 21870.*x.^3.*y.^2 - 19950.*x.*y + 39285.*x.*y.^2 + 39285.*x.^2.*y - 34020.*x.*y.^3 - 34020.*x.^3.*y + 10935.*x.*y.^4 + 10935.*x.^4.*y - 9975.*x.^2 + 13095.*x.^3 - 8505.*x.^4 + 2187.*x.^5 - 9975.*y.^2 + 13095.*y.^3 - 8505.*y.^4 + 2187.*y.^5 - 560))/320;
            case 24
                z = (81.*x.*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8).*(1458.*x.^2.*y.^2 - 550.*y - 550.*x + 2010.*x.*y - 2430.*x.*y.^2 - 2430.*x.^2.*y + 972.*x.*y.^3 + 972.*x.^3.*y + 1005.*x.^2 - 810.*x.^3 + 243.*x.^4 + 1005.*y.^2 - 810.*y.^3 + 243.*y.^4 + 112))/320;
            case 25
                z = -(9.*x.*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/160;
            case 26
                z = (81.*x.*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8).*(14616.*x.^2 - 1764.*x - 59535.*x.^3 + 127575.*x.^4 - 137781.*x.^5 + 59049.*x.^6 + 80))/1120;
            case 27
                z = -(81.*x.*(x + y - 1).*(13068.*x - 118188.*x.^2 + 548289.*x.^3 - 1428840.*x.^4 + 2112642.*x.^5 - 1653372.*x.^6 + 531441.*x.^7 - 560))/4480;
            case 28
                z = -(729.*x.*y.*(26792.*x + 26792.*y - 3470040.*x.^2.*y.^2 + 5409180.*x.^2.*y.^3 + 5409180.*x.^3.*y.^2 - 4133430.*x.^2.*y.^4 - 5511240.*x.^3.*y.^3 - 4133430.*x.^4.*y.^2 + 1240029.*x.^2.*y.^5 + 2066715.*x.^3.*y.^4 + 2066715.*x.^4.*y.^3 + 1240029.*x.^5.*y.^2 - 267876.*x.*y + 1089963.*x.*y.^2 + 1089963.*x.^2.*y - 2313360.*x.*y.^3 - 2313360.*x.^3.*y + 2704590.*x.*y.^4 + 2704590.*x.^4.*y - 1653372.*x.*y.^5 - 1653372.*x.^5.*y + 413343.*x.*y.^6 + 413343.*x.^6.*y - 133938.*x.^2 + 363321.*x.^3 - 578340.*x.^4 + 540918.*x.^5 - 275562.*x.^6 + 59049.*x.^7 - 133938.*y.^2 + 363321.*y.^3 - 578340.*y.^4 + 540918.*y.^5 - 275562.*y.^6 + 59049.*y.^7 - 2240))/560;
            case 29
                z = -(729.*x.*y.*(x + y - 1).*(14616.*x.^2 - 1764.*x - 59535.*x.^3 + 127575.*x.^4 - 137781.*x.^5 + 59049.*x.^6 + 80))/560;
            case 30
                z = -(729.*x.*y.*(x + y - 1).*(14616.*y.^2 - 1764.*y - 59535.*y.^3 + 127575.*y.^4 - 137781.*y.^5 + 59049.*y.^6 + 80))/560;
            case 31
                z = -(243.*x.*y.*(9.*y - 1).*(x + y - 1).*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40))/160;
            case 32
                z = -(243.*x.*y.*(x + y - 1).*(81.*y.^2 - 27.*y + 2).*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8))/80;
            case 33
                z = -(729.*x.*y.*(x + y - 1).*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(33.*y - 162.*y.^2 + 243.*y.^3 - 2))/64;
            case 34
                z = -(243.*x.*y.*(x + y - 1).*(81.*x.^2 - 27.*x + 2).*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8))/80;
            case 35
                z = -(243.*x.*y.*(9.*x - 1).*(x + y - 1).*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40))/160;
            case 36
                z = (243.*x.*y.*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160;
            case 37
                z = -(243.*x.*y.*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/80;
            case 38
                z = (729.*x.*y.*(33.*y - 162.*y.^2 + 243.*y.^3 - 2).*(1458.*x.^2.*y.^2 - 550.*y - 550.*x + 2010.*x.*y - 2430.*x.*y.^2 - 2430.*x.^2.*y + 972.*x.*y.^3 + 972.*x.^3.*y + 1005.*x.^2 - 810.*x.^3 + 243.*x.^4 + 1005.*y.^2 - 810.*y.^3 + 243.*y.^4 + 112))/64;
            case 39
                z = -(243.*x.*y.*(81.*y.^2 - 27.*y + 2).*(3758.*x + 3758.*y - 51030.*x.^2.*y.^2 + 21870.*x.^2.*y.^3 + 21870.*x.^3.*y.^2 - 19950.*x.*y + 39285.*x.*y.^2 + 39285.*x.^2.*y - 34020.*x.*y.^3 - 34020.*x.^3.*y + 10935.*x.*y.^4 + 10935.*x.^4.*y - 9975.*x.^2 + 13095.*x.^3 - 8505.*x.^4 + 2187.*x.^5 - 9975.*y.^2 + 13095.*y.^3 - 8505.*y.^4 + 2187.*y.^5 - 560))/80;
            case 40
                z = (243.*x.*y.*(9.*y - 1).*(911250.*x.^2.*y.^2 - 20072.*y - 20072.*x - 852930.*x.^2.*y.^3 - 852930.*x.^3.*y.^2 + 295245.*x.^2.*y.^4 + 393660.*x.^3.*y.^3 + 295245.*x.^4.*y.^2 + 147444.*x.*y - 426465.*x.*y.^2 - 426465.*x.^2.*y + 607500.*x.*y.^3 + 607500.*x.^3.*y - 426465.*x.*y.^4 - 426465.*x.^4.*y + 118098.*x.*y.^5 + 118098.*x.^5.*y + 73722.*x.^2 - 142155.*x.^3 + 151875.*x.^4 - 85293.*x.^5 + 19683.*x.^6 + 73722.*y.^2 - 142155.*y.^3 + 151875.*y.^4 - 85293.*y.^5 + 19683.*y.^6 + 2240))/160;
            case 41
                z = (243.*x.*y.*(9.*x - 1).*(911250.*x.^2.*y.^2 - 20072.*y - 20072.*x - 852930.*x.^2.*y.^3 - 852930.*x.^3.*y.^2 + 295245.*x.^2.*y.^4 + 393660.*x.^3.*y.^3 + 295245.*x.^4.*y.^2 + 147444.*x.*y - 426465.*x.*y.^2 - 426465.*x.^2.*y + 607500.*x.*y.^3 + 607500.*x.^3.*y - 426465.*x.*y.^4 - 426465.*x.^4.*y + 118098.*x.*y.^5 + 118098.*x.^5.*y + 73722.*x.^2 - 142155.*x.^3 + 151875.*x.^4 - 85293.*x.^5 + 19683.*x.^6 + 73722.*y.^2 - 142155.*y.^3 + 151875.*y.^4 - 85293.*y.^5 + 19683.*y.^6 + 2240))/160;
            case 42
                z = -(243.*x.*y.*(81.*x.^2 - 27.*x + 2).*(3758.*x + 3758.*y - 51030.*x.^2.*y.^2 + 21870.*x.^2.*y.^3 + 21870.*x.^3.*y.^2 - 19950.*x.*y + 39285.*x.*y.^2 + 39285.*x.^2.*y - 34020.*x.*y.^3 - 34020.*x.^3.*y + 10935.*x.*y.^4 + 10935.*x.^4.*y - 9975.*x.^2 + 13095.*x.^3 - 8505.*x.^4 + 2187.*x.^5 - 9975.*y.^2 + 13095.*y.^3 - 8505.*y.^4 + 2187.*y.^5 - 560))/80;
            case 43
                z = (729.*x.*y.*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(1458.*x.^2.*y.^2 - 550.*y - 550.*x + 2010.*x.*y - 2430.*x.*y.^2 - 2430.*x.^2.*y + 972.*x.*y.^3 + 972.*x.^3.*y + 1005.*x.^2 - 810.*x.^3 + 243.*x.^4 + 1005.*y.^2 - 810.*y.^3 + 243.*y.^4 + 112))/64;
            case 44
                z = -(243.*x.*y.*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/80;
            case 45
                z = (243.*x.*y.*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160;
            case 46
                z = -(729.*x.*y.*(9.*x - 1).*(9.*y - 1).*(3758.*x + 3758.*y - 51030.*x.^2.*y.^2 + 21870.*x.^2.*y.^3 + 21870.*x.^3.*y.^2 - 19950.*x.*y + 39285.*x.*y.^2 + 39285.*x.^2.*y - 34020.*x.*y.^3 - 34020.*x.^3.*y + 10935.*x.*y.^4 + 10935.*x.^4.*y - 9975.*x.^2 + 13095.*x.^3 - 8505.*x.^4 + 2187.*x.^5 - 9975.*y.^2 + 13095.*y.^3 - 8505.*y.^4 + 2187.*y.^5 - 560))/160;
            case 47
                z = (729.*x.*y.*(9.*y - 1).*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160;
            case 48
                z = (729.*x.*y.*(9.*x - 1).*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160;
            case 49
                z = (243.*x.*y.*(81.*y.^2 - 27.*y + 2).*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/32;
            case 50
                z = (243.*x.*y.*(81.*x.^2 - 27.*x + 2).*(33.*y - 162.*y.^2 + 243.*y.^3 - 2).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/32;
            case 51
                z = -(243.*x.*y.*(9.*x - 1).*(33.*y - 162.*y.^2 + 243.*y.^3 - 2).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/32;
            case 52
                z = (243.*x.*y.*(9.*x - 1).*(81.*y.^2 - 27.*y + 2).*(1458.*x.^2.*y.^2 - 550.*y - 550.*x + 2010.*x.*y - 2430.*x.*y.^2 - 2430.*x.^2.*y + 972.*x.*y.^3 + 972.*x.^3.*y + 1005.*x.^2 - 810.*x.^3 + 243.*x.^4 + 1005.*y.^2 - 810.*y.^3 + 243.*y.^4 + 112))/32;
            case 53
                z = (243.*x.*y.*(9.*y - 1).*(81.*x.^2 - 27.*x + 2).*(1458.*x.^2.*y.^2 - 550.*y - 550.*x + 2010.*x.*y - 2430.*x.*y.^2 - 2430.*x.^2.*y + 972.*x.*y.^3 + 972.*x.^3.*y + 1005.*x.^2 - 810.*x.^3 + 243.*x.^4 + 1005.*y.^2 - 810.*y.^3 + 243.*y.^4 + 112))/32;
            case 54
                z = -(243.*x.*y.*(9.*y - 1).*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/32;
            case 55
                z = -(27.*x.*y.*(81.*x.^2 - 27.*x + 2).*(81.*y.^2 - 27.*y + 2).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/8;
        end
    case 10
        switch i
            case 1
                z = (42711625.*x.^2.*y.^2)/756 - (7381.*y)/252 - (7381.*x)/252 - (26846875.*x.^2.*y.^3)/108 - (26846875.*x.^3.*y.^2)/108 + (23478125.*x.^2.*y.^4)/36 + (23478125.*x.^3.*y.^3)/27 + (23478125.*x.^4.*y.^2)/36 - (9453125.*x.^2.*y.^5)/9 - (47265625.*x.^3.*y.^4)/27 - (47265625.*x.^4.*y.^3)/27 - (9453125.*x.^5.*y.^2)/9 + (27500000.*x.^2.*y.^6)/27 + (55000000.*x.^3.*y.^5)/27 + (68750000.*x.^4.*y.^4)/27 + (55000000.*x.^5.*y.^3)/27 + (27500000.*x.^6.*y.^2)/27 - (34375000.*x.^2.*y.^7)/63 - (34375000.*x.^3.*y.^6)/27 - (17187500.*x.^4.*y.^5)/9 - (17187500.*x.^5.*y.^4)/9 - (34375000.*x.^6.*y.^3)/27 - (34375000.*x.^7.*y.^2)/63 + (7812500.*x.^2.*y.^8)/63 + (62500000.*x.^3.*y.^7)/189 + (15625000.*x.^4.*y.^6)/27 + (6250000.*x.^5.*y.^5)/9 + (15625000.*x.^6.*y.^4)/27 + (62500000.*x.^7.*y.^3)/189 + (7812500.*x.^8.*y.^2)/63 + (177133.*x.*y)/252 - (10511875.*x.*y.^2)/1512 - (10511875.*x.^2.*y)/1512 + (42711625.*x.*y.^3)/1134 + (42711625.*x.^3.*y)/1134 - (26846875.*x.*y.^4)/216 - (26846875.*x.^4.*y)/216 + (4695625.*x.*y.^5)/18 + (4695625.*x.^5.*y)/18 - (9453125.*x.*y.^6)/27 - (9453125.*x.^6.*y)/27 + (55000000.*x.*y.^7)/189 + (55000000.*x.^7.*y)/189 - (8593750.*x.*y.^8)/63 - (8593750.*x.^8.*y)/63 + (15625000.*x.*y.^9)/567 + (15625000.*x.^9.*y)/567 + (177133.*x.^2)/504 - (10511875.*x.^3)/4536 + (42711625.*x.^4)/4536 - (5369375.*x.^5)/216 + (4695625.*x.^6)/108 - (9453125.*x.^7)/189 + (6875000.*x.^8)/189 - (8593750.*x.^9)/567 + (1562500.*x.^10)/567 + (177133.*y.^2)/504 - (10511875.*y.^3)/4536 + (42711625.*y.^4)/4536 - (5369375.*y.^5)/216 + (4695625.*y.^6)/108 - (9453125.*y.^7)/189 + (6875000.*y.^8)/189 - (8593750.*y.^9)/567 + (1562500.*y.^10)/567 + 1;
            case 2
                z = (7129.*x.^2)/252 - x - (162875.*x.^3)/504 + (1130750.*x.^4)/567 - (59375.*x.^5)/8 + (1883125.*x.^6)/108 - (78125.*x.^7)/3 + (4531250.*x.^8)/189 - (781250.*x.^9)/63 + (1562500.*x.^10)/567;
            case 3
                z = (7129.*y.^2)/252 - y - (162875.*y.^3)/504 + (1130750.*y.^4)/567 - (59375.*y.^5)/8 + (1883125.*y.^6)/108 - (78125.*y.^7)/3 + (4531250.*y.^8)/189 - (781250.*y.^9)/63 + (1562500.*y.^10)/567;
            case 4
                z = (25.*x.*y.*(147655.*x.^2 - 13698.*x - 841050.*x.^3 + 2806125.*x.^4 - 5670000.*x.^5 + 6825000.*x.^6 - 4500000.*x.^7 + 1250000.*x.^8 + 504))/1134;
            case 5
                z = (25.*x.*y.*(10.*y - 1).*(3267.*x - 32830.*x.^2 + 169225.*x.^3 - 490000.*x.^4 + 805000.*x.^5 - 700000.*x.^6 + 250000.*x.^7 - 126))/504;
            case 6
                z = (50.*x.*y.*(50.*y.^2 - 15.*y + 1).*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18))/189;
            case 7
                z = (25.*x.*y.*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6))/108;
            case 8
                z = (x.*y.*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6))/9;
            case 9
                z = (25.*x.*y.*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6))/108;
            case 10
                z = (50.*x.*y.*(50.*x.^2 - 15.*x + 1).*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18))/189;
            case 11
                z = (25.*x.*y.*(10.*x - 1).*(3267.*y - 32830.*y.^2 + 169225.*y.^3 - 490000.*y.^4 + 805000.*y.^5 - 700000.*y.^6 + 250000.*y.^7 - 126))/504;
            case 12
                z = (25.*x.*y.*(147655.*y.^2 - 13698.*y - 841050.*y.^3 + 2806125.*y.^4 - 5670000.*y.^5 + 6825000.*y.^6 - 4500000.*y.^7 + 1250000.*y.^8 + 504))/1134;
            case 13
                z = -(25.*y.*(x + y - 1).*(147655.*y.^2 - 13698.*y - 841050.*y.^3 + 2806125.*y.^4 - 5670000.*y.^5 + 6825000.*y.^6 - 4500000.*y.^7 + 1250000.*y.^8 + 504))/1134;
            case 14
                z = (25.*y.*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9).*(3267.*y - 32830.*y.^2 + 169225.*y.^3 - 490000.*y.^4 + 805000.*y.^5 - 700000.*y.^6 + 250000.*y.^7 - 126))/504;
            case 15
                z = -(50.*y.*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/189;
            case 16
                z = (25.*y.*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(3000.*x.^2.*y.^2 - 1207.*y - 1207.*x + 4310.*x.*y - 5100.*x.*y.^2 - 5100.*x.^2.*y + 2000.*x.*y.^3 + 2000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/108;
            case 17
                z = -(y.*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(4881.*x + 4881.*y - 60000.*x.^2.*y.^2 + 25000.*x.^2.*y.^3 + 25000.*x.^3.*y.^2 - 25000.*x.*y + 47625.*x.*y.^2 + 47625.*x.^2.*y - 40000.*x.*y.^3 - 40000.*x.^3.*y + 12500.*x.*y.^4 + 12500.*x.^4.*y - 12500.*x.^2 + 15875.*x.^3 - 10000.*x.^4 + 2500.*x.^5 - 12500.*y.^2 + 15875.*y.^3 - 10000.*y.^4 + 2500.*y.^5 - 756))/9;
            case 18
                z = (25.*y.*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(250500.*x.^2.*y.^2 - 6393.*y - 6393.*x - 225000.*x.^2.*y.^3 - 225000.*x.^3.*y.^2 + 75000.*x.^2.*y.^4 + 100000.*x.^3.*y.^3 + 75000.*x.^4.*y.^2 + 44524.*x.*y - 122625.*x.*y.^2 - 122625.*x.^2.*y + 167000.*x.*y.^3 + 167000.*x.^3.*y - 112500.*x.*y.^4 - 112500.*x.^4.*y + 30000.*x.*y.^5 + 30000.*x.^5.*y + 22262.*x.^2 - 40875.*x.^3 + 41750.*x.^4 - 22500.*x.^5 + 5000.*x.^6 + 22262.*y.^2 - 40875.*y.^3 + 41750.*y.^4 - 22500.*y.^5 + 5000.*y.^6 + 756))/108;
            case 19
                z = -(50.*y.*(50.*y.^2 - 15.*y + 1).*(16566.*x + 16566.*y - 1727250.*x.^2.*y.^2 + 2537500.*x.^2.*y.^3 + 2537500.*x.^3.*y.^2 - 1837500.*x.^2.*y.^4 - 2450000.*x.^3.*y.^3 - 1837500.*x.^4.*y.^2 + 525000.*x.^2.*y.^5 + 875000.*x.^3.*y.^4 + 875000.*x.^4.*y.^3 + 525000.*x.^5.*y.^2 - 152978.*x.*y + 579180.*x.*y.^2 + 579180.*x.^2.*y - 1151500.*x.*y.^3 - 1151500.*x.^3.*y + 1268750.*x.*y.^4 + 1268750.*x.^4.*y - 735000.*x.*y.^5 - 735000.*x.^5.*y + 175000.*x.*y.^6 + 175000.*x.^6.*y - 76489.*x.^2 + 193060.*x.^3 - 287875.*x.^4 + 253750.*x.^5 - 122500.*x.^6 + 25000.*x.^7 - 76489.*y.^2 + 193060.*y.^3 - 287875.*y.^4 + 253750.*y.^5 - 122500.*y.^6 + 25000.*y.^7 - 1512))/189;
            case 20
                z = (25.*y.*(10.*y - 1).*(16765350.*x.^2.*y.^2 - 64818.*y - 64818.*x - 36400000.*x.^2.*y.^3 - 36400000.*x.^3.*y.^2 + 43575000.*x.^2.*y.^4 + 58100000.*x.^3.*y.^3 + 43575000.*x.^4.*y.^2 - 27300000.*x.^2.*y.^5 - 45500000.*x.^3.*y.^4 - 45500000.*x.^4.*y.^3 - 27300000.*x.^5.*y.^2 + 7000000.*x.^2.*y.^6 + 14000000.*x.^3.*y.^5 + 17500000.*x.^4.*y.^4 + 14000000.*x.^5.*y.^3 + 7000000.*x.^6.*y.^2 + 790254.*x.*y - 4032210.*x.*y.^2 - 4032210.*x.^2.*y + 11176900.*x.*y.^3 + 11176900.*x.^3.*y - 18200000.*x.*y.^4 - 18200000.*x.^4.*y + 17430000.*x.*y.^5 + 17430000.*x.^5.*y - 9100000.*x.*y.^6 - 9100000.*x.^6.*y + 2000000.*x.*y.^7 + 2000000.*x.^7.*y + 395127.*x.^2 - 1344070.*x.^3 + 2794225.*x.^4 - 3640000.*x.^5 + 2905000.*x.^6 - 1300000.*x.^7 + 250000.*x.^8 + 395127.*y.^2 - 1344070.*y.^3 + 2794225.*y.^4 - 3640000.*y.^5 + 2905000.*y.^6 - 1300000.*y.^7 + 250000.*y.^8 + 4536))/504;
            case 21
                z = -(25.*y.*(87498.*x + 87498.*y - 57087450.*x.^2.*y.^2 + 176111250.*x.^2.*y.^3 + 176111250.*x.^3.*y.^2 - 316575000.*x.^2.*y.^4 - 422100000.*x.^3.*y.^3 - 316575000.*x.^4.*y.^2 + 332325000.*x.^2.*y.^5 + 553875000.*x.^3.*y.^4 + 553875000.*x.^4.*y.^3 + 332325000.*x.^5.*y.^2 - 189000000.*x.^2.*y.^6 - 378000000.*x.^3.*y.^5 - 472500000.*x.^4.*y.^4 - 378000000.*x.^5.*y.^3 - 189000000.*x.^6.*y.^2 + 45000000.*x.^2.*y.^7 + 105000000.*x.^3.*y.^6 + 157500000.*x.^4.*y.^5 + 157500000.*x.^5.*y.^4 + 105000000.*x.^6.*y.^3 + 45000000.*x.^7.*y.^2 - 1438434.*x.*y + 9959115.*x.*y.^2 + 9959115.*x.^2.*y - 38058300.*x.*y.^3 - 38058300.*x.^3.*y + 88055625.*x.*y.^4 + 88055625.*x.^4.*y - 126630000.*x.*y.^5 - 126630000.*x.^5.*y + 110775000.*x.*y.^6 + 110775000.*x.^6.*y - 54000000.*x.*y.^7 - 54000000.*x.^7.*y + 11250000.*x.*y.^8 + 11250000.*x.^8.*y - 719217.*x.^2 + 3319705.*x.^3 - 9514575.*x.^4 + 17611125.*x.^5 - 21105000.*x.^6 + 15825000.*x.^7 - 6750000.*x.^8 + 1250000.*x.^9 - 719217.*y.^2 + 3319705.*y.^3 - 9514575.*y.^4 + 17611125.*y.^5 - 21105000.*y.^6 + 15825000.*y.^7 - 6750000.*y.^8 + 1250000.*y.^9 - 4536))/1134;
            case 22
                z = -(25.*x.*(87498.*x + 87498.*y - 57087450.*x.^2.*y.^2 + 176111250.*x.^2.*y.^3 + 176111250.*x.^3.*y.^2 - 316575000.*x.^2.*y.^4 - 422100000.*x.^3.*y.^3 - 316575000.*x.^4.*y.^2 + 332325000.*x.^2.*y.^5 + 553875000.*x.^3.*y.^4 + 553875000.*x.^4.*y.^3 + 332325000.*x.^5.*y.^2 - 189000000.*x.^2.*y.^6 - 378000000.*x.^3.*y.^5 - 472500000.*x.^4.*y.^4 - 378000000.*x.^5.*y.^3 - 189000000.*x.^6.*y.^2 + 45000000.*x.^2.*y.^7 + 105000000.*x.^3.*y.^6 + 157500000.*x.^4.*y.^5 + 157500000.*x.^5.*y.^4 + 105000000.*x.^6.*y.^3 + 45000000.*x.^7.*y.^2 - 1438434.*x.*y + 9959115.*x.*y.^2 + 9959115.*x.^2.*y - 38058300.*x.*y.^3 - 38058300.*x.^3.*y + 88055625.*x.*y.^4 + 88055625.*x.^4.*y - 126630000.*x.*y.^5 - 126630000.*x.^5.*y + 110775000.*x.*y.^6 + 110775000.*x.^6.*y - 54000000.*x.*y.^7 - 54000000.*x.^7.*y + 11250000.*x.*y.^8 + 11250000.*x.^8.*y - 719217.*x.^2 + 3319705.*x.^3 - 9514575.*x.^4 + 17611125.*x.^5 - 21105000.*x.^6 + 15825000.*x.^7 - 6750000.*x.^8 + 1250000.*x.^9 - 719217.*y.^2 + 3319705.*y.^3 - 9514575.*y.^4 + 17611125.*y.^5 - 21105000.*y.^6 + 15825000.*y.^7 - 6750000.*y.^8 + 1250000.*y.^9 - 4536))/1134;
            case 23
                z = (25.*x.*(10.*x - 1).*(16765350.*x.^2.*y.^2 - 64818.*y - 64818.*x - 36400000.*x.^2.*y.^3 - 36400000.*x.^3.*y.^2 + 43575000.*x.^2.*y.^4 + 58100000.*x.^3.*y.^3 + 43575000.*x.^4.*y.^2 - 27300000.*x.^2.*y.^5 - 45500000.*x.^3.*y.^4 - 45500000.*x.^4.*y.^3 - 27300000.*x.^5.*y.^2 + 7000000.*x.^2.*y.^6 + 14000000.*x.^3.*y.^5 + 17500000.*x.^4.*y.^4 + 14000000.*x.^5.*y.^3 + 7000000.*x.^6.*y.^2 + 790254.*x.*y - 4032210.*x.*y.^2 - 4032210.*x.^2.*y + 11176900.*x.*y.^3 + 11176900.*x.^3.*y - 18200000.*x.*y.^4 - 18200000.*x.^4.*y + 17430000.*x.*y.^5 + 17430000.*x.^5.*y - 9100000.*x.*y.^6 - 9100000.*x.^6.*y + 2000000.*x.*y.^7 + 2000000.*x.^7.*y + 395127.*x.^2 - 1344070.*x.^3 + 2794225.*x.^4 - 3640000.*x.^5 + 2905000.*x.^6 - 1300000.*x.^7 + 250000.*x.^8 + 395127.*y.^2 - 1344070.*y.^3 + 2794225.*y.^4 - 3640000.*y.^5 + 2905000.*y.^6 - 1300000.*y.^7 + 250000.*y.^8 + 4536))/504;
            case 24
                z = -(50.*x.*(50.*x.^2 - 15.*x + 1).*(16566.*x + 16566.*y - 1727250.*x.^2.*y.^2 + 2537500.*x.^2.*y.^3 + 2537500.*x.^3.*y.^2 - 1837500.*x.^2.*y.^4 - 2450000.*x.^3.*y.^3 - 1837500.*x.^4.*y.^2 + 525000.*x.^2.*y.^5 + 875000.*x.^3.*y.^4 + 875000.*x.^4.*y.^3 + 525000.*x.^5.*y.^2 - 152978.*x.*y + 579180.*x.*y.^2 + 579180.*x.^2.*y - 1151500.*x.*y.^3 - 1151500.*x.^3.*y + 1268750.*x.*y.^4 + 1268750.*x.^4.*y - 735000.*x.*y.^5 - 735000.*x.^5.*y + 175000.*x.*y.^6 + 175000.*x.^6.*y - 76489.*x.^2 + 193060.*x.^3 - 287875.*x.^4 + 253750.*x.^5 - 122500.*x.^6 + 25000.*x.^7 - 76489.*y.^2 + 193060.*y.^3 - 287875.*y.^4 + 253750.*y.^5 - 122500.*y.^6 + 25000.*y.^7 - 1512))/189;
            case 25
                z = (25.*x.*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(250500.*x.^2.*y.^2 - 6393.*y - 6393.*x - 225000.*x.^2.*y.^3 - 225000.*x.^3.*y.^2 + 75000.*x.^2.*y.^4 + 100000.*x.^3.*y.^3 + 75000.*x.^4.*y.^2 + 44524.*x.*y - 122625.*x.*y.^2 - 122625.*x.^2.*y + 167000.*x.*y.^3 + 167000.*x.^3.*y - 112500.*x.*y.^4 - 112500.*x.^4.*y + 30000.*x.*y.^5 + 30000.*x.^5.*y + 22262.*x.^2 - 40875.*x.^3 + 41750.*x.^4 - 22500.*x.^5 + 5000.*x.^6 + 22262.*y.^2 - 40875.*y.^3 + 41750.*y.^4 - 22500.*y.^5 + 5000.*y.^6 + 756))/108;
            case 26
                z = -(x.*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(4881.*x + 4881.*y - 60000.*x.^2.*y.^2 + 25000.*x.^2.*y.^3 + 25000.*x.^3.*y.^2 - 25000.*x.*y + 47625.*x.*y.^2 + 47625.*x.^2.*y - 40000.*x.*y.^3 - 40000.*x.^3.*y + 12500.*x.*y.^4 + 12500.*x.^4.*y - 12500.*x.^2 + 15875.*x.^3 - 10000.*x.^4 + 2500.*x.^5 - 12500.*y.^2 + 15875.*y.^3 - 10000.*y.^4 + 2500.*y.^5 - 756))/9;
            case 27
                z = (25.*x.*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(3000.*x.^2.*y.^2 - 1207.*y - 1207.*x + 4310.*x.*y - 5100.*x.*y.^2 - 5100.*x.^2.*y + 2000.*x.*y.^3 + 2000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/108;
            case 28
                z = -(50.*x.*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/189;
            case 29
                z = (25.*x.*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9).*(3267.*x - 32830.*x.^2 + 169225.*x.^3 - 490000.*x.^4 + 805000.*x.^5 - 700000.*x.^6 + 250000.*x.^7 - 126))/504;
            case 30
                z = -(25.*x.*(x + y - 1).*(147655.*x.^2 - 13698.*x - 841050.*x.^3 + 2806125.*x.^4 - 5670000.*x.^5 + 6825000.*x.^6 - 4500000.*x.^7 + 1250000.*x.^8 + 504))/1134;
            case 31
                z = (125.*x.*y.*(16765350.*x.^2.*y.^2 - 64818.*y - 64818.*x - 36400000.*x.^2.*y.^3 - 36400000.*x.^3.*y.^2 + 43575000.*x.^2.*y.^4 + 58100000.*x.^3.*y.^3 + 43575000.*x.^4.*y.^2 - 27300000.*x.^2.*y.^5 - 45500000.*x.^3.*y.^4 - 45500000.*x.^4.*y.^3 - 27300000.*x.^5.*y.^2 + 7000000.*x.^2.*y.^6 + 14000000.*x.^3.*y.^5 + 17500000.*x.^4.*y.^4 + 14000000.*x.^5.*y.^3 + 7000000.*x.^6.*y.^2 + 790254.*x.*y - 4032210.*x.*y.^2 - 4032210.*x.^2.*y + 11176900.*x.*y.^3 + 11176900.*x.^3.*y - 18200000.*x.*y.^4 - 18200000.*x.^4.*y + 17430000.*x.*y.^5 + 17430000.*x.^5.*y - 9100000.*x.*y.^6 - 9100000.*x.^6.*y + 2000000.*x.*y.^7 + 2000000.*x.^7.*y + 395127.*x.^2 - 1344070.*x.^3 + 2794225.*x.^4 - 3640000.*x.^5 + 2905000.*x.^6 - 1300000.*x.^7 + 250000.*x.^8 + 395127.*y.^2 - 1344070.*y.^3 + 2794225.*y.^4 - 3640000.*y.^5 + 2905000.*y.^6 - 1300000.*y.^7 + 250000.*y.^8 + 4536))/126;
            case 32
                z = -(125.*x.*y.*(x + y - 1).*(3267.*x - 32830.*x.^2 + 169225.*x.^3 - 490000.*x.^4 + 805000.*x.^5 - 700000.*x.^6 + 250000.*x.^7 - 126))/126;
            case 33
                z = -(125.*x.*y.*(x + y - 1).*(3267.*y - 32830.*y.^2 + 169225.*y.^3 - 490000.*y.^4 + 805000.*y.^5 - 700000.*y.^6 + 250000.*y.^7 - 126))/126;
            case 34
                z = -(250.*x.*y.*(10.*y - 1).*(x + y - 1).*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18))/63;
            case 35
                z = -(250.*x.*y.*(x + y - 1).*(50.*y.^2 - 15.*y + 1).*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6))/27;
            case 36
                z = -(25.*x.*y.*(x + y - 1).*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6))/9;
            case 37
                z = -(25.*x.*y.*(x + y - 1).*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6))/9;
            case 38
                z = -(250.*x.*y.*(x + y - 1).*(50.*x.^2 - 15.*x + 1).*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6))/27;
            case 39
                z = -(250.*x.*y.*(10.*x - 1).*(x + y - 1).*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18))/63;
            case 40
                z = (250.*x.*y.*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9).*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18))/63;
            case 41
                z = -(250.*x.*y.*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/27;
            case 42
                z = (25.*x.*y.*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(3000.*x.^2.*y.^2 - 1207.*y - 1207.*x + 4310.*x.*y - 5100.*x.*y.^2 - 5100.*x.^2.*y + 2000.*x.*y.^3 + 2000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/9;
            case 43
                z = -(25.*x.*y.*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(4881.*x + 4881.*y - 60000.*x.^2.*y.^2 + 25000.*x.^2.*y.^3 + 25000.*x.^3.*y.^2 - 25000.*x.*y + 47625.*x.*y.^2 + 47625.*x.^2.*y - 40000.*x.*y.^3 - 40000.*x.^3.*y + 12500.*x.*y.^4 + 12500.*x.^4.*y - 12500.*x.^2 + 15875.*x.^3 - 10000.*x.^4 + 2500.*x.^5 - 12500.*y.^2 + 15875.*y.^3 - 10000.*y.^4 + 2500.*y.^5 - 756))/9;
            case 44
                z = (250.*x.*y.*(50.*y.^2 - 15.*y + 1).*(250500.*x.^2.*y.^2 - 6393.*y - 6393.*x - 225000.*x.^2.*y.^3 - 225000.*x.^3.*y.^2 + 75000.*x.^2.*y.^4 + 100000.*x.^3.*y.^3 + 75000.*x.^4.*y.^2 + 44524.*x.*y - 122625.*x.*y.^2 - 122625.*x.^2.*y + 167000.*x.*y.^3 + 167000.*x.^3.*y - 112500.*x.*y.^4 - 112500.*x.^4.*y + 30000.*x.*y.^5 + 30000.*x.^5.*y + 22262.*x.^2 - 40875.*x.^3 + 41750.*x.^4 - 22500.*x.^5 + 5000.*x.^6 + 22262.*y.^2 - 40875.*y.^3 + 41750.*y.^4 - 22500.*y.^5 + 5000.*y.^6 + 756))/27;
            case 45
                z = -(250.*x.*y.*(10.*y - 1).*(16566.*x + 16566.*y - 1727250.*x.^2.*y.^2 + 2537500.*x.^2.*y.^3 + 2537500.*x.^3.*y.^2 - 1837500.*x.^2.*y.^4 - 2450000.*x.^3.*y.^3 - 1837500.*x.^4.*y.^2 + 525000.*x.^2.*y.^5 + 875000.*x.^3.*y.^4 + 875000.*x.^4.*y.^3 + 525000.*x.^5.*y.^2 - 152978.*x.*y + 579180.*x.*y.^2 + 579180.*x.^2.*y - 1151500.*x.*y.^3 - 1151500.*x.^3.*y + 1268750.*x.*y.^4 + 1268750.*x.^4.*y - 735000.*x.*y.^5 - 735000.*x.^5.*y + 175000.*x.*y.^6 + 175000.*x.^6.*y - 76489.*x.^2 + 193060.*x.^3 - 287875.*x.^4 + 253750.*x.^5 - 122500.*x.^6 + 25000.*x.^7 - 76489.*y.^2 + 193060.*y.^3 - 287875.*y.^4 + 253750.*y.^5 - 122500.*y.^6 + 25000.*y.^7 - 1512))/63;
            case 46
                z = -(250.*x.*y.*(10.*x - 1).*(16566.*x + 16566.*y - 1727250.*x.^2.*y.^2 + 2537500.*x.^2.*y.^3 + 2537500.*x.^3.*y.^2 - 1837500.*x.^2.*y.^4 - 2450000.*x.^3.*y.^3 - 1837500.*x.^4.*y.^2 + 525000.*x.^2.*y.^5 + 875000.*x.^3.*y.^4 + 875000.*x.^4.*y.^3 + 525000.*x.^5.*y.^2 - 152978.*x.*y + 579180.*x.*y.^2 + 579180.*x.^2.*y - 1151500.*x.*y.^3 - 1151500.*x.^3.*y + 1268750.*x.*y.^4 + 1268750.*x.^4.*y - 735000.*x.*y.^5 - 735000.*x.^5.*y + 175000.*x.*y.^6 + 175000.*x.^6.*y - 76489.*x.^2 + 193060.*x.^3 - 287875.*x.^4 + 253750.*x.^5 - 122500.*x.^6 + 25000.*x.^7 - 76489.*y.^2 + 193060.*y.^3 - 287875.*y.^4 + 253750.*y.^5 - 122500.*y.^6 + 25000.*y.^7 - 1512))/63;
            case 47
                z = (250.*x.*y.*(50.*x.^2 - 15.*x + 1).*(250500.*x.^2.*y.^2 - 6393.*y - 6393.*x - 225000.*x.^2.*y.^3 - 225000.*x.^3.*y.^2 + 75000.*x.^2.*y.^4 + 100000.*x.^3.*y.^3 + 75000.*x.^4.*y.^2 + 44524.*x.*y - 122625.*x.*y.^2 - 122625.*x.^2.*y + 167000.*x.*y.^3 + 167000.*x.^3.*y - 112500.*x.*y.^4 - 112500.*x.^4.*y + 30000.*x.*y.^5 + 30000.*x.^5.*y + 22262.*x.^2 - 40875.*x.^3 + 41750.*x.^4 - 22500.*x.^5 + 5000.*x.^6 + 22262.*y.^2 - 40875.*y.^3 + 41750.*y.^4 - 22500.*y.^5 + 5000.*y.^6 + 756))/27;
            case 48
                z = -(25.*x.*y.*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(4881.*x + 4881.*y - 60000.*x.^2.*y.^2 + 25000.*x.^2.*y.^3 + 25000.*x.^3.*y.^2 - 25000.*x.*y + 47625.*x.*y.^2 + 47625.*x.^2.*y - 40000.*x.*y.^3 - 40000.*x.^3.*y + 12500.*x.*y.^4 + 12500.*x.^4.*y - 12500.*x.^2 + 15875.*x.^3 - 10000.*x.^4 + 2500.*x.^5 - 12500.*y.^2 + 15875.*y.^3 - 10000.*y.^4 + 2500.*y.^5 - 756))/9;
            case 49
                z = (25.*x.*y.*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(3000.*x.^2.*y.^2 - 1207.*y - 1207.*x + 4310.*x.*y - 5100.*x.*y.^2 - 5100.*x.^2.*y + 2000.*x.*y.^3 + 2000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/9;
            case 50
                z = -(250.*x.*y.*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/27;
            case 51
                z = (250.*x.*y.*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9).*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18))/63;
            case 52
                z = (125.*x.*y.*(10.*x - 1).*(10.*y - 1).*(250500.*x.^2.*y.^2 - 6393.*y - 6393.*x - 225000.*x.^2.*y.^3 - 225000.*x.^3.*y.^2 + 75000.*x.^2.*y.^4 + 100000.*x.^3.*y.^3 + 75000.*x.^4.*y.^2 + 44524.*x.*y - 122625.*x.*y.^2 - 122625.*x.^2.*y + 167000.*x.*y.^3 + 167000.*x.^3.*y - 112500.*x.*y.^4 - 112500.*x.^4.*y + 30000.*x.*y.^5 + 30000.*x.^5.*y + 22262.*x.^2 - 40875.*x.^3 + 41750.*x.^4 - 22500.*x.^5 + 5000.*x.^6 + 22262.*y.^2 - 40875.*y.^3 + 41750.*y.^4 - 22500.*y.^5 + 5000.*y.^6 + 756))/18;
            case 53
                z = (125.*x.*y.*(10.*y - 1).*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/18;
            case 54
                z = (125.*x.*y.*(10.*x - 1).*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/18;
            case 55
                z = (50.*x.*y.*(50.*y.^2 - 15.*y + 1).*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/9;
            case 56
                z = (125.*x.*y.*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/36;
            case 57
                z = (50.*x.*y.*(50.*x.^2 - 15.*x + 1).*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/9;
            case 58
                z = -(50.*x.*y.*(10.*x - 1).*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/9;
            case 59
                z = (125.*x.*y.*(10.*x - 1).*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(3000.*x.^2.*y.^2 - 1207.*y - 1207.*x + 4310.*x.*y - 5100.*x.*y.^2 - 5100.*x.^2.*y + 2000.*x.*y.^3 + 2000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/36;
            case 60
                z = -(50.*x.*y.*(10.*x - 1).*(50.*y.^2 - 15.*y + 1).*(4881.*x + 4881.*y - 60000.*x.^2.*y.^2 + 25000.*x.^2.*y.^3 + 25000.*x.^3.*y.^2 - 25000.*x.*y + 47625.*x.*y.^2 + 47625.*x.^2.*y - 40000.*x.*y.^3 - 40000.*x.^3.*y + 12500.*x.*y.^4 + 12500.*x.^4.*y - 12500.*x.^2 + 15875.*x.^3 - 10000.*x.^4 + 2500.*x.^5 - 12500.*y.^2 + 15875.*y.^3 - 10000.*y.^4 + 2500.*y.^5 - 756))/9;
            case 61
                z = -(50.*x.*y.*(10.*y - 1).*(50.*x.^2 - 15.*x + 1).*(4881.*x + 4881.*y - 60000.*x.^2.*y.^2 + 25000.*x.^2.*y.^3 + 25000.*x.^3.*y.^2 - 25000.*x.*y + 47625.*x.*y.^2 + 47625.*x.^2.*y - 40000.*x.*y.^3 - 40000.*x.^3.*y + 12500.*x.*y.^4 + 12500.*x.^4.*y - 12500.*x.^2 + 15875.*x.^3 - 10000.*x.^4 + 2500.*x.^5 - 12500.*y.^2 + 15875.*y.^3 - 10000.*y.^4 + 2500.*y.^5 - 756))/9;
            case 62
                z = (125.*x.*y.*(10.*y - 1).*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(3000.*x.^2.*y.^2 - 1207.*y - 1207.*x + 4310.*x.*y - 5100.*x.*y.^2 - 5100.*x.^2.*y + 2000.*x.*y.^3 + 2000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/36;
            case 63
                z = -(50.*x.*y.*(10.*y - 1).*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/9;
            case 64
                z = (250.*x.*y.*(50.*x.^2 - 15.*x + 1).*(50.*y.^2 - 15.*y + 1).*(3000.*x.^2.*y.^2 - 1207.*y - 1207.*x + 4310.*x.*y - 5100.*x.*y.^2 - 5100.*x.^2.*y + 2000.*x.*y.^3 + 2000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/27;
            case 65
                z = -(250.*x.*y.*(50.*y.^2 - 15.*y + 1).*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/27;
            case 66
                z = -(250.*x.*y.*(50.*x.^2 - 15.*x + 1).*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/27;
        end
end


end

function [gx,gy]=gradhphih(i,x,y,k)

switch k
    case 1
        switch i
            case 1
                gx = -1;
                gy = -1;
            case 2
                gx = 1;
                gy = 0;
            case 3
                gx = 0;
                gy = 1;
        end
    case 2
        switch i
            case 1
                gx = 4.*x + 4.*y - 3;
                gy = 4.*x + 4.*y - 3;
            case 2
                gx = 4.*x - 1;
                gy = 0;
            case 3
                gx = 0;
                gy = 4.*y - 1;
            case 4
                gx = 4.*y;
                gy = 4.*x;
            case 5
                gx = -4.*y;
                gy = 4 - 8.*y - 4.*x;
            case 6
                gx = 4 - 4.*y - 8.*x;
                gy = -4.*x;
        end
    case 3
        switch i
            case 1
                gx = 18.*x + 18.*y - 27.*x.*y - (27.*x.^2)/2 - (27.*y.^2)/2 - 11/2;
                gy = 18.*x + 18.*y - 27.*x.*y - (27.*x.^2)/2 - (27.*y.^2)/2 - 11/2;
            case 2
                gx = (27.*x.^2)/2 - 9.*x + 1;
                gy = 0;
            case 3
                gx = 0;
                gy = (27.*y.^2)/2 - 9.*y + 1;
            case 4
                gx = (9.*y.*(6.*x - 1))/2;
                gy = (9.*x.*(3.*x - 1))/2;
            case 5
                gx = (9.*y.*(3.*y - 1))/2;
                gy = (9.*x.*(6.*y - 1))/2;
            case 6
                gx = -(9.*y.*(3.*y - 1))/2;
                gy = (9.*x)/2 + 36.*y - 27.*x.*y - (81.*y.^2)/2 - 9/2;
            case 7
                gx = (9.*y.*(6.*x + 6.*y - 5))/2;
                gy = 54.*x.*y - 45.*y - (45.*x)/2 + (27.*x.^2)/2 + (81.*y.^2)/2 + 9;
            case 8
                gx = 54.*x.*y - (45.*y)/2 - 45.*x + (81.*x.^2)/2 + (27.*y.^2)/2 + 9;
                gy = (9.*x.*(6.*x + 6.*y - 5))/2;
            case 9
                gx = 36.*x + (9.*y)/2 - 27.*x.*y - (81.*x.^2)/2 - 9/2;
                gy = -(9.*x.*(3.*x - 1))/2;
            case 10
                gx = -27.*y.*(2.*x + y - 1);
                gy = -27.*x.*(x + 2.*y - 1);
        end
    case 4
        switch i
            case 1
                gx = (140.*x)/3 + (140.*y)/3 - 160.*x.*y + 128.*x.*y.^2 + 128.*x.^2.*y - 80.*x.^2 + (128.*x.^3)/3 - 80.*y.^2 + (128.*y.^3)/3 - 25/3;
                gy = (140.*x)/3 + (140.*y)/3 - 160.*x.*y + 128.*x.*y.^2 + 128.*x.^2.*y - 80.*x.^2 + (128.*x.^3)/3 - 80.*y.^2 + (128.*y.^3)/3 - 25/3;
            case 2
                gx = (44.*x)/3 - 48.*x.^2 + (128.*x.^3)/3 - 1;
                gy = 0;
            case 3
                gx = 0;
                gy = (44.*y)/3 - 48.*y.^2 + (128.*y.^3)/3 - 1;
            case 4
                gx = (16.*y.*(24.*x.^2 - 12.*x + 1))/3;
                gy = (16.*x.*(8.*x.^2 - 6.*x + 1))/3;
            case 5
                gx = 4.*y.*(8.*x - 1).*(4.*y - 1);
                gy = 4.*x.*(4.*x - 1).*(8.*y - 1);
            case 6
                gx = (16.*y.*(8.*y.^2 - 6.*y + 1))/3;
                gy = (16.*x.*(24.*y.^2 - 12.*y + 1))/3;
            case 7
                gx = -(16.*y.*(8.*y.^2 - 6.*y + 1))/3;
                gy = 64.*x.*y - (224.*y)/3 - (16.*x)/3 - 128.*x.*y.^2 + 224.*y.^2 - (512.*y.^3)/3 + 16/3;
            case 8
                gx = 4.*y.*(4.*y - 1).*(8.*x + 8.*y - 7);
                gy = 28.*x + 152.*y - 288.*x.*y + 384.*x.*y.^2 + 128.*x.^2.*y - 16.*x.^2 - 384.*y.^2 + 256.*y.^3 - 12;
            case 9
                gx = -(16.*y.*(48.*x.*y - 36.*y - 36.*x + 24.*x.^2 + 24.*y.^2 + 13))/3;
                gy = 384.*x.*y - (416.*y)/3 - (208.*x)/3 - 384.*x.*y.^2 - 256.*x.^2.*y + 96.*x.^2 - (128.*x.^3)/3 + 288.*y.^2 - (512.*y.^3)/3 + 16;
            case 10
                gx = 384.*x.*y - (208.*y)/3 - (416.*x)/3 - 256.*x.*y.^2 - 384.*x.^2.*y + 288.*x.^2 - (512.*x.^3)/3 + 96.*y.^2 - (128.*y.^3)/3 + 16;
                gy = -(16.*x.*(48.*x.*y - 36.*y - 36.*x + 24.*x.^2 + 24.*y.^2 + 13))/3;
            case 11
                gx = 152.*x + 28.*y - 288.*x.*y + 128.*x.*y.^2 + 384.*x.^2.*y - 384.*x.^2 + 256.*x.^3 - 16.*y.^2 - 12;
                gy = 4.*x.*(4.*x - 1).*(8.*x + 8.*y - 7);
            case 12
                gx = 64.*x.*y - (16.*y)/3 - (224.*x)/3 - 128.*x.^2.*y + 224.*x.^2 - (512.*x.^3)/3 + 16/3;
                gy = -(16.*x.*(8.*x.^2 - 6.*x + 1))/3;
            case 13
                gx = 32.*y.*(16.*x.*y - 7.*y - 14.*x + 12.*x.^2 + 4.*y.^2 + 3);
                gy = 32.*x.*(16.*x.*y - 14.*y - 7.*x + 4.*x.^2 + 12.*y.^2 + 3);
            case 14
                gx = -32.*y.*(8.*x.*y - y - 10.*x + 12.*x.^2 + 1);
                gy = -32.*x.*(4.*x - 1).*(x + 2.*y - 1);
            case 15
                gx = -32.*y.*(4.*y - 1).*(2.*x + y - 1);
                gy = -32.*x.*(8.*x.*y - 10.*y - x + 12.*y.^2 + 1);
        end
    case 5
        switch i
            case 1
                gx = (375.*x)/4 + (375.*y)/4 - (3125.*x.^2.*y.^2)/4 - (2125.*x.*y)/4 + (1875.*x.*y.^2)/2 + (1875.*x.^2.*y)/2 - (3125.*x.*y.^3)/6 - (3125.*x.^3.*y)/6 - (2125.*x.^2)/8 + (625.*x.^3)/2 - (3125.*x.^4)/24 - (2125.*y.^2)/8 + (625.*y.^3)/2 - (3125.*y.^4)/24 - 137/12;
                gy = (375.*x)/4 + (375.*y)/4 - (3125.*x.^2.*y.^2)/4 - (2125.*x.*y)/4 + (1875.*x.*y.^2)/2 + (1875.*x.^2.*y)/2 - (3125.*x.*y.^3)/6 - (3125.*x.^3.*y)/6 - (2125.*x.^2)/8 + (625.*x.^3)/2 - (3125.*x.^4)/24 - (2125.*y.^2)/8 + (625.*y.^3)/2 - (3125.*y.^4)/24 - 137/12;
            case 2
                gx = (875.*x.^2)/8 - (125.*x)/6 - (625.*x.^3)/3 + (3125.*x.^4)/24 + 1;
                gy = 0;
            case 3
                gx = 0;
                gy = (875.*y.^2)/8 - (125.*y)/6 - (625.*y.^3)/3 + (3125.*y.^4)/24 + 1;
            case 4
                gx = (25.*y.*(55.*x - 225.*x.^2 + 250.*x.^3 - 3))/12;
                gy = (25.*x.*(55.*x - 150.*x.^2 + 125.*x.^3 - 6))/24;
            case 5
                gx = (25.*y.*(5.*y - 1).*(75.*x.^2 - 30.*x + 2))/12;
                gy = (25.*x.*(10.*y - 1).*(25.*x.^2 - 15.*x + 2))/12;
            case 6
                gx = (25.*y.*(10.*x - 1).*(25.*y.^2 - 15.*y + 2))/12;
                gy = (25.*x.*(5.*x - 1).*(75.*y.^2 - 30.*y + 2))/12;
            case 7
                gx = (25.*y.*(55.*y - 150.*y.^2 + 125.*y.^3 - 6))/24;
                gy = (25.*x.*(55.*y - 225.*y.^2 + 250.*y.^3 - 3))/12;
            case 8
                gx = -(25.*y.*(55.*y - 150.*y.^2 + 125.*y.^3 - 6))/24;
                gy = (25.*x)/4 + (1525.*y)/12 - (1375.*x.*y)/12 + (1875.*x.*y.^2)/4 - (3125.*x.*y.^3)/6 - (5125.*y.^2)/8 + (6875.*y.^3)/6 - (15625.*y.^4)/24 - 25/4;
            case 9
                gx = (25.*y.*(10.*x + 10.*y - 9).*(25.*y.^2 - 15.*y + 2))/12;
                gy = (3125.*x.^2.*y.^2)/4 - 325.*y - (75.*x)/2 + (3875.*x.*y)/6 - (9375.*x.*y.^2)/4 - (625.*x.^2.*y)/2 + (6250.*x.*y.^3)/3 + (125.*x.^2)/6 + (6125.*y.^2)/4 - 2500.*y.^3 + (15625.*y.^4)/12 + 50/3;
            case 10
                gx = -(25.*y.*(5.*y - 1).*(150.*x.*y - 120.*y - 120.*x + 75.*x.^2 + 75.*y.^2 + 47))/12;
                gy = (1175.*x)/12 + (2675.*y)/6 - (9375.*x.^2.*y.^2)/4 - (8875.*x.*y)/6 + (16875.*x.*y.^2)/4 + (3125.*x.^2.*y)/2 - 3125.*x.*y.^3 - (3125.*x.^3.*y)/6 - 125.*x.^2 + (625.*x.^3)/12 - (7375.*y.^2)/4 + (8125.*y.^3)/3 - (15625.*y.^4)/12 - 25;
            case 11
                gx = (25.*y.*(710.*x + 710.*y - 2100.*x.*y + 1500.*x.*y.^2 + 1500.*x.^2.*y - 1050.*x.^2 + 500.*x.^3 - 1050.*y.^2 + 500.*y.^3 - 154))/24;
                gy = (9375.*x.^2.*y.^2)/4 - (1925.*y)/6 - (1925.*x)/12 + (8875.*x.*y)/6 - (13125.*x.*y.^2)/4 - (4375.*x.^2.*y)/2 + (6250.*x.*y.^3)/3 + (3125.*x.^3.*y)/3 + (8875.*x.^2)/24 - (4375.*x.^3)/12 + (3125.*x.^4)/24 + (8875.*y.^2)/8 - (4375.*y.^3)/3 + (15625.*y.^4)/24 + 25;
            case 12
                gx = (9375.*x.^2.*y.^2)/4 - (1925.*y)/12 - (1925.*x)/6 + (8875.*x.*y)/6 - (4375.*x.*y.^2)/2 - (13125.*x.^2.*y)/4 + (3125.*x.*y.^3)/3 + (6250.*x.^3.*y)/3 + (8875.*x.^2)/8 - (4375.*x.^3)/3 + (15625.*x.^4)/24 + (8875.*y.^2)/24 - (4375.*y.^3)/12 + (3125.*y.^4)/24 + 25;
                gy = (25.*x.*(710.*x + 710.*y - 2100.*x.*y + 1500.*x.*y.^2 + 1500.*x.^2.*y - 1050.*x.^2 + 500.*x.^3 - 1050.*y.^2 + 500.*y.^3 - 154))/24;
            case 13
                gx = (2675.*x)/6 + (1175.*y)/12 - (9375.*x.^2.*y.^2)/4 - (8875.*x.*y)/6 + (3125.*x.*y.^2)/2 + (16875.*x.^2.*y)/4 - (3125.*x.*y.^3)/6 - 3125.*x.^3.*y - (7375.*x.^2)/4 + (8125.*x.^3)/3 - (15625.*x.^4)/12 - 125.*y.^2 + (625.*y.^3)/12 - 25;
                gy = -(25.*x.*(5.*x - 1).*(150.*x.*y - 120.*y - 120.*x + 75.*x.^2 + 75.*y.^2 + 47))/12;
            case 14
                gx = (3125.*x.^2.*y.^2)/4 - (75.*y)/2 - 325.*x + (3875.*x.*y)/6 - (625.*x.*y.^2)/2 - (9375.*x.^2.*y)/4 + (6250.*x.^3.*y)/3 + (6125.*x.^2)/4 - 2500.*x.^3 + (15625.*x.^4)/12 + (125.*y.^2)/6 + 50/3;
                gy = (25.*x.*(10.*x + 10.*y - 9).*(25.*x.^2 - 15.*x + 2))/12;
            case 15
                gx = (1525.*x)/12 + (25.*y)/4 - (1375.*x.*y)/12 + (1875.*x.^2.*y)/4 - (3125.*x.^3.*y)/6 - (5125.*x.^2)/8 + (6875.*x.^3)/6 - (15625.*x.^4)/24 - 25/4;
                gy = -(25.*x.*(55.*x - 150.*x.^2 + 125.*x.^3 - 6))/24;
            case 16
                gx = -(125.*y.*(94.*x + 47.*y - 240.*x.*y + 150.*x.*y.^2 + 225.*x.^2.*y - 180.*x.^2 + 100.*x.^3 - 60.*y.^2 + 25.*y.^3 - 12))/6;
                gy = -(125.*x.*(47.*x + 94.*y - 240.*x.*y + 225.*x.*y.^2 + 150.*x.^2.*y - 60.*x.^2 + 25.*x.^3 - 180.*y.^2 + 100.*y.^3 - 12))/6;
            case 17
                gx = -(125.*y.*(34.*x + 2.*y - 30.*x.*y + 75.*x.^2.*y - 120.*x.^2 + 100.*x.^3 - 2))/6;
                gy = -(125.*x.*(25.*x.^2 - 15.*x + 2).*(x + 2.*y - 1))/6;
            case 18
                gx = -(125.*y.*(25.*y.^2 - 15.*y + 2).*(2.*x + y - 1))/6;
                gy = -(125.*x.*(2.*x + 34.*y - 30.*x.*y + 75.*x.*y.^2 - 120.*y.^2 + 100.*y.^3 - 2))/6;
            case 19
                gx = -(125.*y.*(5.*y - 1).*(10.*x.*y - y - 12.*x + 15.*x.^2 + 1))/4;
                gy = -(125.*x.*(5.*x - 1).*(10.*x.*y - 12.*y - x + 15.*y.^2 + 1))/4;
            case 20
                gx = (125.*y.*(5.*y - 1).*(20.*x.*y - 9.*y - 18.*x + 15.*x.^2 + 5.*y.^2 + 4))/4;
                gy = (125.*x.*(9.*x + 58.*y - 110.*x.*y + 150.*x.*y.^2 + 50.*x.^2.*y - 5.*x.^2 - 150.*y.^2 + 100.*y.^3 - 4))/4;
            case 21
                gx = (125.*y.*(58.*x + 9.*y - 110.*x.*y + 50.*x.*y.^2 + 150.*x.^2.*y - 150.*x.^2 + 100.*x.^3 - 5.*y.^2 - 4))/4;
                gy = (125.*x.*(5.*x - 1).*(20.*x.*y - 18.*y - 9.*x + 5.*x.^2 + 15.*y.^2 + 4))/4;
        end
    case 6
        switch i
            case 1
                gx = (812.*x)/5 + (812.*y)/5 - 6804.*x.^2.*y.^2 + 3888.*x.^2.*y.^3 + 3888.*x.^3.*y.^2 - 1323.*x.*y + 3780.*x.*y.^2 + 3780.*x.^2.*y - 4536.*x.*y.^3 - 4536.*x.^3.*y + 1944.*x.*y.^4 + 1944.*x.^4.*y - (1323.*x.^2)/2 + 1260.*x.^3 - 1134.*x.^4 + (1944.*x.^5)/5 - (1323.*y.^2)/2 + 1260.*y.^3 - 1134.*y.^4 + (1944.*y.^5)/5 - 147/10;
                gy = (812.*x)/5 + (812.*y)/5 - 6804.*x.^2.*y.^2 + 3888.*x.^2.*y.^3 + 3888.*x.^3.*y.^2 - 1323.*x.*y + 3780.*x.*y.^2 + 3780.*x.^2.*y - 4536.*x.*y.^3 - 4536.*x.^3.*y + 1944.*x.*y.^4 + 1944.*x.^4.*y - (1323.*x.^2)/2 + 1260.*x.^3 - 1134.*x.^4 + (1944.*x.^5)/5 - (1323.*y.^2)/2 + 1260.*y.^3 - 1134.*y.^4 + (1944.*y.^5)/5 - 147/10;
            case 2
                gx = (137.*x)/5 - (405.*x.^2)/2 + 612.*x.^3 - 810.*x.^4 + (1944.*x.^5)/5 - 1;
                gy = 0;
            case 3
                gx = 0;
                gy = (137.*y)/5 - (405.*y.^2)/2 + 612.*y.^3 - 810.*y.^4 + (1944.*y.^5)/5 - 1;
            case 4
                gx = (18.*y.*(315.*x.^2 - 50.*x - 720.*x.^3 + 540.*x.^4 + 2))/5;
                gy = (18.*x.*(105.*x.^2 - 25.*x - 180.*x.^3 + 108.*x.^4 + 2))/5;
            case 5
                gx = (9.*y.*(6.*y - 1).*(22.*x - 108.*x.^2 + 144.*x.^3 - 1))/2;
                gy = (9.*x.*(12.*y - 1).*(11.*x - 36.*x.^2 + 36.*x.^3 - 1))/2;
            case 6
                gx = 4.*y.*(54.*x.^2 - 18.*x + 1).*(18.*y.^2 - 9.*y + 1);
                gy = 4.*x.*(18.*x.^2 - 9.*x + 1).*(54.*y.^2 - 18.*y + 1);
            case 7
                gx = (9.*y.*(12.*x - 1).*(11.*y - 36.*y.^2 + 36.*y.^3 - 1))/2;
                gy = (9.*x.*(6.*x - 1).*(22.*y - 108.*y.^2 + 144.*y.^3 - 1))/2;
            case 8
                gx = (18.*y.*(105.*y.^2 - 25.*y - 180.*y.^3 + 108.*y.^4 + 2))/5;
                gy = (18.*x.*(315.*y.^2 - 50.*y - 720.*y.^3 + 540.*y.^4 + 2))/5;
            case 9
                gx = -(18.*y.*(105.*y.^2 - 25.*y - 180.*y.^3 + 108.*y.^4 + 2))/5;
                gy = 180.*x.*y - (972.*y)/5 - (36.*x)/5 - 1134.*x.*y.^2 + 2592.*x.*y.^3 - 1944.*x.*y.^4 + 1404.*y.^2 - 4104.*y.^3 + 5184.*y.^4 - (11664.*y.^5)/5 + 36/5;
            case 10
                gx = (9.*y.*(12.*x + 12.*y - 11).*(11.*y - 36.*y.^2 + 36.*y.^3 - 1))/2;
                gy = (9.*(11.*y - 36.*y.^2 + 36.*y.^3 - 1).*(12.*x.*y - 11.*y - 11.*x + 6.*x.^2 + 6.*y.^2 + 5))/2 + (9.*y.*(108.*y.^2 - 72.*y + 11).*(12.*x.*y - 11.*y - 11.*x + 6.*x.^2 + 6.*y.^2 + 5))/2 + (9.*y.*(12.*x + 12.*y - 11).*(11.*y - 36.*y.^2 + 36.*y.^3 - 1))/2;
            case 11
                gx = -4.*y.*(18.*y.^2 - 9.*y + 1).*(108.*x.*y - 90.*y - 90.*x + 54.*x.^2 + 54.*y.^2 + 37);
                gy = 15552.*x.^2.*y.^2 - 1016.*y - 148.*x - 15552.*x.^2.*y.^3 - 3888.*x.^3.*y.^2 + 3384.*x.*y - 18360.*x.*y.^2 - 3672.*x.^2.*y + 33696.*x.*y.^3 + 1296.*x.^3.*y - 19440.*x.*y.^4 + 180.*x.^2 - 72.*x.^3 + 6696.*y.^2 - 17424.*y.^3 + 19440.*y.^4 - 7776.*y.^5 + 40;
            case 12
                gx = (9.*y.*(6.*y - 1).*(238.*x + 238.*y - 648.*x.*y + 432.*x.*y.^2 + 432.*x.^2.*y - 324.*x.^2 + 144.*x.^3 - 324.*y.^2 + 144.*y.^3 - 57))/2;
                gy = (513.*x)/2 + 1053.*y - 29160.*x.^2.*y.^2 + 23328.*x.^2.*y.^3 + 11664.*x.^3.*y.^2 - 5220.*x.*y + 23652.*x.*y.^2 + 9342.*x.^2.*y - 37584.*x.*y.^3 - 7128.*x.^3.*y + 19440.*x.*y.^4 + 1944.*x.^4.*y - (1071.*x.^2)/2 + 486.*x.^3 - 162.*x.^4 - (12447.*y.^2)/2 + 14796.*y.^3 - 15390.*y.^4 + 5832.*y.^5 - 45;
            case 13
                gx = -(18.*y.*(3240.*x.^2.*y.^2 - 580.*y - 580.*x + 2790.*x.*y - 4320.*x.*y.^2 - 4320.*x.^2.*y + 2160.*x.*y.^3 + 2160.*x.^3.*y + 1395.*x.^2 - 1440.*x.^3 + 540.*x.^4 + 1395.*y.^2 - 1440.*y.^3 + 540.*y.^4 + 87))/5;
                gy = 23328.*x.^2.*y.^2 - (3132.*y)/5 - (1566.*x)/5 - 15552.*x.^2.*y.^3 - 11664.*x.^3.*y.^2 + 4176.*x.*y - 15066.*x.*y.^2 - 10044.*x.^2.*y + 20736.*x.*y.^3 + 10368.*x.^3.*y - 9720.*x.*y.^4 - 3888.*x.^4.*y + 1044.*x.^2 - 1674.*x.^3 + 1296.*x.^4 - (1944.*x.^5)/5 + 3132.*y.^2 - 6696.*y.^3 + 6480.*y.^4 - (11664.*y.^5)/5 + 36;
            case 14
                gx = 23328.*x.^2.*y.^2 - (1566.*y)/5 - (3132.*x)/5 - 11664.*x.^2.*y.^3 - 15552.*x.^3.*y.^2 + 4176.*x.*y - 10044.*x.*y.^2 - 15066.*x.^2.*y + 10368.*x.*y.^3 + 20736.*x.^3.*y - 3888.*x.*y.^4 - 9720.*x.^4.*y + 3132.*x.^2 - 6696.*x.^3 + 6480.*x.^4 - (11664.*x.^5)/5 + 1044.*y.^2 - 1674.*y.^3 + 1296.*y.^4 - (1944.*y.^5)/5 + 36;
                gy = -(18.*x.*(3240.*x.^2.*y.^2 - 580.*y - 580.*x + 2790.*x.*y - 4320.*x.*y.^2 - 4320.*x.^2.*y + 2160.*x.*y.^3 + 2160.*x.^3.*y + 1395.*x.^2 - 1440.*x.^3 + 540.*x.^4 + 1395.*y.^2 - 1440.*y.^3 + 540.*y.^4 + 87))/5;
            case 15
                gx = 1053.*x + (513.*y)/2 - 29160.*x.^2.*y.^2 + 11664.*x.^2.*y.^3 + 23328.*x.^3.*y.^2 - 5220.*x.*y + 9342.*x.*y.^2 + 23652.*x.^2.*y - 7128.*x.*y.^3 - 37584.*x.^3.*y + 1944.*x.*y.^4 + 19440.*x.^4.*y - (12447.*x.^2)/2 + 14796.*x.^3 - 15390.*x.^4 + 5832.*x.^5 - (1071.*y.^2)/2 + 486.*y.^3 - 162.*y.^4 - 45;
                gy = (9.*x.*(6.*x - 1).*(238.*x + 238.*y - 648.*x.*y + 432.*x.*y.^2 + 432.*x.^2.*y - 324.*x.^2 + 144.*x.^3 - 324.*y.^2 + 144.*y.^3 - 57))/2;
            case 16
                gx = 15552.*x.^2.*y.^2 - 148.*y - 1016.*x - 3888.*x.^2.*y.^3 - 15552.*x.^3.*y.^2 + 3384.*x.*y - 3672.*x.*y.^2 - 18360.*x.^2.*y + 1296.*x.*y.^3 + 33696.*x.^3.*y - 19440.*x.^4.*y + 6696.*x.^2 - 17424.*x.^3 + 19440.*x.^4 - 7776.*x.^5 + 180.*y.^2 - 72.*y.^3 + 40;
                gy = -4.*x.*(18.*x.^2 - 9.*x + 1).*(108.*x.*y - 90.*y - 90.*x + 54.*x.^2 + 54.*y.^2 + 37);
            case 17
                gx = (9.*(11.*x - 36.*x.^2 + 36.*x.^3 - 1).*(12.*x.*y - 11.*y - 11.*x + 6.*x.^2 + 6.*y.^2 + 5))/2 + (9.*x.*(108.*x.^2 - 72.*x + 11).*(12.*x.*y - 11.*y - 11.*x + 6.*x.^2 + 6.*y.^2 + 5))/2 + (9.*x.*(12.*x + 12.*y - 11).*(11.*x - 36.*x.^2 + 36.*x.^3 - 1))/2;
                gy = (9.*x.*(12.*x + 12.*y - 11).*(11.*x - 36.*x.^2 + 36.*x.^3 - 1))/2;
            case 18
                gx = 180.*x.*y - (36.*y)/5 - (972.*x)/5 - 1134.*x.^2.*y + 2592.*x.^3.*y - 1944.*x.^4.*y + 1404.*x.^2 - 4104.*x.^3 + 5184.*x.^4 - (11664.*x.^5)/5 + 36/5;
                gy = -(18.*x.*(105.*x.^2 - 25.*x - 180.*x.^3 + 108.*x.^4 + 2))/5;
            case 19
                gx = 54.*y.*(648.*x.^2.*y.^2 - 57.*y - 114.*x + 476.*x.*y - 648.*x.*y.^2 - 972.*x.^2.*y + 288.*x.*y.^3 + 576.*x.^3.*y + 357.*x.^2 - 432.*x.^3 + 180.*x.^4 + 119.*y.^2 - 108.*y.^3 + 36.*y.^4 + 10);
                gy = 54.*x.*(648.*x.^2.*y.^2 - 114.*y - 57.*x + 476.*x.*y - 972.*x.*y.^2 - 648.*x.^2.*y + 576.*x.*y.^3 + 288.*x.^3.*y + 119.*x.^2 - 108.*x.^3 + 36.*x.^4 + 357.*y.^2 - 432.*y.^3 + 180.*y.^4 + 10);
            case 20
                gx = -54.*y.*(22.*x.*y - y - 24.*x - 108.*x.^2.*y + 144.*x.^3.*y + 141.*x.^2 - 288.*x.^3 + 180.*x.^4 + 1);
                gy = -54.*x.*(x + 2.*y - 1).*(11.*x - 36.*x.^2 + 36.*x.^3 - 1);
            case 21
                gx = -54.*y.*(2.*x + y - 1).*(11.*y - 36.*y.^2 + 36.*y.^3 - 1);
                gy = -54.*x.*(22.*x.*y - 24.*y - x - 108.*x.*y.^2 + 144.*x.*y.^3 + 141.*y.^2 - 288.*y.^3 + 180.*y.^4 + 1);
            case 22
                gx = -36.*y.*(6.*y - 1).*(20.*x + y - 18.*x.*y + 54.*x.^2.*y - 81.*x.^2 + 72.*x.^3 - 1);
                gy = -36.*x.*(18.*x.^2 - 9.*x + 1).*(12.*x.*y - 14.*y - x + 18.*y.^2 + 1);
            case 23
                gx = -36.*y.*(18.*y.^2 - 9.*y + 1).*(12.*x.*y - y - 14.*x + 18.*x.^2 + 1);
                gy = -36.*x.*(6.*x - 1).*(x + 20.*y - 18.*x.*y + 54.*x.*y.^2 - 81.*y.^2 + 72.*y.^3 - 1);
            case 24
                gx = 36.*y.*(18.*y.^2 - 9.*y + 1).*(24.*x.*y - 11.*y - 22.*x + 18.*x.^2 + 6.*y.^2 + 5);
                gy = 36.*x.*(324.*x.^2.*y.^2 - 112.*y - 11.*x + 222.*x.*y - 918.*x.*y.^2 - 108.*x.^2.*y + 864.*x.*y.^3 + 6.*x.^2 + 585.*y.^2 - 1008.*y.^3 + 540.*y.^4 + 5);
            case 25
                gx = -36.*y.*(6.*y - 1).*(74.*x + 37.*y - 180.*x.*y + 108.*x.*y.^2 + 162.*x.^2.*y - 135.*x.^2 + 72.*x.^3 - 45.*y.^2 + 18.*y.^3 - 10);
                gy = -36.*x.*(972.*x.^2.*y.^2 - 194.*y - 37.*x + 624.*x.*y - 1782.*x.*y.^2 - 648.*x.^2.*y + 1296.*x.*y.^3 + 216.*x.^3.*y + 45.*x.^2 - 18.*x.^3 + 801.*y.^2 - 1152.*y.^3 + 540.*y.^4 + 10);
            case 26
                gx = -36.*y.*(972.*x.^2.*y.^2 - 37.*y - 194.*x + 624.*x.*y - 648.*x.*y.^2 - 1782.*x.^2.*y + 216.*x.*y.^3 + 1296.*x.^3.*y + 801.*x.^2 - 1152.*x.^3 + 540.*x.^4 + 45.*y.^2 - 18.*y.^3 + 10);
                gy = -36.*x.*(6.*x - 1).*(37.*x + 74.*y - 180.*x.*y + 162.*x.*y.^2 + 108.*x.^2.*y - 45.*x.^2 + 18.*x.^3 - 135.*y.^2 + 72.*y.^3 - 10);
            case 27
                gx = 36.*y.*(324.*x.^2.*y.^2 - 11.*y - 112.*x + 222.*x.*y - 108.*x.*y.^2 - 918.*x.^2.*y + 864.*x.^3.*y + 585.*x.^2 - 1008.*x.^3 + 540.*x.^4 + 6.*y.^2 + 5);
                gy = 36.*x.*(18.*x.^2 - 9.*x + 1).*(24.*x.*y - 22.*y - 11.*x + 6.*x.^2 + 18.*y.^2 + 5);
            case 28
                gx = 27.*y.*(6.*y - 1).*(82.*x + 11.*y - 156.*x.*y + 72.*x.*y.^2 + 216.*x.^2.*y - 216.*x.^2 + 144.*x.^3 - 6.*y.^2 - 5);
                gy = 27.*x.*(6.*x - 1).*(11.*x + 82.*y - 156.*x.*y + 216.*x.*y.^2 + 72.*x.^2.*y - 6.*x.^2 - 216.*y.^2 + 144.*y.^3 - 5);
        end
    case 7
        switch i
            case 1
                gx = (22981.*x)/90 + (22981.*y)/90 - (386561.*x.^2.*y.^2)/12 + (117649.*x.^2.*y.^3)/3 + (117649.*x.^3.*y.^2)/3 - (823543.*x.^2.*y.^4)/48 - (823543.*x.^3.*y.^3)/36 - (823543.*x.^4.*y.^2)/48 - (331681.*x.*y)/120 + (33614.*x.*y.^2)/3 + (33614.*x.^2.*y)/3 - (386561.*x.*y.^3)/18 - (386561.*x.^3.*y)/18 + (117649.*x.*y.^4)/6 + (117649.*x.^4.*y)/6 - (823543.*x.*y.^5)/120 - (823543.*x.^5.*y)/120 - (331681.*x.^2)/240 + (33614.*x.^3)/9 - (386561.*x.^4)/72 + (117649.*x.^5)/30 - (823543.*x.^6)/720 - (331681.*y.^2)/240 + (33614.*y.^3)/9 - (386561.*y.^4)/72 + (117649.*y.^5)/30 - (823543.*y.^6)/720 - 363/20;
                gy = (22981.*x)/90 + (22981.*y)/90 - (386561.*x.^2.*y.^2)/12 + (117649.*x.^2.*y.^3)/3 + (117649.*x.^3.*y.^2)/3 - (823543.*x.^2.*y.^4)/48 - (823543.*x.^3.*y.^3)/36 - (823543.*x.^4.*y.^2)/48 - (331681.*x.*y)/120 + (33614.*x.*y.^2)/3 + (33614.*x.^2.*y)/3 - (386561.*x.*y.^3)/18 - (386561.*x.^3.*y)/18 + (117649.*x.*y.^4)/6 + (117649.*x.^4.*y)/6 - (823543.*x.*y.^5)/120 - (823543.*x.^5.*y)/120 - (331681.*x.^2)/240 + (33614.*x.^3)/9 - (386561.*x.^4)/72 + (117649.*x.^5)/30 - (823543.*x.^6)/720 - (331681.*y.^2)/240 + (33614.*y.^3)/9 - (386561.*y.^4)/72 + (117649.*y.^5)/30 - (823543.*y.^6)/720 - 363/20;
            case 2
                gx = (9947.*x.^2)/30 - (343.*x)/10 - (16807.*x.^3)/12 + (420175.*x.^4)/144 - (117649.*x.^5)/40 + (823543.*x.^6)/720 + 1;
                gy = 0;
            case 3
                gx = 0;
                gy = (9947.*y.^2)/30 - (343.*y)/10 - (16807.*y.^3)/12 + (420175.*y.^4)/144 - (117649.*y.^5)/40 + (823543.*y.^6)/720 + 1;
            case 4
                gx = (49.*y.*(3836.*x - 33075.*x.^2 + 116620.*x.^3 - 180075.*x.^4 + 100842.*x.^5 - 120))/720;
                gy = (49.*x.*(1918.*x - 11025.*x.^2 + 29155.*x.^3 - 36015.*x.^4 + 16807.*x.^5 - 120))/720;
            case 5
                gx = (49.*y.*(7.*y - 1).*(5145.*x.^2 - 700.*x - 13720.*x.^3 + 12005.*x.^4 + 24))/240;
                gy = (49.*x.*(14.*y - 1).*(1715.*x.^2 - 350.*x - 3430.*x.^3 + 2401.*x.^4 + 24))/240;
            case 6
                gx = (49.*y.*(49.*y.^2 - 21.*y + 2).*(77.*x - 441.*x.^2 + 686.*x.^3 - 3))/72;
                gy = (49.*x.*(147.*y.^2 - 42.*y + 2).*(77.*x - 294.*x.^2 + 343.*x.^3 - 6))/144;
            case 7
                gx = (49.*y.*(147.*x.^2 - 42.*x + 2).*(77.*y - 294.*y.^2 + 343.*y.^3 - 6))/144;
                gy = (49.*x.*(49.*x.^2 - 21.*x + 2).*(77.*y - 441.*y.^2 + 686.*y.^3 - 3))/72;
            case 8
                gx = (49.*y.*(14.*x - 1).*(1715.*y.^2 - 350.*y - 3430.*y.^3 + 2401.*y.^4 + 24))/240;
                gy = (49.*x.*(7.*x - 1).*(5145.*y.^2 - 700.*y - 13720.*y.^3 + 12005.*y.^4 + 24))/240;
            case 9
                gx = (49.*y.*(1918.*y - 11025.*y.^2 + 29155.*y.^3 - 36015.*y.^4 + 16807.*y.^5 - 120))/720;
                gy = (49.*x.*(3836.*y - 33075.*y.^2 + 116620.*y.^3 - 180075.*y.^4 + 100842.*y.^5 - 120))/720;
            case 10
                gx = -(49.*y.*(1918.*y - 11025.*y.^2 + 29155.*y.^3 - 36015.*y.^4 + 16807.*y.^5 - 120))/720;
                gy = (49.*x)/6 + (49931.*y)/180 - (46991.*x.*y)/180 + (36015.*x.*y.^2)/16 - (285719.*x.*y.^3)/36 + (588245.*x.*y.^4)/48 - (823543.*x.*y.^5)/120 - (634207.*y.^2)/240 + (98441.*y.^3)/9 - (1596665.*y.^4)/72 + (1294139.*y.^5)/60 - (5764801.*y.^6)/720 - 49/6;
            case 11
                gx = (49.*y.*(14.*x + 14.*y - 13).*(1715.*y.^2 - 350.*y - 3430.*y.^3 + 2401.*y.^4 + 24))/240;
                gy = (49.*(1715.*y.^2 - 350.*y - 3430.*y.^3 + 2401.*y.^4 + 24).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/240 + (49.*y.*(3430.*y - 10290.*y.^2 + 9604.*y.^3 - 350).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/240 + (49.*y.*(14.*x + 14.*y - 13).*(1715.*y.^2 - 350.*y - 3430.*y.^3 + 2401.*y.^4 + 24))/240;
            case 12
                gx = -(49.*y.*(77.*y - 294.*y.^2 + 343.*y.^3 - 6).*(294.*x.*y - 252.*y - 252.*x + 147.*x.^2 + 147.*y.^2 + 107))/144;
                gy = (5243.*x)/24 + 2009.*y - (789929.*x.^2.*y.^2)/16 + 117649.*x.^2.*y.^3 + (117649.*x.^3.*y.^2)/8 - (4117715.*x.^2.*y.^4)/48 - (823543.*x.^3.*y.^3)/36 - (477799.*x.*y)/72 + 52822.*x.*y.^2 + 7203.*x.^2.*y - (1495823.*x.*y.^3)/9 - (184877.*x.^3.*y)/72 + (1764735.*x.*y.^4)/8 - (823543.*x.*y.^5)/8 - (1029.*x.^2)/4 + (2401.*x.^3)/24 - (872935.*y.^2)/48 + (211288.*y.^3)/3 - (9495955.*y.^4)/72 + 117649.*y.^5 - (5764801.*y.^6)/144 - 245/4;
            case 13
                gx = (49.*y.*(49.*y.^2 - 21.*y + 2).*(2506.*x + 2506.*y - 6468.*x.*y + 4116.*x.*y.^2 + 4116.*x.^2.*y - 3234.*x.^2 + 1372.*x.^3 - 3234.*y.^2 + 1372.*y.^3 - 638))/144;
                gy = (6537923.*x.^2.*y.^2)/48 - (46501.*y)/18 - (15631.*x)/36 - (823543.*x.^2.*y.^3)/3 - (2000033.*x.^3.*y.^2)/24 + (4117715.*x.^2.*y.^4)/24 + (823543.*x.^3.*y.^3)/9 + (823543.*x.^4.*y.^2)/48 + (451045.*x.*y)/36 - (1106861.*x.*y.^2)/12 - (535423.*x.^2.*y)/24 + (789929.*x.*y.^3)/3 + (621859.*x.^3.*y)/36 - (7647185.*x.*y.^4)/24 - (117649.*x.^4.*y)/24 + (823543.*x.*y.^5)/6 + (61397.*x.^2)/72 - (26411.*x.^3)/36 + (16807.*x.^4)/72 + (133427.*y.^2)/6 - (2926819.*y.^3)/36 + (20756645.*y.^4)/144 - (2941225.*y.^5)/24 + (5764801.*y.^6)/144 + 245/3;
            case 14
                gx = -(49.*y.*(7.*y - 1).*(72030.*x.^2.*y.^2 - 16450.*y - 16450.*x + 72030.*x.*y - 102900.*x.*y.^2 - 102900.*x.^2.*y + 48020.*x.*y.^3 + 48020.*x.^3.*y + 36015.*x.^2 - 34300.*x.^3 + 12005.*x.^4 + 36015.*y.^2 - 34300.*y.^3 + 12005.*y.^4 + 2754))/240;
                gy = (22491.*x)/40 + (43071.*y)/20 - (2974839.*x.^2.*y.^2)/16 + (941192.*x.^2.*y.^3)/3 + (1294139.*x.^3.*y.^2)/8 - (4117715.*x.^2.*y.^4)/24 - (823543.*x.^3.*y.^3)/6 - (823543.*x.^4.*y.^2)/16 - (218834.*x.*y)/15 + (1481417.*x.*y.^2)/16 + (458591.*x.^2.*y)/12 - (2806769.*x.*y.^3)/12 - (386561.*x.^3.*y)/8 + (4117715.*x.*y.^4)/16 + (117649.*x.^4.*y)/4 - (823543.*x.*y.^5)/8 - (823543.*x.^5.*y)/120 - (80605.*x.^2)/48 + (117649.*x.^3)/48 - (84035.*x.^4)/48 + (117649.*x.^5)/240 - (1347647.*y.^2)/80 + (170471.*y.^3)/3 - (756315.*y.^4)/8 + (1529437.*y.^5)/20 - (5764801.*y.^6)/240 - 147/2;
            case 15
                gx = (49.*y.*(71456.*x + 71456.*y - 1944810.*x.^2.*y.^2 + 1008420.*x.^2.*y.^3 + 1008420.*x.^3.*y.^2 - 489510.*x.*y + 1214220.*x.*y.^2 + 1214220.*x.^2.*y - 1296540.*x.*y.^3 - 1296540.*x.^3.*y + 504210.*x.*y.^4 + 504210.*x.^4.*y - 244755.*x.^2 + 404740.*x.^3 - 324135.*x.^4 + 100842.*x.^5 - 244755.*y.^2 + 404740.*y.^3 - 324135.*y.^4 + 100842.*y.^5 - 8028))/720;
                gy = (991613.*x.^2.*y.^2)/8 - (10927.*y)/10 - (10927.*x)/20 - (352947.*x.^2.*y.^3)/2 - (1058841.*x.^3.*y.^2)/8 + (4117715.*x.^2.*y.^4)/48 + (823543.*x.^3.*y.^3)/9 + (823543.*x.^4.*y.^2)/16 + (437668.*x.*y)/45 - (799533.*x.*y.^2)/16 - (266511.*x.^2.*y)/8 + (991613.*x.*y.^3)/9 + (991613.*x.^3.*y)/18 - (1764735.*x.*y.^4)/16 - (352947.*x.^4.*y)/8 + (823543.*x.*y.^5)/20 + (823543.*x.^5.*y)/60 + (109417.*x.^2)/45 - (88837.*x.^3)/16 + (991613.*x.^4)/144 - (352947.*x.^5)/80 + (823543.*x.^6)/720 + (109417.*y.^2)/15 - (88837.*y.^3)/4 + (4958065.*y.^4)/144 - (1058841.*y.^5)/40 + (5764801.*y.^6)/720 + 49;
            case 16
                gx = (991613.*x.^2.*y.^2)/8 - (10927.*y)/20 - (10927.*x)/10 - (1058841.*x.^2.*y.^3)/8 - (352947.*x.^3.*y.^2)/2 + (823543.*x.^2.*y.^4)/16 + (823543.*x.^3.*y.^3)/9 + (4117715.*x.^4.*y.^2)/48 + (437668.*x.*y)/45 - (266511.*x.*y.^2)/8 - (799533.*x.^2.*y)/16 + (991613.*x.*y.^3)/18 + (991613.*x.^3.*y)/9 - (352947.*x.*y.^4)/8 - (1764735.*x.^4.*y)/16 + (823543.*x.*y.^5)/60 + (823543.*x.^5.*y)/20 + (109417.*x.^2)/15 - (88837.*x.^3)/4 + (4958065.*x.^4)/144 - (1058841.*x.^5)/40 + (5764801.*x.^6)/720 + (109417.*y.^2)/45 - (88837.*y.^3)/16 + (991613.*y.^4)/144 - (352947.*y.^5)/80 + (823543.*y.^6)/720 + 49;
                gy = (49.*x.*(71456.*x + 71456.*y - 1944810.*x.^2.*y.^2 + 1008420.*x.^2.*y.^3 + 1008420.*x.^3.*y.^2 - 489510.*x.*y + 1214220.*x.*y.^2 + 1214220.*x.^2.*y - 1296540.*x.*y.^3 - 1296540.*x.^3.*y + 504210.*x.*y.^4 + 504210.*x.^4.*y - 244755.*x.^2 + 404740.*x.^3 - 324135.*x.^4 + 100842.*x.^5 - 244755.*y.^2 + 404740.*y.^3 - 324135.*y.^4 + 100842.*y.^5 - 8028))/720;
            case 17
                gx = (43071.*x)/20 + (22491.*y)/40 - (2974839.*x.^2.*y.^2)/16 + (1294139.*x.^2.*y.^3)/8 + (941192.*x.^3.*y.^2)/3 - (823543.*x.^2.*y.^4)/16 - (823543.*x.^3.*y.^3)/6 - (4117715.*x.^4.*y.^2)/24 - (218834.*x.*y)/15 + (458591.*x.*y.^2)/12 + (1481417.*x.^2.*y)/16 - (386561.*x.*y.^3)/8 - (2806769.*x.^3.*y)/12 + (117649.*x.*y.^4)/4 + (4117715.*x.^4.*y)/16 - (823543.*x.*y.^5)/120 - (823543.*x.^5.*y)/8 - (1347647.*x.^2)/80 + (170471.*x.^3)/3 - (756315.*x.^4)/8 + (1529437.*x.^5)/20 - (5764801.*x.^6)/240 - (80605.*y.^2)/48 + (117649.*y.^3)/48 - (84035.*y.^4)/48 + (117649.*y.^5)/240 - 147/2;
                gy = -(49.*x.*(7.*x - 1).*(72030.*x.^2.*y.^2 - 16450.*y - 16450.*x + 72030.*x.*y - 102900.*x.*y.^2 - 102900.*x.^2.*y + 48020.*x.*y.^3 + 48020.*x.^3.*y + 36015.*x.^2 - 34300.*x.^3 + 12005.*x.^4 + 36015.*y.^2 - 34300.*y.^3 + 12005.*y.^4 + 2754))/240;
            case 18
                gx = (6537923.*x.^2.*y.^2)/48 - (15631.*y)/36 - (46501.*x)/18 - (2000033.*x.^2.*y.^3)/24 - (823543.*x.^3.*y.^2)/3 + (823543.*x.^2.*y.^4)/48 + (823543.*x.^3.*y.^3)/9 + (4117715.*x.^4.*y.^2)/24 + (451045.*x.*y)/36 - (535423.*x.*y.^2)/24 - (1106861.*x.^2.*y)/12 + (621859.*x.*y.^3)/36 + (789929.*x.^3.*y)/3 - (117649.*x.*y.^4)/24 - (7647185.*x.^4.*y)/24 + (823543.*x.^5.*y)/6 + (133427.*x.^2)/6 - (2926819.*x.^3)/36 + (20756645.*x.^4)/144 - (2941225.*x.^5)/24 + (5764801.*x.^6)/144 + (61397.*y.^2)/72 - (26411.*y.^3)/36 + (16807.*y.^4)/72 + 245/3;
                gy = (49.*x.*(49.*x.^2 - 21.*x + 2).*(2506.*x + 2506.*y - 6468.*x.*y + 4116.*x.*y.^2 + 4116.*x.^2.*y - 3234.*x.^2 + 1372.*x.^3 - 3234.*y.^2 + 1372.*y.^3 - 638))/144;
            case 19
                gx = 2009.*x + (5243.*y)/24 - (789929.*x.^2.*y.^2)/16 + (117649.*x.^2.*y.^3)/8 + 117649.*x.^3.*y.^2 - (823543.*x.^3.*y.^3)/36 - (4117715.*x.^4.*y.^2)/48 - (477799.*x.*y)/72 + 7203.*x.*y.^2 + 52822.*x.^2.*y - (184877.*x.*y.^3)/72 - (1495823.*x.^3.*y)/9 + (1764735.*x.^4.*y)/8 - (823543.*x.^5.*y)/8 - (872935.*x.^2)/48 + (211288.*x.^3)/3 - (9495955.*x.^4)/72 + 117649.*x.^5 - (5764801.*x.^6)/144 - (1029.*y.^2)/4 + (2401.*y.^3)/24 - 245/4;
                gy = -(49.*x.*(77.*x - 294.*x.^2 + 343.*x.^3 - 6).*(294.*x.*y - 252.*y - 252.*x + 147.*x.^2 + 147.*y.^2 + 107))/144;
            case 20
                gx = (49.*(1715.*x.^2 - 350.*x - 3430.*x.^3 + 2401.*x.^4 + 24).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/240 + (49.*x.*(3430.*x - 10290.*x.^2 + 9604.*x.^3 - 350).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/240 + (49.*x.*(14.*x + 14.*y - 13).*(1715.*x.^2 - 350.*x - 3430.*x.^3 + 2401.*x.^4 + 24))/240;
                gy = (49.*x.*(14.*x + 14.*y - 13).*(1715.*x.^2 - 350.*x - 3430.*x.^3 + 2401.*x.^4 + 24))/240;
            case 21
                gx = (49931.*x)/180 + (49.*y)/6 - (46991.*x.*y)/180 + (36015.*x.^2.*y)/16 - (285719.*x.^3.*y)/36 + (588245.*x.^4.*y)/48 - (823543.*x.^5.*y)/120 - (634207.*x.^2)/240 + (98441.*x.^3)/9 - (1596665.*x.^4)/72 + (1294139.*x.^5)/60 - (5764801.*x.^6)/720 - 49/6;
                gy = -(49.*x.*(1918.*x - 11025.*x.^2 + 29155.*x.^3 - 36015.*x.^4 + 16807.*x.^5 - 120))/720;
            case 22
                gx = -(343.*y.*(5508.*x + 2754.*y - 154350.*x.^2.*y.^2 + 72030.*x.^2.*y.^3 + 96040.*x.^3.*y.^2 - 32900.*x.*y + 72030.*x.*y.^2 + 108045.*x.^2.*y - 68600.*x.*y.^3 - 137200.*x.^3.*y + 24010.*x.*y.^4 + 60025.*x.^4.*y - 24675.*x.^2 + 48020.*x.^3 - 42875.*x.^4 + 14406.*x.^5 - 8225.*y.^2 + 12005.*y.^3 - 8575.*y.^4 + 2401.*y.^5 - 360))/120;
                gy = -(343.*x.*(2754.*x + 5508.*y - 154350.*x.^2.*y.^2 + 96040.*x.^2.*y.^3 + 72030.*x.^3.*y.^2 - 32900.*x.*y + 108045.*x.*y.^2 + 72030.*x.^2.*y - 137200.*x.*y.^3 - 68600.*x.^3.*y + 60025.*x.*y.^4 + 24010.*x.^4.*y - 8225.*x.^2 + 12005.*x.^3 - 8575.*x.^4 + 2401.*x.^5 - 24675.*y.^2 + 48020.*y.^3 - 42875.*y.^4 + 14406.*y.^5 - 360))/120;
            case 23
                gx = -(343.*y.*(748.*x + 24.*y - 700.*x.*y + 5145.*x.^2.*y - 13720.*x.^3.*y + 12005.*x.^4.*y - 6195.*x.^2 + 20580.*x.^3 - 29155.*x.^4 + 14406.*x.^5 - 24))/120;
                gy = -(343.*x.*(x + 2.*y - 1).*(1715.*x.^2 - 350.*x - 3430.*x.^3 + 2401.*x.^4 + 24))/120;
            case 24
                gx = -(343.*y.*(2.*x + y - 1).*(1715.*y.^2 - 350.*y - 3430.*y.^3 + 2401.*y.^4 + 24))/120;
                gy = -(343.*x.*(24.*x + 748.*y - 700.*x.*y + 5145.*x.*y.^2 - 13720.*x.*y.^3 + 12005.*x.*y.^4 - 6195.*y.^2 + 20580.*y.^3 - 29155.*y.^4 + 14406.*y.^5 - 24))/120;
            case 25
                gx = -(343.*y.*(7.*y - 1).*(154.*x.*y - 6.*y - 166.*x - 882.*x.^2.*y + 1372.*x.^3.*y + 1113.*x.^2 - 2548.*x.^3 + 1715.*x.^4 + 6))/48;
                gy = -(343.*x.*(77.*x - 294.*x.^2 + 343.*x.^3 - 6).*(14.*x.*y - 16.*y - x + 21.*y.^2 + 1))/48;
            case 26
                gx = -(343.*y.*(49.*y.^2 - 21.*y + 2).*(46.*x + 2.*y - 42.*x.*y + 147.*x.^2.*y - 210.*x.^2 + 196.*x.^3 - 2))/36;
                gy = -(343.*x.*(49.*x.^2 - 21.*x + 2).*(2.*x + 46.*y - 42.*x.*y + 147.*x.*y.^2 - 210.*y.^2 + 196.*y.^3 - 2))/36;
            case 27
                gx = -(343.*y.*(77.*y - 294.*y.^2 + 343.*y.^3 - 6).*(14.*x.*y - y - 16.*x + 21.*x.^2 + 1))/48;
                gy = -(343.*x.*(7.*x - 1).*(154.*x.*y - 166.*y - 6.*x - 882.*x.*y.^2 + 1372.*x.*y.^3 + 1113.*y.^2 - 2548.*y.^3 + 1715.*y.^4 + 6))/48;
            case 28
                gx = (343.*y.*(77.*y - 294.*y.^2 + 343.*y.^3 - 6).*(28.*x.*y - 13.*y - 26.*x + 21.*x.^2 + 7.*y.^2 + 6))/48;
                gy = (343.*x.*(77.*y - 294.*y.^2 + 343.*y.^3 - 6).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/48 + (343.*x.*y.*(14.*x + 14.*y - 13).*(77.*y - 294.*y.^2 + 343.*y.^3 - 6))/48 + (343.*x.*y.*(1029.*y.^2 - 588.*y + 77).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/48;
            case 29
                gx = -(343.*y.*(49.*y.^2 - 21.*y + 2).*(214.*x + 107.*y - 504.*x.*y + 294.*x.*y.^2 + 441.*x.^2.*y - 378.*x.^2 + 196.*x.^3 - 126.*y.^2 + 49.*y.^3 - 30))/36;
                gy = -(343.*x.*(214.*x + 1688.*y - 27783.*x.^2.*y.^2 + 28812.*x.^2.*y.^3 + 7203.*x.^3.*y.^2 - 5502.*x.*y + 32487.*x.*y.^2 + 5880.*x.^2.*y - 61740.*x.*y.^3 - 2058.*x.^3.*y + 36015.*x.*y.^4 - 252.*x.^2 + 98.*x.^3 - 11907.*y.^2 + 31948.*y.^3 - 36015.*y.^4 + 14406.*y.^5 - 60))/36;
            case 30
                gx = (343.*y.*(7.*y - 1).*(6174.*x.^2.*y.^2 - 638.*y - 1276.*x + 5012.*x.*y - 6468.*x.*y.^2 - 9702.*x.^2.*y + 2744.*x.*y.^3 + 5488.*x.^3.*y + 3759.*x.^2 - 4312.*x.^3 + 1715.*x.^4 + 1253.*y.^2 - 1078.*y.^3 + 343.*y.^4 + 120))/48;
                gy = (343.*x.*(638.*x + 2956.*y - 74088.*x.^2.*y.^2 + 57624.*x.^2.*y.^3 + 28812.*x.^3.*y.^2 - 13944.*x.*y + 62328.*x.*y.^2 + 24010.*x.^2.*y - 96040.*x.*y.^3 - 17836.*x.^3.*y + 48020.*x.*y.^4 + 4802.*x.^4.*y - 1253.*x.^2 + 1078.*x.^3 - 343.*x.^4 - 17157.*y.^2 + 39396.*y.^3 - 39445.*y.^4 + 14406.*y.^5 - 120))/48;
            case 31
                gx = (343.*y.*(2956.*x + 638.*y - 74088.*x.^2.*y.^2 + 28812.*x.^2.*y.^3 + 57624.*x.^3.*y.^2 - 13944.*x.*y + 24010.*x.*y.^2 + 62328.*x.^2.*y - 17836.*x.*y.^3 - 96040.*x.^3.*y + 4802.*x.*y.^4 + 48020.*x.^4.*y - 17157.*x.^2 + 39396.*x.^3 - 39445.*x.^4 + 14406.*x.^5 - 1253.*y.^2 + 1078.*y.^3 - 343.*y.^4 - 120))/48;
                gy = (343.*x.*(7.*x - 1).*(6174.*x.^2.*y.^2 - 1276.*y - 638.*x + 5012.*x.*y - 9702.*x.*y.^2 - 6468.*x.^2.*y + 5488.*x.*y.^3 + 2744.*x.^3.*y + 1253.*x.^2 - 1078.*x.^3 + 343.*x.^4 + 3759.*y.^2 - 4312.*y.^3 + 1715.*y.^4 + 120))/48;
            case 32
                gx = -(343.*y.*(1688.*x + 214.*y - 27783.*x.^2.*y.^2 + 7203.*x.^2.*y.^3 + 28812.*x.^3.*y.^2 - 5502.*x.*y + 5880.*x.*y.^2 + 32487.*x.^2.*y - 2058.*x.*y.^3 - 61740.*x.^3.*y + 36015.*x.^4.*y - 11907.*x.^2 + 31948.*x.^3 - 36015.*x.^4 + 14406.*x.^5 - 252.*y.^2 + 98.*y.^3 - 60))/36;
                gy = -(343.*x.*(49.*x.^2 - 21.*x + 2).*(107.*x + 214.*y - 504.*x.*y + 441.*x.*y.^2 + 294.*x.^2.*y - 126.*x.^2 + 49.*x.^3 - 378.*y.^2 + 196.*y.^3 - 30))/36;
            case 33
                gx = (343.*y.*(77.*x - 294.*x.^2 + 343.*x.^3 - 6).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/48 + (343.*x.*y.*(14.*x + 14.*y - 13).*(77.*x - 294.*x.^2 + 343.*x.^3 - 6))/48 + (343.*x.*y.*(1029.*x.^2 - 588.*x + 77).*(14.*x.*y - 13.*y - 13.*x + 7.*x.^2 + 7.*y.^2 + 6))/48;
                gy = (343.*x.*(77.*x - 294.*x.^2 + 343.*x.^3 - 6).*(28.*x.*y - 26.*y - 13.*x + 7.*x.^2 + 21.*y.^2 + 6))/48;
            case 34
                gx = -(343.*y.*(7.*y - 1).*(3087.*x.^2.*y.^2 - 107.*y - 634.*x + 2002.*x.*y - 2058.*x.*y.^2 - 5733.*x.^2.*y + 686.*x.*y.^3 + 4116.*x.^3.*y + 2625.*x.^2 - 3724.*x.^3 + 1715.*x.^4 + 126.*y.^2 - 49.*y.^3 + 30))/24;
                gy = -(343.*x.*(7.*x - 1).*(3087.*x.^2.*y.^2 - 634.*y - 107.*x + 2002.*x.*y - 5733.*x.*y.^2 - 2058.*x.^2.*y + 4116.*x.*y.^3 + 686.*x.^3.*y + 126.*x.^2 - 49.*x.^3 + 2625.*y.^2 - 3724.*y.^3 + 1715.*y.^4 + 30))/24;
            case 35
                gx = (343.*y.*(7.*y - 1).*(1029.*x.^2.*y.^2 - 26.*y - 304.*x + 602.*x.*y - 294.*x.*y.^2 - 2793.*x.^2.*y + 2744.*x.^3.*y + 1743.*x.^2 - 3136.*x.^3 + 1715.*x.^4 + 14.*y.^2 + 12))/24;
                gy = (343.*x.*(49.*x.^2 - 21.*x + 2).*(13.*x + 110.*y - 210.*x.*y + 294.*x.*y.^2 + 98.*x.^2.*y - 7.*x.^2 - 294.*y.^2 + 196.*y.^3 - 6))/24;
            case 36
                gx = (343.*y.*(49.*y.^2 - 21.*y + 2).*(110.*x + 13.*y - 210.*x.*y + 98.*x.*y.^2 + 294.*x.^2.*y - 294.*x.^2 + 196.*x.^3 - 7.*y.^2 - 6))/24;
                gy = (343.*x.*(7.*x - 1).*(1029.*x.^2.*y.^2 - 304.*y - 26.*x + 602.*x.*y - 2793.*x.*y.^2 - 294.*x.^2.*y + 2744.*x.*y.^3 + 14.*x.^2 + 1743.*y.^2 - 3136.*y.^3 + 1715.*y.^4 + 12))/24;
        end
    case 8
        switch i
            case 1
                gx = (118124.*x)/315 + (118124.*y)/315 - 110592.*x.^2.*y.^2 + 212992.*x.^2.*y.^3 + 212992.*x.^3.*y.^2 - 196608.*x.^2.*y.^4 - 262144.*x.^3.*y.^3 - 196608.*x.^4.*y.^2 + (1048576.*x.^2.*y.^5)/15 + (1048576.*x.^3.*y.^4)/9 + (1048576.*x.^4.*y.^3)/9 + (1048576.*x.^5.*y.^2)/15 - (25632.*x.*y)/5 + (136832.*x.*y.^2)/5 + (136832.*x.^2.*y)/5 - 73728.*x.*y.^3 - 73728.*x.^3.*y + 106496.*x.*y.^4 + 106496.*x.^4.*y - (393216.*x.*y.^5)/5 - (393216.*x.^5.*y)/5 + (1048576.*x.*y.^6)/45 + (1048576.*x.^6.*y)/45 - (12816.*x.^2)/5 + (136832.*x.^3)/15 - 18432.*x.^4 + (106496.*x.^5)/5 - (65536.*x.^6)/5 + (1048576.*x.^7)/315 - (12816.*y.^2)/5 + (136832.*y.^3)/15 - 18432.*y.^4 + (106496.*y.^5)/5 - (65536.*y.^6)/5 + (1048576.*y.^7)/315 - 761/35;
                gy = (118124.*x)/315 + (118124.*y)/315 - 110592.*x.^2.*y.^2 + 212992.*x.^2.*y.^3 + 212992.*x.^3.*y.^2 - 196608.*x.^2.*y.^4 - 262144.*x.^3.*y.^3 - 196608.*x.^4.*y.^2 + (1048576.*x.^2.*y.^5)/15 + (1048576.*x.^3.*y.^4)/9 + (1048576.*x.^4.*y.^3)/9 + (1048576.*x.^5.*y.^2)/15 - (25632.*x.*y)/5 + (136832.*x.*y.^2)/5 + (136832.*x.^2.*y)/5 - 73728.*x.*y.^3 - 73728.*x.^3.*y + 106496.*x.*y.^4 + 106496.*x.^4.*y - (393216.*x.*y.^5)/5 - (393216.*x.^5.*y)/5 + (1048576.*x.*y.^6)/45 + (1048576.*x.^6.*y)/45 - (12816.*x.^2)/5 + (136832.*x.^3)/15 - 18432.*x.^4 + (106496.*x.^5)/5 - (65536.*x.^6)/5 + (1048576.*x.^7)/315 - (12816.*y.^2)/5 + (136832.*y.^3)/15 - 18432.*y.^4 + (106496.*y.^5)/5 - (65536.*y.^6)/5 + (1048576.*y.^7)/315 - 761/35;
            case 2
                gx = (1452.*x)/35 - (7504.*x.^2)/15 + (123776.*x.^3)/45 - (71680.*x.^4)/9 + (188416.*x.^5)/15 - (458752.*x.^6)/45 + (1048576.*x.^7)/315 - 1;
                gy = 0;
            case 3
                gx = 0;
                gy = (1452.*y)/35 - (7504.*y.^2)/15 + (123776.*y.^3)/45 - (71680.*y.^4)/9 + (188416.*y.^5)/15 - (458752.*y.^6)/45 + (1048576.*y.^7)/315 - 1;
            case 4
                gx = (64.*y.*(19488.*x.^2 - 1764.*x - 94080.*x.^3 + 224000.*x.^4 - 258048.*x.^5 + 114688.*x.^6 + 45))/315;
                gy = (64.*x.*(6496.*x.^2 - 882.*x - 23520.*x.^3 + 44800.*x.^4 - 43008.*x.^5 + 16384.*x.^6 + 45))/315;
            case 5
                gx = (16.*y.*(8.*y - 1).*(548.*x - 5400.*x.^2 + 21760.*x.^3 - 38400.*x.^4 + 24576.*x.^5 - 15))/45;
                gy = (16.*x.*(16.*y - 1).*(274.*x - 1800.*x.^2 + 5440.*x.^3 - 7680.*x.^4 + 4096.*x.^5 - 15))/45;
            case 6
                gx = (64.*y.*(32.*y.^2 - 12.*y + 1).*(840.*x.^2 - 100.*x - 2560.*x.^3 + 2560.*x.^4 + 3))/45;
                gy = (64.*x.*(96.*y.^2 - 24.*y + 1).*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3))/45;
            case 7
                gx = (4.*y.*(88.*x - 576.*x.^2 + 1024.*x.^3 - 3).*(44.*y - 192.*y.^2 + 256.*y.^3 - 3))/9;
                gy = (4.*x.*(44.*x - 192.*x.^2 + 256.*x.^3 - 3).*(88.*y - 576.*y.^2 + 1024.*y.^3 - 3))/9;
            case 8
                gx = (64.*y.*(96.*x.^2 - 24.*x + 1).*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3))/45;
                gy = (64.*x.*(32.*x.^2 - 12.*x + 1).*(840.*y.^2 - 100.*y - 2560.*y.^3 + 2560.*y.^4 + 3))/45;
            case 9
                gx = (16.*y.*(16.*x - 1).*(274.*y - 1800.*y.^2 + 5440.*y.^3 - 7680.*y.^4 + 4096.*y.^5 - 15))/45;
                gy = (16.*x.*(8.*x - 1).*(548.*y - 5400.*y.^2 + 21760.*y.^3 - 38400.*y.^4 + 24576.*y.^5 - 15))/45;
            case 10
                gx = (64.*y.*(6496.*y.^2 - 882.*y - 23520.*y.^3 + 44800.*y.^4 - 43008.*y.^5 + 16384.*y.^6 + 45))/315;
                gy = (64.*x.*(19488.*y.^2 - 1764.*y - 94080.*y.^3 + 224000.*y.^4 - 258048.*y.^5 + 114688.*y.^6 + 45))/315;
            case 11
                gx = -(64.*y.*(6496.*y.^2 - 882.*y - 23520.*y.^3 + 44800.*y.^4 - 43008.*y.^5 + 16384.*y.^6 + 45))/315;
                gy = (1792.*x.*y)/5 - (13184.*y)/35 - (64.*x)/7 - (59392.*x.*y.^2)/15 + (57344.*x.*y.^3)/3 - (409600.*x.*y.^4)/9 + (262144.*x.*y.^5)/5 - (1048576.*x.*y.^6)/45 + (67456.*y.^2)/15 - (1097728.*y.^3)/45 + (624640.*y.^4)/9 - (1605632.*y.^5)/15 + (3801088.*y.^6)/45 - (8388608.*y.^7)/315 + 64/7;
            case 12
                gx = (16.*y.*(16.*x + 16.*y - 15).*(274.*y - 1800.*y.^2 + 5440.*y.^3 - 7680.*y.^4 + 4096.*y.^5 - 15))/45;
                gy = (16.*(274.*y - 1800.*y.^2 + 5440.*y.^3 - 7680.*y.^4 + 4096.*y.^5 - 15).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/45 + (16.*y.*(16320.*y.^2 - 3600.*y - 30720.*y.^3 + 20480.*y.^4 + 274).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/45 + (16.*y.*(16.*x + 16.*y - 15).*(274.*y - 1800.*y.^2 + 5440.*y.^3 - 7680.*y.^4 + 4096.*y.^5 - 15))/45;
            case 13
                gx = -(64.*y.*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3).*(192.*x.*y - 168.*y - 168.*x + 96.*x.^2 + 96.*y.^2 + 73))/45;
                gy = - (64.*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3).*(73.*x + 73.*y - 168.*x.*y + 96.*x.*y.^2 + 96.*x.^2.*y - 84.*x.^2 + 32.*x.^3 - 84.*y.^2 + 32.*y.^3 - 21))/45 - (64.*y.*(560.*y - 1920.*y.^2 + 2048.*y.^3 - 50).*(73.*x + 73.*y - 168.*x.*y + 96.*x.*y.^2 + 96.*x.^2.*y - 84.*x.^2 + 32.*x.^3 - 84.*y.^2 + 32.*y.^3 - 21))/45 - (64.*y.*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3).*(192.*x.*y - 168.*y - 168.*x + 96.*x.^2 + 96.*y.^2 + 73))/45;
            case 14
                gx = (4.*y.*(44.*y - 192.*y.^2 + 256.*y.^3 - 3).*(2008.*x + 2008.*y - 4992.*x.*y + 3072.*x.*y.^2 + 3072.*x.^2.*y - 2496.*x.^2 + 1024.*x.^3 - 2496.*y.^2 + 1024.*y.^3 - 533))/9;
                gy = (2132.*x)/3 + 5528.*y - 409600.*x.^2.*y.^2 + (12861440.*x.^2.*y.^3)/9 + (819200.*x.^3.*y.^2)/3 - (6225920.*x.^2.*y.^4)/3 - (6553600.*x.^3.*y.^3)/9 - 65536.*x.^4.*y.^2 + 1048576.*x.^2.*y.^5 + (5242880.*x.^3.*y.^4)/9 + (1048576.*x.^4.*y.^3)/9 - (235808.*x.*y)/9 + (792704.*x.*y.^2)/3 + (413312.*x.^2.*y)/9 - (10158080.*x.*y.^3)/9 - (317440.*x.^3.*y)/9 + (6922240.*x.*y.^4)/3 + (90112.*x.^4.*y)/9 - 2228224.*x.*y.^5 + (7340032.*x.*y.^6)/9 - (4016.*x.^2)/3 + (3328.*x.^3)/3 - (1024.*x.^4)/3 - (186496.*y.^2)/3 + (2814208.*y.^3)/9 - (7331840.*y.^4)/9 + (3424256.*y.^5)/3 - (7340032.*y.^6)/9 + (2097152.*y.^7)/9 - 140;
            case 15
                gx = -(64.*y.*(32.*y.^2 - 12.*y + 1).*(15360.*x.^2.*y.^2 - 4140.*y - 4140.*x + 17040.*x.*y - 23040.*x.*y.^2 - 23040.*x.^2.*y + 10240.*x.*y.^3 + 10240.*x.^3.*y + 8520.*x.^2 - 7680.*x.^3 + 2560.*x.^4 + 8520.*y.^2 - 7680.*y.^3 + 2560.*y.^4 + 743))/45;
                gy = 768000.*x.^2.*y.^2 - (256384.*y)/45 - (47552.*x)/45 - (21299200.*x.^2.*y.^3)/9 - 802816.*x.^3.*y.^2 + (9175040.*x.^2.*y.^4)/3 + (5242880.*x.^3.*y.^3)/3 + 393216.*x.^4.*y.^2 - (4194304.*x.^2.*y.^5)/3 - (10485760.*x.^3.*y.^4)/9 - (4194304.*x.^4.*y.^3)/9 - (1048576.*x.^5.*y.^2)/15 + (557056.*x.*y)/15 - (5246464.*x.*y.^2)/15 - (284672.*x.^2.*y)/3 + (4136960.*x.*y.^3)/3 + 118784.*x.^3.*y - (23511040.*x.*y.^4)/9 - (655360.*x.^4.*y)/9 + 2359296.*x.*y.^5 + (262144.*x.^5.*y)/15 - (7340032.*x.*y.^6)/9 + 2944.*x.^2 - (36352.*x.^3)/9 + (8192.*x.^4)/3 - (32768.*x.^5)/45 + (306048.*y.^2)/5 - (4390912.*y.^3)/15 + 727040.*y.^4 - (4882432.*y.^5)/5 + (10092544.*y.^6)/15 - (8388608.*y.^7)/45 + 448/3;
            case 16
                gx = (16.*y.*(8.*y - 1).*(24308.*x + 24308.*y - 506880.*x.^2.*y.^2 + 245760.*x.^2.*y.^3 + 245760.*x.^3.*y.^2 - 150480.*x.*y + 341760.*x.*y.^2 + 341760.*x.^2.*y - 337920.*x.*y.^3 - 337920.*x.^3.*y + 122880.*x.*y.^4 + 122880.*x.^4.*y - 75240.*x.^2 + 113920.*x.^3 - 84480.*x.^4 + 24576.*x.^5 - 75240.*y.^2 + 113920.*y.^3 - 84480.*y.^4 + 24576.*y.^5 - 3069))/45;
                gy = (5456.*x)/5 + (19872.*y)/5 - 824320.*x.^2.*y.^2 + (6553600.*x.^2.*y.^3)/3 + (3457024.*x.^3.*y.^2)/3 - (7536640.*x.^2.*y.^4)/3 - (18350080.*x.^3.*y.^3)/9 - 786432.*x.^4.*y.^2 + 1048576.*x.^2.*y.^5 + (10485760.*x.^3.*y.^4)/9 + (2097152.*x.^4.*y.^3)/3 + (1048576.*x.^5.*y.^2)/5 - (312704.*x.*y)/9 + (4315264.*x.*y.^2)/15 + (5519104.*x.^2.*y)/45 - (9162752.*x.*y.^3)/9 - (2013184.*x.^3.*y)/9 + (15933440.*x.*y.^4)/9 + (1998848.*x.^4.*y)/9 - (7471104.*x.*y.^5)/5 - (1703936.*x.^5.*y)/15 + (7340032.*x.*y.^6)/15 + (1048576.*x.^6.*y)/45 - (194464.*x.^2)/45 + (26752.*x.^3)/3 - (91136.*x.^4)/9 + (90112.*x.^5)/15 - (65536.*x.^6)/45 - (587296.*y.^2)/15 + (7827968.*y.^3)/45 - (3665920.*y.^4)/9 + (7831552.*y.^5)/15 - (15597568.*y.^6)/45 + (4194304.*y.^7)/45 - 112;
            case 17
                gx = -(64.*y.*(3924480.*x.^2.*y.^2 - 48860.*y - 48860.*x - 4300800.*x.^2.*y.^3 - 4300800.*x.^3.*y.^2 + 1720320.*x.^2.*y.^4 + 2293760.*x.^3.*y.^3 + 1720320.*x.^4.*y.^2 + 442176.*x.*y - 1545600.*x.*y.^2 - 1545600.*x.^2.*y + 2616320.*x.*y.^3 + 2616320.*x.^3.*y - 2150400.*x.*y.^4 - 2150400.*x.^4.*y + 688128.*x.*y.^5 + 688128.*x.^5.*y + 221088.*x.^2 - 515200.*x.^3 + 654080.*x.^4 - 430080.*x.^5 + 114688.*x.^6 + 221088.*y.^2 - 515200.*y.^3 + 654080.*y.^4 - 430080.*y.^5 + 114688.*y.^6 + 4329))/315;
                gy = 471040.*x.^2.*y.^2 - (61568.*y)/35 - (30784.*x)/35 - (9568256.*x.^2.*y.^3)/9 - (2392064.*x.^3.*y.^2)/3 + (3276800.*x.^2.*y.^4)/3 + (10485760.*x.^3.*y.^3)/9 + 655360.*x.^4.*y.^2 - (2097152.*x.^2.*y.^5)/5 - (5242880.*x.^3.*y.^4)/9 - (4194304.*x.^4.*y.^3)/9 - (1048576.*x.^5.*y.^2)/5 + (178688.*x.*y)/9 - (673792.*x.*y.^2)/5 - (1347584.*x.^2.*y)/15 + (3768320.*x.*y.^3)/9 + (1884160.*x.^3.*y)/9 - (5980160.*x.*y.^4)/9 - (2392064.*x.^4.*y)/9 + 524288.*x.*y.^5 + (524288.*x.^5.*y)/3 - (7340032.*x.*y.^6)/45 - (2097152.*x.^6.*y)/45 + (44672.*x.^2)/9 - (673792.*x.^3)/45 + (235520.*x.^4)/9 - (1196032.*x.^5)/45 + (131072.*x.^6)/9 - (1048576.*x.^7)/315 + (44672.*y.^2)/3 - (2695168.*y.^3)/45 + (1177600.*y.^4)/9 - (2392064.*y.^5)/15 + (917504.*y.^6)/9 - (8388608.*y.^7)/315 + 64;
            case 18
                gx = 471040.*x.^2.*y.^2 - (30784.*y)/35 - (61568.*x)/35 - (2392064.*x.^2.*y.^3)/3 - (9568256.*x.^3.*y.^2)/9 + 655360.*x.^2.*y.^4 + (10485760.*x.^3.*y.^3)/9 + (3276800.*x.^4.*y.^2)/3 - (1048576.*x.^2.*y.^5)/5 - (4194304.*x.^3.*y.^4)/9 - (5242880.*x.^4.*y.^3)/9 - (2097152.*x.^5.*y.^2)/5 + (178688.*x.*y)/9 - (1347584.*x.*y.^2)/15 - (673792.*x.^2.*y)/5 + (1884160.*x.*y.^3)/9 + (3768320.*x.^3.*y)/9 - (2392064.*x.*y.^4)/9 - (5980160.*x.^4.*y)/9 + (524288.*x.*y.^5)/3 + 524288.*x.^5.*y - (2097152.*x.*y.^6)/45 - (7340032.*x.^6.*y)/45 + (44672.*x.^2)/3 - (2695168.*x.^3)/45 + (1177600.*x.^4)/9 - (2392064.*x.^5)/15 + (917504.*x.^6)/9 - (8388608.*x.^7)/315 + (44672.*y.^2)/9 - (673792.*y.^3)/45 + (235520.*y.^4)/9 - (1196032.*y.^5)/45 + (131072.*y.^6)/9 - (1048576.*y.^7)/315 + 64;
                gy = -(64.*x.*(3924480.*x.^2.*y.^2 - 48860.*y - 48860.*x - 4300800.*x.^2.*y.^3 - 4300800.*x.^3.*y.^2 + 1720320.*x.^2.*y.^4 + 2293760.*x.^3.*y.^3 + 1720320.*x.^4.*y.^2 + 442176.*x.*y - 1545600.*x.*y.^2 - 1545600.*x.^2.*y + 2616320.*x.*y.^3 + 2616320.*x.^3.*y - 2150400.*x.*y.^4 - 2150400.*x.^4.*y + 688128.*x.*y.^5 + 688128.*x.^5.*y + 221088.*x.^2 - 515200.*x.^3 + 654080.*x.^4 - 430080.*x.^5 + 114688.*x.^6 + 221088.*y.^2 - 515200.*y.^3 + 654080.*y.^4 - 430080.*y.^5 + 114688.*y.^6 + 4329))/315;
            case 19
                gx = (19872.*x)/5 + (5456.*y)/5 - 824320.*x.^2.*y.^2 + (3457024.*x.^2.*y.^3)/3 + (6553600.*x.^3.*y.^2)/3 - 786432.*x.^2.*y.^4 - (18350080.*x.^3.*y.^3)/9 - (7536640.*x.^4.*y.^2)/3 + (1048576.*x.^2.*y.^5)/5 + (2097152.*x.^3.*y.^4)/3 + (10485760.*x.^4.*y.^3)/9 + 1048576.*x.^5.*y.^2 - (312704.*x.*y)/9 + (5519104.*x.*y.^2)/45 + (4315264.*x.^2.*y)/15 - (2013184.*x.*y.^3)/9 - (9162752.*x.^3.*y)/9 + (1998848.*x.*y.^4)/9 + (15933440.*x.^4.*y)/9 - (1703936.*x.*y.^5)/15 - (7471104.*x.^5.*y)/5 + (1048576.*x.*y.^6)/45 + (7340032.*x.^6.*y)/15 - (587296.*x.^2)/15 + (7827968.*x.^3)/45 - (3665920.*x.^4)/9 + (7831552.*x.^5)/15 - (15597568.*x.^6)/45 + (4194304.*x.^7)/45 - (194464.*y.^2)/45 + (26752.*y.^3)/3 - (91136.*y.^4)/9 + (90112.*y.^5)/15 - (65536.*y.^6)/45 - 112;
                gy = (16.*x.*(8.*x - 1).*(24308.*x + 24308.*y - 506880.*x.^2.*y.^2 + 245760.*x.^2.*y.^3 + 245760.*x.^3.*y.^2 - 150480.*x.*y + 341760.*x.*y.^2 + 341760.*x.^2.*y - 337920.*x.*y.^3 - 337920.*x.^3.*y + 122880.*x.*y.^4 + 122880.*x.^4.*y - 75240.*x.^2 + 113920.*x.^3 - 84480.*x.^4 + 24576.*x.^5 - 75240.*y.^2 + 113920.*y.^3 - 84480.*y.^4 + 24576.*y.^5 - 3069))/45;
            case 20
                gx = 768000.*x.^2.*y.^2 - (47552.*y)/45 - (256384.*x)/45 - 802816.*x.^2.*y.^3 - (21299200.*x.^3.*y.^2)/9 + 393216.*x.^2.*y.^4 + (5242880.*x.^3.*y.^3)/3 + (9175040.*x.^4.*y.^2)/3 - (1048576.*x.^2.*y.^5)/15 - (4194304.*x.^3.*y.^4)/9 - (10485760.*x.^4.*y.^3)/9 - (4194304.*x.^5.*y.^2)/3 + (557056.*x.*y)/15 - (284672.*x.*y.^2)/3 - (5246464.*x.^2.*y)/15 + 118784.*x.*y.^3 + (4136960.*x.^3.*y)/3 - (655360.*x.*y.^4)/9 - (23511040.*x.^4.*y)/9 + (262144.*x.*y.^5)/15 + 2359296.*x.^5.*y - (7340032.*x.^6.*y)/9 + (306048.*x.^2)/5 - (4390912.*x.^3)/15 + 727040.*x.^4 - (4882432.*x.^5)/5 + (10092544.*x.^6)/15 - (8388608.*x.^7)/45 + 2944.*y.^2 - (36352.*y.^3)/9 + (8192.*y.^4)/3 - (32768.*y.^5)/45 + 448/3;
                gy = -(64.*x.*(32.*x.^2 - 12.*x + 1).*(15360.*x.^2.*y.^2 - 4140.*y - 4140.*x + 17040.*x.*y - 23040.*x.*y.^2 - 23040.*x.^2.*y + 10240.*x.*y.^3 + 10240.*x.^3.*y + 8520.*x.^2 - 7680.*x.^3 + 2560.*x.^4 + 8520.*y.^2 - 7680.*y.^3 + 2560.*y.^4 + 743))/45;
            case 21
                gx = 5528.*x + (2132.*y)/3 - 409600.*x.^2.*y.^2 + (819200.*x.^2.*y.^3)/3 + (12861440.*x.^3.*y.^2)/9 - 65536.*x.^2.*y.^4 - (6553600.*x.^3.*y.^3)/9 - (6225920.*x.^4.*y.^2)/3 + (1048576.*x.^3.*y.^4)/9 + (5242880.*x.^4.*y.^3)/9 + 1048576.*x.^5.*y.^2 - (235808.*x.*y)/9 + (413312.*x.*y.^2)/9 + (792704.*x.^2.*y)/3 - (317440.*x.*y.^3)/9 - (10158080.*x.^3.*y)/9 + (90112.*x.*y.^4)/9 + (6922240.*x.^4.*y)/3 - 2228224.*x.^5.*y + (7340032.*x.^6.*y)/9 - (186496.*x.^2)/3 + (2814208.*x.^3)/9 - (7331840.*x.^4)/9 + (3424256.*x.^5)/3 - (7340032.*x.^6)/9 + (2097152.*x.^7)/9 - (4016.*y.^2)/3 + (3328.*y.^3)/3 - (1024.*y.^4)/3 - 140;
                gy = (4.*x.*(44.*x - 192.*x.^2 + 256.*x.^3 - 3).*(2008.*x + 2008.*y - 4992.*x.*y + 3072.*x.*y.^2 + 3072.*x.^2.*y - 2496.*x.^2 + 1024.*x.^3 - 2496.*y.^2 + 1024.*y.^3 - 533))/9;
            case 22
                gx = - (64.*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3).*(73.*x + 73.*y - 168.*x.*y + 96.*x.*y.^2 + 96.*x.^2.*y - 84.*x.^2 + 32.*x.^3 - 84.*y.^2 + 32.*y.^3 - 21))/45 - (64.*x.*(560.*x - 1920.*x.^2 + 2048.*x.^3 - 50).*(73.*x + 73.*y - 168.*x.*y + 96.*x.*y.^2 + 96.*x.^2.*y - 84.*x.^2 + 32.*x.^3 - 84.*y.^2 + 32.*y.^3 - 21))/45 - (64.*x.*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3).*(192.*x.*y - 168.*y - 168.*x + 96.*x.^2 + 96.*y.^2 + 73))/45;
                gy = -(64.*x.*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3).*(192.*x.*y - 168.*y - 168.*x + 96.*x.^2 + 96.*y.^2 + 73))/45;
            case 23
                gx = (16.*(274.*x - 1800.*x.^2 + 5440.*x.^3 - 7680.*x.^4 + 4096.*x.^5 - 15).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/45 + (16.*x.*(16320.*x.^2 - 3600.*x - 30720.*x.^3 + 20480.*x.^4 + 274).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/45 + (16.*x.*(16.*x + 16.*y - 15).*(274.*x - 1800.*x.^2 + 5440.*x.^3 - 7680.*x.^4 + 4096.*x.^5 - 15))/45;
                gy = (16.*x.*(16.*x + 16.*y - 15).*(274.*x - 1800.*x.^2 + 5440.*x.^3 - 7680.*x.^4 + 4096.*x.^5 - 15))/45;
            case 24
                gx = (1792.*x.*y)/5 - (64.*y)/7 - (13184.*x)/35 - (59392.*x.^2.*y)/15 + (57344.*x.^3.*y)/3 - (409600.*x.^4.*y)/9 + (262144.*x.^5.*y)/5 - (1048576.*x.^6.*y)/45 + (67456.*x.^2)/15 - (1097728.*x.^3)/45 + (624640.*x.^4)/9 - (1605632.*x.^5)/15 + (3801088.*x.^6)/45 - (8388608.*x.^7)/315 + 64/7;
                gy = -(64.*x.*(6496.*x.^2 - 882.*x - 23520.*x.^3 + 44800.*x.^4 - 43008.*x.^5 + 16384.*x.^6 + 45))/315;
            case 25
                gx = (256.*y.*(512640.*x.^2.*y.^2 - 3069.*y - 6138.*x - 506880.*x.^2.*y.^3 - 675840.*x.^3.*y.^2 + 184320.*x.^2.*y.^4 + 327680.*x.^3.*y.^3 + 307200.*x.^4.*y.^2 + 48616.*x.*y - 150480.*x.*y.^2 - 225720.*x.^2.*y + 227840.*x.*y.^3 + 455680.*x.^3.*y - 168960.*x.*y.^4 - 422400.*x.^4.*y + 49152.*x.*y.^5 + 147456.*x.^5.*y + 36462.*x.^2 - 100320.*x.^3 + 142400.*x.^4 - 101376.*x.^5 + 28672.*x.^6 + 12154.*y.^2 - 25080.*y.^3 + 28480.*y.^4 - 16896.*y.^5 + 4096.*y.^6 + 315))/45;
                gy = (256.*x.*(512640.*x.^2.*y.^2 - 6138.*y - 3069.*x - 675840.*x.^2.*y.^3 - 506880.*x.^3.*y.^2 + 307200.*x.^2.*y.^4 + 327680.*x.^3.*y.^3 + 184320.*x.^4.*y.^2 + 48616.*x.*y - 225720.*x.*y.^2 - 150480.*x.^2.*y + 455680.*x.*y.^3 + 227840.*x.^3.*y - 422400.*x.*y.^4 - 168960.*x.^4.*y + 147456.*x.*y.^5 + 49152.*x.^5.*y + 12154.*x.^2 - 25080.*x.^3 + 28480.*x.^4 - 16896.*x.^5 + 4096.*x.^6 + 36462.*y.^2 - 100320.*y.^3 + 142400.*y.^4 - 101376.*y.^5 + 28672.*y.^6 + 315))/45;
            case 26
                gx = -(256.*y.*(548.*x.*y - 15.*y - 578.*x - 5400.*x.^2.*y + 21760.*x.^3.*y - 38400.*x.^4.*y + 24576.*x.^5.*y + 6222.*x.^2 - 28960.*x.^3 + 65600.*x.^4 - 70656.*x.^5 + 28672.*x.^6 + 15))/45;
                gy = -(256.*x.*(x + 2.*y - 1).*(274.*x - 1800.*x.^2 + 5440.*x.^3 - 7680.*x.^4 + 4096.*x.^5 - 15))/45;
            case 27
                gx = -(256.*y.*(2.*x + y - 1).*(274.*y - 1800.*y.^2 + 5440.*y.^3 - 7680.*y.^4 + 4096.*y.^5 - 15))/45;
                gy = -(256.*x.*(548.*x.*y - 578.*y - 15.*x - 5400.*x.*y.^2 + 21760.*x.*y.^3 - 38400.*x.*y.^4 + 24576.*x.*y.^5 + 6222.*y.^2 - 28960.*y.^3 + 65600.*y.^4 - 70656.*y.^5 + 28672.*y.^6 + 15))/45;
            case 28
                gx = -(256.*y.*(8.*y - 1).*(106.*x + 3.*y - 100.*x.*y + 840.*x.^2.*y - 2560.*x.^3.*y + 2560.*x.^4.*y - 990.*x.^2 + 3680.*x.^3 - 5760.*x.^4 + 3072.*x.^5 - 3))/15;
                gy = -(256.*x.*(16.*x.*y - 18.*y - x + 24.*y.^2 + 1).*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3))/15;
            case 29
                gx = -(128.*y.*(32.*y.^2 - 12.*y + 1).*(88.*x.*y - 3.*y - 94.*x - 576.*x.^2.*y + 1024.*x.^3.*y + 708.*x.^2 - 1792.*x.^3 + 1280.*x.^4 + 3))/9;
                gy = -(128.*x.*(44.*x - 192.*x.^2 + 256.*x.^3 - 3).*(x + 26.*y - 24.*x.*y + 96.*x.*y.^2 - 132.*y.^2 + 128.*y.^3 - 1))/9;
            case 30
                gx = -(128.*y.*(44.*y - 192.*y.^2 + 256.*y.^3 - 3).*(26.*x + y - 24.*x.*y + 96.*x.^2.*y - 132.*x.^2 + 128.*x.^3 - 1))/9;
                gy = -(128.*x.*(32.*x.^2 - 12.*x + 1).*(88.*x.*y - 94.*y - 3.*x - 576.*x.*y.^2 + 1024.*x.*y.^3 + 708.*y.^2 - 1792.*y.^3 + 1280.*y.^4 + 3))/9;
            case 31
                gx = -(256.*y.*(16.*x.*y - y - 18.*x + 24.*x.^2 + 1).*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3))/15;
                gy = -(256.*x.*(8.*x - 1).*(3.*x + 106.*y - 100.*x.*y + 840.*x.*y.^2 - 2560.*x.*y.^3 + 2560.*x.*y.^4 - 990.*y.^2 + 3680.*y.^3 - 5760.*y.^4 + 3072.*y.^5 - 3))/15;
            case 32
                gx = (256.*y.*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3).*(32.*x.*y - 15.*y - 30.*x + 24.*x.^2 + 8.*y.^2 + 7))/15;
                gy = (256.*x.*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/15 + (256.*x.*y.*(16.*x + 16.*y - 15).*(280.*y.^2 - 50.*y - 640.*y.^3 + 512.*y.^4 + 3))/15 + (256.*x.*y.*(560.*y - 1920.*y.^2 + 2048.*y.^3 - 50).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/15;
            case 33
                gx = -(128.*y.*(44.*y - 192.*y.^2 + 256.*y.^3 - 3).*(146.*x + 73.*y - 336.*x.*y + 192.*x.*y.^2 + 288.*x.^2.*y - 252.*x.^2 + 128.*x.^3 - 84.*y.^2 + 32.*y.^3 - 21))/9;
                gy = -(128.*x.*(61056.*x.^2.*y.^2 - 2286.*y - 219.*x - 159744.*x.^2.*y.^3 - 18432.*x.^3.*y.^2 + 122880.*x.^2.*y.^4 + 32768.*x.^3.*y.^3 + 7432.*x.*y - 65088.*x.*y.^2 - 7968.*x.^2.*y + 220672.*x.*y.^3 + 2816.*x.^3.*y - 307200.*x.*y.^4 + 147456.*x.*y.^5 + 252.*x.^2 - 96.*x.^3 + 22488.*y.^2 - 92736.*y.^3 + 181120.*y.^4 - 165888.*y.^5 + 57344.*y.^6 + 63))/9;
            case 34
                gx = (128.*y.*(32.*y.^2 - 12.*y + 1).*(4608.*x.^2.*y.^2 - 533.*y - 1066.*x + 4016.*x.*y - 4992.*x.*y.^2 - 7488.*x.^2.*y + 2048.*x.*y.^3 + 4096.*x.^3.*y + 3012.*x.^2 - 3328.*x.^3 + 1280.*x.^4 + 1004.*y.^2 - 832.*y.^3 + 256.*y.^4 + 105))/9;
                gy = (128.*x.*(190848.*x.^2.*y.^2 - 3586.*y - 533.*x - 393216.*x.^2.*y.^3 - 116736.*x.^3.*y.^2 + 245760.*x.^2.*y.^4 + 131072.*x.^3.*y.^3 + 24576.*x.^4.*y.^2 + 16808.*x.*y - 130944.*x.*y.^2 - 29088.*x.^2.*y + 380928.*x.*y.^3 + 22016.*x.^3.*y - 460800.*x.*y.^4 - 6144.*x.^4.*y + 196608.*x.*y.^5 + 1004.*x.^2 - 832.*x.^3 + 256.*x.^4 + 32280.*y.^2 - 119744.*y.^3 + 211840.*y.^4 - 178176.*y.^5 + 57344.*y.^6 + 105))/9;
            case 35
                gx = -(256.*y.*(8.*y - 1).*(1486.*x + 743.*y - 34560.*x.^2.*y.^2 + 15360.*x.^2.*y.^3 + 20480.*x.^3.*y.^2 - 8280.*x.*y + 17040.*x.*y.^2 + 25560.*x.^2.*y - 15360.*x.*y.^3 - 30720.*x.^3.*y + 5120.*x.*y.^4 + 12800.*x.^4.*y - 6210.*x.^2 + 11360.*x.^3 - 9600.*x.^4 + 3072.*x.^5 - 2070.*y.^2 + 2840.*y.^3 - 1920.*y.^4 + 512.*y.^5 - 105))/15;
                gy = -(256.*x.*(239040.*x.^2.*y.^2 - 3166.*y - 743.*x - 389120.*x.^2.*y.^3 - 199680.*x.^3.*y.^2 + 204800.*x.^2.*y.^4 + 163840.*x.^3.*y.^3 + 61440.*x.^4.*y.^2 + 20168.*x.*y - 124920.*x.*y.^2 - 50160.*x.^2.*y + 303360.*x.*y.^3 + 60800.*x.^3.*y - 320000.*x.*y.^4 - 35840.*x.^4.*y + 122880.*x.*y.^5 + 8192.*x.^5.*y + 2070.*x.^2 - 2840.*x.^3 + 1920.*x.^4 - 512.*x.^5 + 24042.*y.^2 - 77600.*y.^3 + 123200.*y.^4 - 95232.*y.^5 + 28672.*y.^6 + 105))/15;
            case 36
                gx = -(256.*y.*(239040.*x.^2.*y.^2 - 743.*y - 3166.*x - 199680.*x.^2.*y.^3 - 389120.*x.^3.*y.^2 + 61440.*x.^2.*y.^4 + 163840.*x.^3.*y.^3 + 204800.*x.^4.*y.^2 + 20168.*x.*y - 50160.*x.*y.^2 - 124920.*x.^2.*y + 60800.*x.*y.^3 + 303360.*x.^3.*y - 35840.*x.*y.^4 - 320000.*x.^4.*y + 8192.*x.*y.^5 + 122880.*x.^5.*y + 24042.*x.^2 - 77600.*x.^3 + 123200.*x.^4 - 95232.*x.^5 + 28672.*x.^6 + 2070.*y.^2 - 2840.*y.^3 + 1920.*y.^4 - 512.*y.^5 + 105))/15;
                gy = -(256.*x.*(8.*x - 1).*(743.*x + 1486.*y - 34560.*x.^2.*y.^2 + 20480.*x.^2.*y.^3 + 15360.*x.^3.*y.^2 - 8280.*x.*y + 25560.*x.*y.^2 + 17040.*x.^2.*y - 30720.*x.*y.^3 - 15360.*x.^3.*y + 12800.*x.*y.^4 + 5120.*x.^4.*y - 2070.*x.^2 + 2840.*x.^3 - 1920.*x.^4 + 512.*x.^5 - 6210.*y.^2 + 11360.*y.^3 - 9600.*y.^4 + 3072.*y.^5 - 105))/15;
            case 37
                gx = (128.*y.*(190848.*x.^2.*y.^2 - 533.*y - 3586.*x - 116736.*x.^2.*y.^3 - 393216.*x.^3.*y.^2 + 24576.*x.^2.*y.^4 + 131072.*x.^3.*y.^3 + 245760.*x.^4.*y.^2 + 16808.*x.*y - 29088.*x.*y.^2 - 130944.*x.^2.*y + 22016.*x.*y.^3 + 380928.*x.^3.*y - 6144.*x.*y.^4 - 460800.*x.^4.*y + 196608.*x.^5.*y + 32280.*x.^2 - 119744.*x.^3 + 211840.*x.^4 - 178176.*x.^5 + 57344.*x.^6 + 1004.*y.^2 - 832.*y.^3 + 256.*y.^4 + 105))/9;
                gy = (128.*x.*(32.*x.^2 - 12.*x + 1).*(4608.*x.^2.*y.^2 - 1066.*y - 533.*x + 4016.*x.*y - 7488.*x.*y.^2 - 4992.*x.^2.*y + 4096.*x.*y.^3 + 2048.*x.^3.*y + 1004.*x.^2 - 832.*x.^3 + 256.*x.^4 + 3012.*y.^2 - 3328.*y.^3 + 1280.*y.^4 + 105))/9;
            case 38
                gx = -(128.*y.*(61056.*x.^2.*y.^2 - 219.*y - 2286.*x - 18432.*x.^2.*y.^3 - 159744.*x.^3.*y.^2 + 32768.*x.^3.*y.^3 + 122880.*x.^4.*y.^2 + 7432.*x.*y - 7968.*x.*y.^2 - 65088.*x.^2.*y + 2816.*x.*y.^3 + 220672.*x.^3.*y - 307200.*x.^4.*y + 147456.*x.^5.*y + 22488.*x.^2 - 92736.*x.^3 + 181120.*x.^4 - 165888.*x.^5 + 57344.*x.^6 + 252.*y.^2 - 96.*y.^3 + 63))/9;
                gy = -(128.*x.*(44.*x - 192.*x.^2 + 256.*x.^3 - 3).*(73.*x + 146.*y - 336.*x.*y + 288.*x.*y.^2 + 192.*x.^2.*y - 84.*x.^2 + 32.*x.^3 - 252.*y.^2 + 128.*y.^3 - 21))/9;
            case 39
                gx = (256.*y.*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/15 + (256.*x.*y.*(16.*x + 16.*y - 15).*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3))/15 + (256.*x.*y.*(560.*x - 1920.*x.^2 + 2048.*x.^3 - 50).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/15;
                gy = (256.*x.*(280.*x.^2 - 50.*x - 640.*x.^3 + 512.*x.^4 + 3).*(32.*x.*y - 30.*y - 15.*x + 8.*x.^2 + 24.*y.^2 + 7))/15;
            case 40
                gx = (32.*y.*(8.*y - 1).*(2746.*x + 533.*y - 64512.*x.^2.*y.^2 + 24576.*x.^2.*y.^3 + 49152.*x.^3.*y.^2 - 12544.*x.*y + 21056.*x.*y.^2 + 55680.*x.^2.*y - 15360.*x.*y.^3 - 83968.*x.^3.*y + 4096.*x.*y.^4 + 40960.*x.^4.*y - 15804.*x.^2 + 35456.*x.^3 - 34560.*x.^4 + 12288.*x.^5 - 1004.*y.^2 + 832.*y.^3 - 256.*y.^4 - 105))/3;
                gy = (32.*x.*(8.*x - 1).*(533.*x + 2746.*y - 64512.*x.^2.*y.^2 + 49152.*x.^2.*y.^3 + 24576.*x.^3.*y.^2 - 12544.*x.*y + 55680.*x.*y.^2 + 21056.*x.^2.*y - 83968.*x.*y.^3 - 15360.*x.^3.*y + 40960.*x.*y.^4 + 4096.*x.^4.*y - 1004.*x.^2 + 832.*x.^3 - 256.*x.^4 - 15804.*y.^2 + 35456.*y.^3 - 34560.*y.^4 + 12288.*y.^5 - 105))/3;
            case 41
                gx = (32.*y.*(8.*y - 1).*(44.*x - 192.*x.^2 + 256.*x.^3 - 3).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/3 + (32.*x.*y.*(8.*y - 1).*(16.*x + 16.*y - 15).*(44.*x - 192.*x.^2 + 256.*x.^3 - 3))/3 + (32.*x.*y.*(8.*y - 1).*(768.*x.^2 - 384.*x + 44).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/3;
                gy = (32.*x.*(44.*x - 192.*x.^2 + 256.*x.^3 - 3).*(15.*x + 142.*y - 272.*x.*y + 384.*x.*y.^2 + 128.*x.^2.*y - 8.*x.^2 - 384.*y.^2 + 256.*y.^3 - 7))/3;
            case 42
                gx = (32.*y.*(44.*y - 192.*y.^2 + 256.*y.^3 - 3).*(142.*x + 15.*y - 272.*x.*y + 128.*x.*y.^2 + 384.*x.^2.*y - 384.*x.^2 + 256.*x.^3 - 8.*y.^2 - 7))/3;
                gy = (32.*x.*(8.*x - 1).*(44.*y - 192.*y.^2 + 256.*y.^3 - 3).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/3 + (32.*x.*y.*(8.*x - 1).*(16.*x + 16.*y - 15).*(44.*y - 192.*y.^2 + 256.*y.^3 - 3))/3 + (32.*x.*y.*(8.*x - 1).*(768.*y.^2 - 384.*y + 44).*(16.*x.*y - 15.*y - 15.*x + 8.*x.^2 + 8.*y.^2 + 7))/3;
            case 43
                gx = (256.*y.*(32.*y.^2 - 12.*y + 1).*(768.*x.^2.*y.^2 - 15.*y - 198.*x + 392.*x.*y - 192.*x.*y.^2 - 2016.*x.^2.*y + 2048.*x.^3.*y + 1236.*x.^2 - 2304.*x.^3 + 1280.*x.^4 + 8.*y.^2 + 7))/9;
                gy = (256.*x.*(32.*x.^2 - 12.*x + 1).*(768.*x.^2.*y.^2 - 198.*y - 15.*x + 392.*x.*y - 2016.*x.*y.^2 - 192.*x.^2.*y + 2048.*x.*y.^3 + 8.*x.^2 + 1236.*y.^2 - 2304.*y.^3 + 1280.*y.^4 + 7))/9;
            case 44
                gx = -(256.*y.*(32.*y.^2 - 12.*y + 1).*(2304.*x.^2.*y.^2 - 73.*y - 482.*x + 1504.*x.*y - 1536.*x.*y.^2 - 4320.*x.^2.*y + 512.*x.*y.^3 + 3072.*x.^3.*y + 2004.*x.^2 - 2816.*x.^3 + 1280.*x.^4 + 84.*y.^2 - 32.*y.^3 + 21))/9;
                gy = -(256.*x.*(8.*x - 1).*(73.*x + 650.*y - 11520.*x.^2.*y.^2 + 12288.*x.^2.*y.^3 + 3072.*x.^3.*y.^2 - 2088.*x.*y + 13344.*x.*y.^2 + 2208.*x.^2.*y - 26112.*x.*y.^3 - 768.*x.^3.*y + 15360.*x.*y.^4 - 84.*x.^2 + 32.*x.^3 - 4896.*y.^2 + 13504.*y.^3 - 15360.*y.^4 + 6144.*y.^5 - 21))/9;
            case 45
                gx = -(256.*y.*(8.*y - 1).*(650.*x + 73.*y - 11520.*x.^2.*y.^2 + 3072.*x.^2.*y.^3 + 12288.*x.^3.*y.^2 - 2088.*x.*y + 2208.*x.*y.^2 + 13344.*x.^2.*y - 768.*x.*y.^3 - 26112.*x.^3.*y + 15360.*x.^4.*y - 4896.*x.^2 + 13504.*x.^3 - 15360.*x.^4 + 6144.*x.^5 - 84.*y.^2 + 32.*y.^3 - 21))/9;
                gy = -(256.*x.*(32.*x.^2 - 12.*x + 1).*(2304.*x.^2.*y.^2 - 482.*y - 73.*x + 1504.*x.*y - 4320.*x.*y.^2 - 1536.*x.^2.*y + 3072.*x.*y.^3 + 512.*x.^3.*y + 84.*x.^2 - 32.*x.^3 + 2004.*y.^2 - 2816.*y.^3 + 1280.*y.^4 + 21))/9;
        end
    case 9
        switch i
            case 1
                gx = (58635.*x)/112 + (58635.*y)/112 - (19768293.*x.^2.*y.^2)/64 + (13286025.*x.^2.*y.^3)/16 + (13286025.*x.^3.*y.^2)/16 - (77058945.*x.^2.*y.^4)/64 - (25686315.*x.^3.*y.^3)/16 - (77058945.*x.^4.*y.^2)/64 + (14348907.*x.^2.*y.^5)/16 + (23914845.*x.^3.*y.^4)/16 + (23914845.*x.^4.*y.^3)/16 + (14348907.*x.^5.*y.^2)/16 - (43046721.*x.^2.*y.^6)/160 - (43046721.*x.^3.*y.^5)/80 - (43046721.*x.^4.*y.^4)/64 - (43046721.*x.^5.*y.^3)/80 - (43046721.*x.^6.*y.^2)/160 - (122121.*x.*y)/14 + (1869885.*x.*y.^2)/32 + (1869885.*x.^2.*y)/32 - (6589431.*x.*y.^3)/32 - (6589431.*x.^3.*y)/32 + (13286025.*x.*y.^4)/32 + (13286025.*x.^4.*y)/32 - (15411789.*x.*y.^5)/32 - (15411789.*x.^5.*y)/32 + (4782969.*x.*y.^6)/16 + (4782969.*x.^6.*y)/16 - (43046721.*x.*y.^7)/560 - (43046721.*x.^7.*y)/560 - (122121.*x.^2)/28 + (623295.*x.^3)/32 - (6589431.*x.^4)/128 + (2657205.*x.^5)/32 - (5137263.*x.^6)/64 + (4782969.*x.^7)/112 - (43046721.*x.^8)/4480 - (122121.*y.^2)/28 + (623295.*y.^3)/32 - (6589431.*y.^4)/128 + (2657205.*y.^5)/32 - (5137263.*y.^6)/64 + (4782969.*y.^7)/112 - (43046721.*y.^8)/4480 - 7129/280;
                gy = (58635.*x)/112 + (58635.*y)/112 - (19768293.*x.^2.*y.^2)/64 + (13286025.*x.^2.*y.^3)/16 + (13286025.*x.^3.*y.^2)/16 - (77058945.*x.^2.*y.^4)/64 - (25686315.*x.^3.*y.^3)/16 - (77058945.*x.^4.*y.^2)/64 + (14348907.*x.^2.*y.^5)/16 + (23914845.*x.^3.*y.^4)/16 + (23914845.*x.^4.*y.^3)/16 + (14348907.*x.^5.*y.^2)/16 - (43046721.*x.^2.*y.^6)/160 - (43046721.*x.^3.*y.^5)/80 - (43046721.*x.^4.*y.^4)/64 - (43046721.*x.^5.*y.^3)/80 - (43046721.*x.^6.*y.^2)/160 - (122121.*x.*y)/14 + (1869885.*x.*y.^2)/32 + (1869885.*x.^2.*y)/32 - (6589431.*x.*y.^3)/32 - (6589431.*x.^3.*y)/32 + (13286025.*x.*y.^4)/32 + (13286025.*x.^4.*y)/32 - (15411789.*x.*y.^5)/32 - (15411789.*x.^5.*y)/32 + (4782969.*x.*y.^6)/16 + (4782969.*x.^6.*y)/16 - (43046721.*x.*y.^7)/560 - (43046721.*x.^7.*y)/560 - (122121.*x.^2)/28 + (623295.*x.^3)/32 - (6589431.*x.^4)/128 + (2657205.*x.^5)/32 - (5137263.*x.^6)/64 + (4782969.*x.^7)/112 - (43046721.*x.^8)/4480 - (122121.*y.^2)/28 + (623295.*y.^3)/32 - (6589431.*y.^4)/128 + (2657205.*y.^5)/32 - (5137263.*y.^6)/64 + (4782969.*y.^7)/112 - (43046721.*y.^8)/4480 - 7129/280;
            case 2
                gx = (797337.*x.^2)/1120 - (6849.*x)/140 - (194643.*x.^3)/40 + (2337903.*x.^4)/128 - (1594323.*x.^5)/40 + (16120377.*x.^6)/320 - (4782969.*x.^7)/140 + (43046721.*x.^8)/4480 + 1;
                gy = 0;
            case 3
                gx = 0;
                gy = (797337.*y.^2)/1120 - (6849.*y)/140 - (194643.*y.^3)/40 + (2337903.*y.^4)/128 - (1594323.*y.^5)/40 + (16120377.*y.^6)/320 - (4782969.*y.^7)/140 + (43046721.*y.^8)/4480 + 1;
            case 4
                gx = (81.*y.*(6534.*x - 88641.*x.^2 + 548289.*x.^3 - 1786050.*x.^4 + 3168963.*x.^5 - 2893401.*x.^6 + 1062882.*x.^7 - 140))/1120;
                gy = (81.*x.*(13068.*x - 118188.*x.^2 + 548289.*x.^3 - 1428840.*x.^4 + 2112642.*x.^5 - 1653372.*x.^6 + 531441.*x.^7 - 560))/4480;
            case 5
                gx = (81.*y.*(9.*y - 1).*(43848.*x.^2 - 3528.*x - 238140.*x.^3 + 637875.*x.^4 - 826686.*x.^5 + 413343.*x.^6 + 80))/1120;
                gy = (81.*x.*(18.*y - 1).*(14616.*x.^2 - 1764.*x - 59535.*x.^3 + 127575.*x.^4 - 137781.*x.^5 + 59049.*x.^6 + 80))/1120;
            case 6
                gx = (9.*y.*(81.*y.^2 - 27.*y + 2).*(1644.*x - 18225.*x.^2 + 82620.*x.^3 - 164025.*x.^4 + 118098.*x.^5 - 40))/160;
                gy = (9.*x.*(243.*y.^2 - 54.*y + 2).*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40))/160;
            case 7
                gx = (81.*y.*(33.*y - 162.*y.^2 + 243.*y.^3 - 2).*(2835.*x.^2 - 300.*x - 9720.*x.^3 + 10935.*x.^4 + 8))/320;
                gy = (81.*x.*(33.*y - 243.*y.^2 + 486.*y.^3 - 1).*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8))/160;
            case 8
                gx = (81.*y.*(33.*x - 243.*x.^2 + 486.*x.^3 - 1).*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8))/160;
                gy = (81.*x.*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(2835.*y.^2 - 300.*y - 9720.*y.^3 + 10935.*y.^4 + 8))/320;
            case 9
                gx = (9.*y.*(243.*x.^2 - 54.*x + 2).*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40))/160;
                gy = (9.*x.*(81.*x.^2 - 27.*x + 2).*(1644.*y - 18225.*y.^2 + 82620.*y.^3 - 164025.*y.^4 + 118098.*y.^5 - 40))/160;
            case 10
                gx = (81.*y.*(18.*x - 1).*(14616.*y.^2 - 1764.*y - 59535.*y.^3 + 127575.*y.^4 - 137781.*y.^5 + 59049.*y.^6 + 80))/1120;
                gy = (81.*x.*(9.*x - 1).*(43848.*y.^2 - 3528.*y - 238140.*y.^3 + 637875.*y.^4 - 826686.*y.^5 + 413343.*y.^6 + 80))/1120;
            case 11
                gx = (81.*y.*(13068.*y - 118188.*y.^2 + 548289.*y.^3 - 1428840.*y.^4 + 2112642.*y.^5 - 1653372.*y.^6 + 531441.*y.^7 - 560))/4480;
                gy = (81.*x.*(6534.*y - 88641.*y.^2 + 548289.*y.^3 - 1786050.*y.^4 + 3168963.*y.^5 - 2893401.*y.^6 + 1062882.*y.^7 - 140))/1120;
            case 12
                gx = -(81.*y.*(13068.*y - 118188.*y.^2 + 548289.*y.^3 - 1428840.*y.^4 + 2112642.*y.^5 - 1653372.*y.^6 + 531441.*y.^7 - 560))/4480;
                gy = (81.*x)/8 + (275967.*y)/560 - (264627.*x.*y)/560 + (1025703.*x.*y.^2)/160 - (6344487.*x.*y.^3)/160 + (2066715.*x.*y.^4)/16 - (36669429.*x.*y.^5)/160 + (33480783.*x.*y.^6)/160 - (43046721.*x.*y.^7)/560 - (3986901.*y.^2)/560 + (7712091.*y.^3)/160 - (22878207.*y.^4)/128 + (61470009.*y.^5)/160 - (152523567.*y.^6)/320 + (176969853.*y.^7)/560 - (387420489.*y.^8)/4480 - 81/8;
            case 13
                gx = (81.*y.*(18.*x + 18.*y - 17).*(14616.*y.^2 - 1764.*y - 59535.*y.^3 + 127575.*y.^4 - 137781.*y.^5 + 59049.*y.^6 + 80))/1120;
                gy = (81.*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8).*(14616.*y.^2 - 1764.*y - 59535.*y.^3 + 127575.*y.^4 - 137781.*y.^5 + 59049.*y.^6 + 80))/1120 + (81.*y.*(18.*x + 18.*y - 17).*(14616.*y.^2 - 1764.*y - 59535.*y.^3 + 127575.*y.^4 - 137781.*y.^5 + 59049.*y.^6 + 80))/1120 + (81.*y.*(29232.*y - 178605.*y.^2 + 510300.*y.^3 - 688905.*y.^4 + 354294.*y.^5 - 1764).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/1120;
            case 14
                gx = -(9.*y.*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40).*(486.*x.*y - 432.*y - 432.*x + 243.*x.^2 + 243.*y.^2 + 191))/160;
                gy = - (9.*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/160 - (9.*y.*(61965.*y.^2 - 12150.*y - 131220.*y.^3 + 98415.*y.^4 + 822).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/160 - (9.*y.*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40).*(486.*x.*y - 432.*y - 432.*x + 243.*x.^2 + 243.*y.^2 + 191))/160;
            case 15
                gx = (81.*y.*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8).*(2010.*x + 2010.*y - 4860.*x.*y + 2916.*x.*y.^2 + 2916.*x.^2.*y - 2430.*x.^2 + 972.*x.^3 - 2430.*y.^2 + 972.*y.^3 - 550))/320;
                gy = (322191027.*x.^2.*y.^2)/320 - (21465.*y)/2 - (4455.*x)/4 - 5019165.*x.^2.*y.^3 - (22143375.*x.^3.*y.^2)/32 + (767932245.*x.^2.*y.^4)/64 + (5845851.*x.^3.*y.^3)/2 + (11160261.*x.^4.*y.^2)/64 - (215233605.*x.^2.*y.^5)/16 - (167403915.*x.^3.*y.^4)/32 - (4782969.*x.^4.*y.^3)/8 + (903981141.*x.^2.*y.^6)/160 + (129140163.*x.^3.*y.^5)/40 + (43046721.*x.^4.*y.^4)/64 + (399249.*x.*y)/8 - (20428767.*x.*y.^2)/32 - (1378539.*x.^2.*y)/16 + (146133153.*x.*y.^3)/40 + (2617839.*x.^3.*y)/40 - (172718325.*x.*y.^4)/16 - (295245.*x.^4.*y)/16 + (272629233.*x.*y.^5)/16 - (435250179.*x.*y.^6)/32 + (43046721.*x.*y.^7)/10 + (16281.*x.^2)/8 - (6561.*x.^3)/4 + (19683.*x.^4)/40 + (2386017.*y.^2)/16 - (3844017.*y.^3)/4 + (215023653.*y.^4)/64 - (54029835.*y.^5)/8 + (249245829.*y.^6)/32 - 4782969.*y.^7 + (387420489.*y.^8)/320 + 1134/5;
            case 16
                gx = -(81.*y.*(33.*y - 162.*y.^2 + 243.*y.^3 - 2).*(65610.*x.^2.*y.^2 - 19950.*y - 19950.*x + 78570.*x.*y - 102060.*x.*y.^2 - 102060.*x.^2.*y + 43740.*x.*y.^3 + 43740.*x.^3.*y + 39285.*x.^2 - 34020.*x.^3 + 10935.*x.^4 + 39285.*y.^2 - 34020.*y.^3 + 10935.*y.^4 + 3758))/320;
                gy = (152199.*x)/80 + (526419.*y)/40 - (146500569.*x.^2.*y.^2)/64 + (170356365.*x.^2.*y.^3)/16 + (79893297.*x.^3.*y.^2)/32 - (1501320825.*x.^2.*y.^4)/64 - (152523567.*x.^3.*y.^3)/16 - (84499119.*x.^4.*y.^2)/64 + (387420489.*x.^2.*y.^5)/16 + (119574225.*x.^3.*y.^4)/8 + (62178597.*x.^4.*y.^3)/16 + (43046721.*x.^5.*y.^2)/160 - (301327047.*x.^2.*y.^6)/32 - (129140163.*x.^3.*y.^5)/16 - (215233605.*x.^4.*y.^4)/64 - (43046721.*x.^5.*y.^3)/80 - (6638517.*x.*y)/80 + (81752247.*x.*y.^2)/80 + (6605469.*x.^2.*y)/32 - (446272659.*x.*y.^3)/80 - (8102835.*x.^3.*y)/32 + (500440275.*x.*y.^4)/32 + (4901067.*x.^4.*y)/32 - (374665905.*x.*y.^5)/16 - (5845851.*x.^5.*y)/160 + (569173311.*x.*y.^6)/32 - (43046721.*x.*y.^7)/8 - (161595.*x.^2)/32 + (212139.*x.^3)/32 - (137781.*x.^4)/32 + (177147.*x.^5)/160 - (14257053.*y.^2)/80 + (89119521.*y.^3)/80 - (241241409.*y.^4)/64 + (586888011.*y.^5)/80 - (1313190711.*y.^6)/160 + (196101729.*y.^7)/40 - (387420489.*y.^8)/320 - 567/2;
            case 17
                gx = (9.*y.*(81.*y.^2 - 27.*y + 2).*(147444.*x + 147444.*y - 2558790.*x.^2.*y.^2 + 1180980.*x.^2.*y.^3 + 1180980.*x.^3.*y.^2 - 852930.*x.*y + 1822500.*x.*y.^2 + 1822500.*x.^2.*y - 1705860.*x.*y.^3 - 1705860.*x.^3.*y + 590490.*x.*y.^4 + 590490.*x.^4.*y - 426465.*x.^2 + 607500.*x.^3 - 426465.*x.^4 + 118098.*x.^5 - 426465.*y.^2 + 607500.*y.^3 - 426465.*y.^4 + 118098.*y.^5 - 20072))/160;
                gy = (521330499.*x.^2.*y.^2)/160 - (56601.*y)/5 - (22581.*x)/10 - (109535895.*x.^2.*y.^3)/8 - (159963741.*x.^3.*y.^2)/32 + (438438825.*x.^2.*y.^4)/16 + (65721537.*x.^3.*y.^3)/4 + (16474671.*x.^4.*y.^2)/4 - (416118303.*x.^2.*y.^5)/16 - (358722675.*x.^3.*y.^4)/16 - 9565938.*x.^4.*y.^3 - (272629233.*x.^5.*y.^2)/160 + (301327047.*x.^2.*y.^6)/32 + (43046721.*x.^3.*y.^5)/4 + (215233605.*x.^4.*y.^4)/32 + (43046721.*x.^5.*y.^3)/20 + (43046721.*x.^6.*y.^2)/160 + (470718.*x.*y)/5 - (17441325.*x.*y.^2)/16 - (1599426.*x.^2.*y)/5 + (222052671.*x.*y.^3)/40 + (9095733.*x.^3.*y)/16 - (466191855.*x.*y.^4)/32 - (8916399.*x.^4.*y)/16 + (1645872777.*x.*y.^5)/80 + (22851963.*x.^5.*y)/80 - (2377135593.*x.*y.^6)/160 - (4782969.*x.^6.*y)/80 + (43046721.*x.*y.^7)/10 + (331749.*x.^2)/40 - (255879.*x.^3)/16 + (273375.*x.^4)/16 - (767637.*x.^5)/80 + (177147.*x.^6)/80 + (5878089.*y.^2)/40 - (8776431.*y.^3)/10 + (91020753.*y.^4)/32 - (213107841.*y.^5)/40 + (115322697.*y.^6)/20 - (33480783.*y.^7)/10 + (129140163.*y.^8)/160 + 252;
            case 18
                gx = -(81.*y.*(9.*y - 1).*(16227540.*x.^2.*y.^2 - 267876.*y - 267876.*x - 16533720.*x.^2.*y.^3 - 16533720.*x.^3.*y.^2 + 6200145.*x.^2.*y.^4 + 8266860.*x.^3.*y.^3 + 6200145.*x.^4.*y.^2 + 2179926.*x.*y - 6940080.*x.*y.^2 - 6940080.*x.^2.*y + 10818360.*x.*y.^3 + 10818360.*x.^3.*y - 8266860.*x.*y.^4 - 8266860.*x.^4.*y + 2480058.*x.*y.^5 + 2480058.*x.^5.*y + 1089963.*x.^2 - 2313360.*x.^3 + 2704590.*x.^4 - 1653372.*x.^5 + 413343.*x.^6 + 1089963.*y.^2 - 2313360.*y.^3 + 2704590.*y.^4 - 1653372.*y.^5 + 413343.*y.^6 + 26792))/1120;
                gy = (271269.*x)/140 + (475389.*y)/70 - (460995543.*x.^2.*y.^2)/160 + (21198591.*x.^2.*y.^3)/2 + (45526779.*x.^3.*y.^2)/8 - (305578575.*x.^2.*y.^4)/16 - (31355019.*x.^3.*y.^3)/2 - (49424013.*x.^4.*y.^2)/8 + (1334448351.*x.^2.*y.^5)/80 + (597871125.*x.^3.*y.^4)/32 + (90876411.*x.^4.*y.^3)/8 + (559607373.*x.^5.*y.^2)/160 - (903981141.*x.^2.*y.^6)/160 - (129140163.*x.^3.*y.^5)/16 - (215233605.*x.^4.*y.^4)/32 - (129140163.*x.^5.*y.^3)/40 - (129140163.*x.^6.*y.^2)/160 - (10307331.*x.*y)/140 + (121529403.*x.*y.^2)/160 + (5312223.*x.^2.*y)/16 - (140280741.*x.*y.^3)/40 - (64606167.*x.^3.*y)/80 + (136107945.*x.*y.^4)/16 + (9152595.*x.^4.*y)/8 - (451193409.*x.*y.^5)/40 - (37732311.*x.^5.*y)/40 + (1238788971.*x.*y.^6)/160 + (33480783.*x.^6.*y)/80 - (43046721.*x.*y.^7)/20 - (43046721.*x.^7.*y)/560 - (774927.*x.^2)/80 + (4204143.*x.^3)/160 - (334611.*x.^4)/8 + (3129597.*x.^5)/80 - (1594323.*x.^6)/80 + (4782969.*x.^7)/1120 - (45570519.*y.^2)/560 + (18152829.*y.^3)/40 - (44529507.*y.^4)/32 + (99733761.*y.^5)/40 - (26040609.*y.^6)/10 + (205667667.*y.^7)/140 - (387420489.*y.^8)/1120 - 162;
            case 19
                gx = (81.*y.*(1018008.*x + 1018008.*y - 188606880.*x.^2.*y.^2 + 325163160.*x.^2.*y.^3 + 325163160.*x.^3.*y.^2 - 272806380.*x.^2.*y.^4 - 363741840.*x.^3.*y.^3 - 272806380.*x.^4.*y.^2 + 89282088.*x.^2.*y.^5 + 148803480.*x.^3.*y.^4 + 148803480.*x.^4.*y.^3 + 89282088.*x.^5.*y.^2 - 11592504.*x.*y + 53118828.*x.*y.^2 + 53118828.*x.^2.*y - 125737920.*x.*y.^3 - 125737920.*x.^3.*y + 162581580.*x.*y.^4 + 162581580.*x.^4.*y - 109122552.*x.*y.^5 - 109122552.*x.^5.*y + 29760696.*x.*y.^6 + 29760696.*x.^6.*y - 5796252.*x.^2 + 17706276.*x.^3 - 31434480.*x.^4 + 32516316.*x.^5 - 18187092.*x.^6 + 4251528.*x.^7 - 5796252.*y.^2 + 17706276.*y.^3 - 31434480.*y.^4 + 32516316.*y.^5 - 18187092.*y.^6 + 4251528.*y.^7 - 73744))/4480;
                gy = (460995543.*x.^2.*y.^2)/320 - (373329.*y)/140 - (373329.*x)/280 - 4546773.*x.^2.*y.^3 - (13640319.*x.^3.*y.^2)/4 + (470325285.*x.^2.*y.^4)/64 + (31355019.*x.^3.*y.^3)/4 + (282195171.*x.^4.*y.^2)/64 - (473513931.*x.^2.*y.^5)/80 - (263063295.*x.^3.*y.^4)/32 - (52612659.*x.^4.*y.^3)/8 - (473513931.*x.^5.*y.^2)/160 + (301327047.*x.^2.*y.^6)/160 + (129140163.*x.^3.*y.^5)/40 + (215233605.*x.^4.*y.^4)/64 + (43046721.*x.^5.*y.^3)/20 + (129140163.*x.^6.*y.^2)/160 + (10307331.*x.*y)/280 - (50303187.*x.*y.^2)/160 - (16767729.*x.^2.*y)/80 + (51221727.*x.*y.^3)/40 + (51221727.*x.^3.*y)/80 - (22733865.*x.*y.^4)/8 - (4546773.*x.^4.*y)/4 + (282195171.*x.*y.^5)/80 + (94065057.*x.^5.*y)/80 - (368288613.*x.*y.^6)/160 - (52612659.*x.^6.*y)/80 + (43046721.*x.*y.^7)/70 + (43046721.*x.^7.*y)/280 + (10307331.*x.^2)/1120 - (5589243.*x.^3)/160 + (51221727.*x.^4)/640 - (4546773.*x.^5)/40 + (31355019.*x.^6)/320 - (52612659.*x.^7)/1120 + (43046721.*x.^8)/4480 + (30921993.*y.^2)/1120 - (5589243.*y.^3)/40 + (51221727.*y.^4)/128 - (13640319.*y.^5)/20 + (219485133.*y.^6)/320 - (52612659.*y.^7)/140 + (387420489.*y.^8)/4480 + 81;
            case 20
                gx = (460995543.*x.^2.*y.^2)/320 - (373329.*y)/280 - (373329.*x)/140 - (13640319.*x.^2.*y.^3)/4 - 4546773.*x.^3.*y.^2 + (282195171.*x.^2.*y.^4)/64 + (31355019.*x.^3.*y.^3)/4 + (470325285.*x.^4.*y.^2)/64 - (473513931.*x.^2.*y.^5)/160 - (52612659.*x.^3.*y.^4)/8 - (263063295.*x.^4.*y.^3)/32 - (473513931.*x.^5.*y.^2)/80 + (129140163.*x.^2.*y.^6)/160 + (43046721.*x.^3.*y.^5)/20 + (215233605.*x.^4.*y.^4)/64 + (129140163.*x.^5.*y.^3)/40 + (301327047.*x.^6.*y.^2)/160 + (10307331.*x.*y)/280 - (16767729.*x.*y.^2)/80 - (50303187.*x.^2.*y)/160 + (51221727.*x.*y.^3)/80 + (51221727.*x.^3.*y)/40 - (4546773.*x.*y.^4)/4 - (22733865.*x.^4.*y)/8 + (94065057.*x.*y.^5)/80 + (282195171.*x.^5.*y)/80 - (52612659.*x.*y.^6)/80 - (368288613.*x.^6.*y)/160 + (43046721.*x.*y.^7)/280 + (43046721.*x.^7.*y)/70 + (30921993.*x.^2)/1120 - (5589243.*x.^3)/40 + (51221727.*x.^4)/128 - (13640319.*x.^5)/20 + (219485133.*x.^6)/320 - (52612659.*x.^7)/140 + (387420489.*x.^8)/4480 + (10307331.*y.^2)/1120 - (5589243.*y.^3)/160 + (51221727.*y.^4)/640 - (4546773.*y.^5)/40 + (31355019.*y.^6)/320 - (52612659.*y.^7)/1120 + (43046721.*y.^8)/4480 + 81;
                gy = (81.*x.*(1018008.*x + 1018008.*y - 188606880.*x.^2.*y.^2 + 325163160.*x.^2.*y.^3 + 325163160.*x.^3.*y.^2 - 272806380.*x.^2.*y.^4 - 363741840.*x.^3.*y.^3 - 272806380.*x.^4.*y.^2 + 89282088.*x.^2.*y.^5 + 148803480.*x.^3.*y.^4 + 148803480.*x.^4.*y.^3 + 89282088.*x.^5.*y.^2 - 11592504.*x.*y + 53118828.*x.*y.^2 + 53118828.*x.^2.*y - 125737920.*x.*y.^3 - 125737920.*x.^3.*y + 162581580.*x.*y.^4 + 162581580.*x.^4.*y - 109122552.*x.*y.^5 - 109122552.*x.^5.*y + 29760696.*x.*y.^6 + 29760696.*x.^6.*y - 5796252.*x.^2 + 17706276.*x.^3 - 31434480.*x.^4 + 32516316.*x.^5 - 18187092.*x.^6 + 4251528.*x.^7 - 5796252.*y.^2 + 17706276.*y.^3 - 31434480.*y.^4 + 32516316.*y.^5 - 18187092.*y.^6 + 4251528.*y.^7 - 73744))/4480;
            case 21
                gx = (475389.*x)/70 + (271269.*y)/140 - (460995543.*x.^2.*y.^2)/160 + (45526779.*x.^2.*y.^3)/8 + (21198591.*x.^3.*y.^2)/2 - (49424013.*x.^2.*y.^4)/8 - (31355019.*x.^3.*y.^3)/2 - (305578575.*x.^4.*y.^2)/16 + (559607373.*x.^2.*y.^5)/160 + (90876411.*x.^3.*y.^4)/8 + (597871125.*x.^4.*y.^3)/32 + (1334448351.*x.^5.*y.^2)/80 - (129140163.*x.^2.*y.^6)/160 - (129140163.*x.^3.*y.^5)/40 - (215233605.*x.^4.*y.^4)/32 - (129140163.*x.^5.*y.^3)/16 - (903981141.*x.^6.*y.^2)/160 - (10307331.*x.*y)/140 + (5312223.*x.*y.^2)/16 + (121529403.*x.^2.*y)/160 - (64606167.*x.*y.^3)/80 - (140280741.*x.^3.*y)/40 + (9152595.*x.*y.^4)/8 + (136107945.*x.^4.*y)/16 - (37732311.*x.*y.^5)/40 - (451193409.*x.^5.*y)/40 + (33480783.*x.*y.^6)/80 + (1238788971.*x.^6.*y)/160 - (43046721.*x.*y.^7)/560 - (43046721.*x.^7.*y)/20 - (45570519.*x.^2)/560 + (18152829.*x.^3)/40 - (44529507.*x.^4)/32 + (99733761.*x.^5)/40 - (26040609.*x.^6)/10 + (205667667.*x.^7)/140 - (387420489.*x.^8)/1120 - (774927.*y.^2)/80 + (4204143.*y.^3)/160 - (334611.*y.^4)/8 + (3129597.*y.^5)/80 - (1594323.*y.^6)/80 + (4782969.*y.^7)/1120 - 162;
                gy = -(81.*x.*(9.*x - 1).*(16227540.*x.^2.*y.^2 - 267876.*y - 267876.*x - 16533720.*x.^2.*y.^3 - 16533720.*x.^3.*y.^2 + 6200145.*x.^2.*y.^4 + 8266860.*x.^3.*y.^3 + 6200145.*x.^4.*y.^2 + 2179926.*x.*y - 6940080.*x.*y.^2 - 6940080.*x.^2.*y + 10818360.*x.*y.^3 + 10818360.*x.^3.*y - 8266860.*x.*y.^4 - 8266860.*x.^4.*y + 2480058.*x.*y.^5 + 2480058.*x.^5.*y + 1089963.*x.^2 - 2313360.*x.^3 + 2704590.*x.^4 - 1653372.*x.^5 + 413343.*x.^6 + 1089963.*y.^2 - 2313360.*y.^3 + 2704590.*y.^4 - 1653372.*y.^5 + 413343.*y.^6 + 26792))/1120;
            case 22
                gx = (521330499.*x.^2.*y.^2)/160 - (22581.*y)/10 - (56601.*x)/5 - (159963741.*x.^2.*y.^3)/32 - (109535895.*x.^3.*y.^2)/8 + (16474671.*x.^2.*y.^4)/4 + (65721537.*x.^3.*y.^3)/4 + (438438825.*x.^4.*y.^2)/16 - (272629233.*x.^2.*y.^5)/160 - 9565938.*x.^3.*y.^4 - (358722675.*x.^4.*y.^3)/16 - (416118303.*x.^5.*y.^2)/16 + (43046721.*x.^2.*y.^6)/160 + (43046721.*x.^3.*y.^5)/20 + (215233605.*x.^4.*y.^4)/32 + (43046721.*x.^5.*y.^3)/4 + (301327047.*x.^6.*y.^2)/32 + (470718.*x.*y)/5 - (1599426.*x.*y.^2)/5 - (17441325.*x.^2.*y)/16 + (9095733.*x.*y.^3)/16 + (222052671.*x.^3.*y)/40 - (8916399.*x.*y.^4)/16 - (466191855.*x.^4.*y)/32 + (22851963.*x.*y.^5)/80 + (1645872777.*x.^5.*y)/80 - (4782969.*x.*y.^6)/80 - (2377135593.*x.^6.*y)/160 + (43046721.*x.^7.*y)/10 + (5878089.*x.^2)/40 - (8776431.*x.^3)/10 + (91020753.*x.^4)/32 - (213107841.*x.^5)/40 + (115322697.*x.^6)/20 - (33480783.*x.^7)/10 + (129140163.*x.^8)/160 + (331749.*y.^2)/40 - (255879.*y.^3)/16 + (273375.*y.^4)/16 - (767637.*y.^5)/80 + (177147.*y.^6)/80 + 252;
                gy = (9.*x.*(81.*x.^2 - 27.*x + 2).*(147444.*x + 147444.*y - 2558790.*x.^2.*y.^2 + 1180980.*x.^2.*y.^3 + 1180980.*x.^3.*y.^2 - 852930.*x.*y + 1822500.*x.*y.^2 + 1822500.*x.^2.*y - 1705860.*x.*y.^3 - 1705860.*x.^3.*y + 590490.*x.*y.^4 + 590490.*x.^4.*y - 426465.*x.^2 + 607500.*x.^3 - 426465.*x.^4 + 118098.*x.^5 - 426465.*y.^2 + 607500.*y.^3 - 426465.*y.^4 + 118098.*y.^5 - 20072))/160;
            case 23
                gx = (526419.*x)/40 + (152199.*y)/80 - (146500569.*x.^2.*y.^2)/64 + (79893297.*x.^2.*y.^3)/32 + (170356365.*x.^3.*y.^2)/16 - (84499119.*x.^2.*y.^4)/64 - (152523567.*x.^3.*y.^3)/16 - (1501320825.*x.^4.*y.^2)/64 + (43046721.*x.^2.*y.^5)/160 + (62178597.*x.^3.*y.^4)/16 + (119574225.*x.^4.*y.^3)/8 + (387420489.*x.^5.*y.^2)/16 - (43046721.*x.^3.*y.^5)/80 - (215233605.*x.^4.*y.^4)/64 - (129140163.*x.^5.*y.^3)/16 - (301327047.*x.^6.*y.^2)/32 - (6638517.*x.*y)/80 + (6605469.*x.*y.^2)/32 + (81752247.*x.^2.*y)/80 - (8102835.*x.*y.^3)/32 - (446272659.*x.^3.*y)/80 + (4901067.*x.*y.^4)/32 + (500440275.*x.^4.*y)/32 - (5845851.*x.*y.^5)/160 - (374665905.*x.^5.*y)/16 + (569173311.*x.^6.*y)/32 - (43046721.*x.^7.*y)/8 - (14257053.*x.^2)/80 + (89119521.*x.^3)/80 - (241241409.*x.^4)/64 + (586888011.*x.^5)/80 - (1313190711.*x.^6)/160 + (196101729.*x.^7)/40 - (387420489.*x.^8)/320 - (161595.*y.^2)/32 + (212139.*y.^3)/32 - (137781.*y.^4)/32 + (177147.*y.^5)/160 - 567/2;
                gy = -(81.*x.*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(65610.*x.^2.*y.^2 - 19950.*y - 19950.*x + 78570.*x.*y - 102060.*x.*y.^2 - 102060.*x.^2.*y + 43740.*x.*y.^3 + 43740.*x.^3.*y + 39285.*x.^2 - 34020.*x.^3 + 10935.*x.^4 + 39285.*y.^2 - 34020.*y.^3 + 10935.*y.^4 + 3758))/320;
            case 24
                gx = (322191027.*x.^2.*y.^2)/320 - (4455.*y)/4 - (21465.*x)/2 - (22143375.*x.^2.*y.^3)/32 - 5019165.*x.^3.*y.^2 + (11160261.*x.^2.*y.^4)/64 + (5845851.*x.^3.*y.^3)/2 + (767932245.*x.^4.*y.^2)/64 - (4782969.*x.^3.*y.^4)/8 - (167403915.*x.^4.*y.^3)/32 - (215233605.*x.^5.*y.^2)/16 + (43046721.*x.^4.*y.^4)/64 + (129140163.*x.^5.*y.^3)/40 + (903981141.*x.^6.*y.^2)/160 + (399249.*x.*y)/8 - (1378539.*x.*y.^2)/16 - (20428767.*x.^2.*y)/32 + (2617839.*x.*y.^3)/40 + (146133153.*x.^3.*y)/40 - (295245.*x.*y.^4)/16 - (172718325.*x.^4.*y)/16 + (272629233.*x.^5.*y)/16 - (435250179.*x.^6.*y)/32 + (43046721.*x.^7.*y)/10 + (2386017.*x.^2)/16 - (3844017.*x.^3)/4 + (215023653.*x.^4)/64 - (54029835.*x.^5)/8 + (249245829.*x.^6)/32 - 4782969.*x.^7 + (387420489.*x.^8)/320 + (16281.*y.^2)/8 - (6561.*y.^3)/4 + (19683.*y.^4)/40 + 1134/5;
                gy = (81.*x.*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8).*(2010.*x + 2010.*y - 4860.*x.*y + 2916.*x.*y.^2 + 2916.*x.^2.*y - 2430.*x.^2 + 972.*x.^3 - 2430.*y.^2 + 972.*y.^3 - 550))/320;
            case 25
                gx = - (9.*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/160 - (9.*x.*(61965.*x.^2 - 12150.*x - 131220.*x.^3 + 98415.*x.^4 + 822).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/160 - (9.*x.*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40).*(486.*x.*y - 432.*y - 432.*x + 243.*x.^2 + 243.*y.^2 + 191))/160;
                gy = -(9.*x.*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40).*(486.*x.*y - 432.*y - 432.*x + 243.*x.^2 + 243.*y.^2 + 191))/160;
            case 26
                gx = (81.*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8).*(14616.*x.^2 - 1764.*x - 59535.*x.^3 + 127575.*x.^4 - 137781.*x.^5 + 59049.*x.^6 + 80))/1120 + (81.*x.*(18.*x + 18.*y - 17).*(14616.*x.^2 - 1764.*x - 59535.*x.^3 + 127575.*x.^4 - 137781.*x.^5 + 59049.*x.^6 + 80))/1120 + (81.*x.*(29232.*x - 178605.*x.^2 + 510300.*x.^3 - 688905.*x.^4 + 354294.*x.^5 - 1764).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/1120;
                gy = (81.*x.*(18.*x + 18.*y - 17).*(14616.*x.^2 - 1764.*x - 59535.*x.^3 + 127575.*x.^4 - 137781.*x.^5 + 59049.*x.^6 + 80))/1120;
            case 27
                gx = (275967.*x)/560 + (81.*y)/8 - (264627.*x.*y)/560 + (1025703.*x.^2.*y)/160 - (6344487.*x.^3.*y)/160 + (2066715.*x.^4.*y)/16 - (36669429.*x.^5.*y)/160 + (33480783.*x.^6.*y)/160 - (43046721.*x.^7.*y)/560 - (3986901.*x.^2)/560 + (7712091.*x.^3)/160 - (22878207.*x.^4)/128 + (61470009.*x.^5)/160 - (152523567.*x.^6)/320 + (176969853.*x.^7)/560 - (387420489.*x.^8)/4480 - 81/8;
                gy = -(81.*x.*(13068.*x - 118188.*x.^2 + 548289.*x.^3 - 1428840.*x.^4 + 2112642.*x.^5 - 1653372.*x.^6 + 531441.*x.^7 - 560))/4480;
            case 28
                gx = -(729.*y.*(53584.*x + 26792.*y - 10410120.*x.^2.*y.^2 + 16227540.*x.^2.*y.^3 + 21636720.*x.^3.*y.^2 - 12400290.*x.^2.*y.^4 - 22044960.*x.^3.*y.^3 - 20667150.*x.^4.*y.^2 + 3720087.*x.^2.*y.^5 + 8266860.*x.^3.*y.^4 + 10333575.*x.^4.*y.^3 + 7440174.*x.^5.*y.^2 - 535752.*x.*y + 2179926.*x.*y.^2 + 3269889.*x.^2.*y - 4626720.*x.*y.^3 - 9253440.*x.^3.*y + 5409180.*x.*y.^4 + 13522950.*x.^4.*y - 3306744.*x.*y.^5 - 9920232.*x.^5.*y + 826686.*x.*y.^6 + 2893401.*x.^6.*y - 401814.*x.^2 + 1453284.*x.^3 - 2891700.*x.^4 + 3245508.*x.^5 - 1928934.*x.^6 + 472392.*x.^7 - 133938.*y.^2 + 363321.*y.^3 - 578340.*y.^4 + 540918.*y.^5 - 275562.*y.^6 + 59049.*y.^7 - 2240))/560;
                gy = -(729.*x.*(26792.*x + 53584.*y - 10410120.*x.^2.*y.^2 + 21636720.*x.^2.*y.^3 + 16227540.*x.^3.*y.^2 - 20667150.*x.^2.*y.^4 - 22044960.*x.^3.*y.^3 - 12400290.*x.^4.*y.^2 + 7440174.*x.^2.*y.^5 + 10333575.*x.^3.*y.^4 + 8266860.*x.^4.*y.^3 + 3720087.*x.^5.*y.^2 - 535752.*x.*y + 3269889.*x.*y.^2 + 2179926.*x.^2.*y - 9253440.*x.*y.^3 - 4626720.*x.^3.*y + 13522950.*x.*y.^4 + 5409180.*x.^4.*y - 9920232.*x.*y.^5 - 3306744.*x.^5.*y + 2893401.*x.*y.^6 + 826686.*x.^6.*y - 133938.*x.^2 + 363321.*x.^3 - 578340.*x.^4 + 540918.*x.^5 - 275562.*x.^6 + 59049.*x.^7 - 401814.*y.^2 + 1453284.*y.^3 - 2891700.*y.^4 + 3245508.*y.^5 - 1928934.*y.^6 + 472392.*y.^7 - 2240))/560;
            case 29
                gx = -(729.*y.*(3688.*x + 80.*y - 3528.*x.*y + 43848.*x.^2.*y - 238140.*x.^3.*y + 637875.*x.^4.*y - 826686.*x.^5.*y + 413343.*x.^6.*y - 49140.*x.^2 + 296604.*x.^3 - 935550.*x.^4 + 1592136.*x.^5 - 1377810.*x.^6 + 472392.*x.^7 - 80))/560;
                gy = -(729.*x.*(x + 2.*y - 1).*(14616.*x.^2 - 1764.*x - 59535.*x.^3 + 127575.*x.^4 - 137781.*x.^5 + 59049.*x.^6 + 80))/560;
            case 30
                gx = -(729.*y.*(2.*x + y - 1).*(14616.*y.^2 - 1764.*y - 59535.*y.^3 + 127575.*y.^4 - 137781.*y.^5 + 59049.*y.^6 + 80))/560;
                gy = -(729.*x.*(80.*x + 3688.*y - 3528.*x.*y + 43848.*x.*y.^2 - 238140.*x.*y.^3 + 637875.*x.*y.^4 - 826686.*x.*y.^5 + 413343.*x.*y.^6 - 49140.*y.^2 + 296604.*y.^3 - 935550.*y.^4 + 1592136.*y.^5 - 1377810.*y.^6 + 472392.*y.^7 - 80))/560;
            case 31
                gx = -(243.*y.*(9.*y - 1).*(1644.*x.*y - 40.*y - 1724.*x - 18225.*x.^2.*y + 82620.*x.^3.*y - 164025.*x.^4.*y + 118098.*x.^5.*y + 20691.*x.^2 - 106920.*x.^3 + 267300.*x.^4 - 314928.*x.^5 + 137781.*x.^6 + 40))/160;
                gy = -(243.*x.*(18.*x.*y - 20.*y - x + 27.*y.^2 + 1).*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40))/160;
            case 32
                gx = -(243.*y.*(81.*y.^2 - 27.*y + 2).*(316.*x + 8.*y - 300.*x.*y + 2835.*x.^2.*y - 9720.*x.^3.*y + 10935.*x.^4.*y - 3285.*x.^2 + 13500.*x.^3 - 23085.*x.^4 + 13122.*x.^5 - 8))/80;
                gy = -(243.*x.*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8).*(2.*x + 58.*y - 54.*x.*y + 243.*x.*y.^2 - 324.*y.^2 + 324.*y.^3 - 2))/80;
            case 33
                gx = -(729.*y.*(33.*y - 162.*y.^2 + 243.*y.^3 - 2).*(66.*x.*y - 2.*y - 70.*x - 486.*x.^2.*y + 972.*x.^3.*y + 585.*x.^2 - 1620.*x.^3 + 1215.*x.^4 + 2))/64;
                gy = -(729.*x.*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(66.*x.*y - 70.*y - 2.*x - 486.*x.*y.^2 + 972.*x.*y.^3 + 585.*y.^2 - 1620.*y.^3 + 1215.*y.^4 + 2))/64;
            case 34
                gx = -(243.*y.*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8).*(58.*x + 2.*y - 54.*x.*y + 243.*x.^2.*y - 324.*x.^2 + 324.*x.^3 - 2))/80;
                gy = -(243.*x.*(81.*x.^2 - 27.*x + 2).*(8.*x + 316.*y - 300.*x.*y + 2835.*x.*y.^2 - 9720.*x.*y.^3 + 10935.*x.*y.^4 - 3285.*y.^2 + 13500.*y.^3 - 23085.*y.^4 + 13122.*y.^5 - 8))/80;
            case 35
                gx = -(243.*y.*(18.*x.*y - y - 20.*x + 27.*x.^2 + 1).*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40))/160;
                gy = -(243.*x.*(9.*x - 1).*(1644.*x.*y - 1724.*y - 40.*x - 18225.*x.*y.^2 + 82620.*x.*y.^3 - 164025.*x.*y.^4 + 118098.*x.*y.^5 + 20691.*y.^2 - 106920.*y.^3 + 267300.*y.^4 - 314928.*y.^5 + 137781.*y.^6 + 40))/160;
            case 36
                gx = (243.*y.*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40).*(36.*x.*y - 17.*y - 34.*x + 27.*x.^2 + 9.*y.^2 + 8))/160;
                gy = (243.*x.*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160 + (243.*x.*y.*(18.*x + 18.*y - 17).*(822.*y - 6075.*y.^2 + 20655.*y.^3 - 32805.*y.^4 + 19683.*y.^5 - 40))/160 + (243.*x.*y.*(61965.*y.^2 - 12150.*y - 131220.*y.^3 + 98415.*y.^4 + 822).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160;
            case 37
                gx = -(243.*y.*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8).*(382.*x + 191.*y - 864.*x.*y + 486.*x.*y.^2 + 729.*x.^2.*y - 648.*x.^2 + 324.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/80;
                gy = - (243.*x.*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/80 - (243.*x.*y.*(1890.*y - 7290.*y.^2 + 8748.*y.^3 - 150).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/80 - (243.*x.*y.*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8).*(486.*x.*y - 432.*y - 432.*x + 243.*x.^2 + 243.*y.^2 + 191))/80;
            case 38
                gx = (729.*y.*(33.*y - 162.*y.^2 + 243.*y.^3 - 2).*(4374.*x.^2.*y.^2 - 550.*y - 1100.*x + 4020.*x.*y - 4860.*x.*y.^2 - 7290.*x.^2.*y + 1944.*x.*y.^3 + 3888.*x.^3.*y + 3015.*x.^2 - 3240.*x.^3 + 1215.*x.^4 + 1005.*y.^2 - 810.*y.^3 + 243.*y.^4 + 112))/64;
                gy = (729.*x.*(550.*x + 4796.*y - 368874.*x.^2.*y.^2 + 1371978.*x.^2.*y.^3 + 244944.*x.^3.*y.^2 - 2066715.*x.^2.*y.^4 - 708588.*x.^3.*y.^3 - 59049.*x.^4.*y.^2 + 1062882.*x.^2.*y.^5 + 590490.*x.^3.*y.^4 + 118098.*x.^4.*y.^3 - 22170.*x.*y + 240435.*x.*y.^2 + 38025.*x.^2.*y - 1082808.*x.*y.^3 - 28674.*x.^3.*y + 2285415.*x.*y.^4 + 8019.*x.^4.*y - 2243862.*x.*y.^5 + 826686.*x.*y.^6 - 1005.*x.^2 + 810.*x.^3 - 243.*x.^4 - 57456.*y.^2 + 302202.*y.^3 - 809190.*y.^4 + 1150362.*y.^5 - 826686.*y.^6 + 236196.*y.^7 - 112))/32;
            case 39
                gx = -(243.*y.*(81.*y.^2 - 27.*y + 2).*(7516.*x + 3758.*y - 153090.*x.^2.*y.^2 + 65610.*x.^2.*y.^3 + 87480.*x.^3.*y.^2 - 39900.*x.*y + 78570.*x.*y.^2 + 117855.*x.^2.*y - 68040.*x.*y.^3 - 136080.*x.^3.*y + 21870.*x.*y.^4 + 54675.*x.^4.*y - 29925.*x.^2 + 52380.*x.^3 - 42525.*x.^4 + 13122.*x.^5 - 9975.*y.^2 + 13095.*y.^3 - 8505.*y.^4 + 2187.*y.^5 - 560))/80;
                gy = -(243.*x.*(7516.*x + 45272.*y - 5912190.*x.^2.*y.^2 + 18414540.*x.^2.*y.^3 + 6068925.*x.^3.*y.^2 - 23619600.*x.^2.*y.^4 - 13384440.*x.^3.*y.^3 - 2952450.*x.^4.*y.^2 + 10628820.*x.^2.*y.^5 + 8857350.*x.^3.*y.^4 + 3542940.*x.^4.*y.^3 + 531441.*x.^5.*y.^2 - 282732.*x.*y + 2764854.*x.*y.^2 + 695790.*x.^2.*y - 10978740.*x.*y.^3 - 843210.*x.^3.*y + 20612475.*x.*y.^4 + 503010.*x.^4.*y - 18305190.*x.*y.^5 - 118098.*x.^5.*y + 6200145.*x.*y.^6 - 19950.*x.^2 + 26190.*x.^3 - 17010.*x.^4 + 4374.*x.^5 - 500328.*y.^2 + 2399652.*y.^3 - 5892750.*y.^4 + 7768224.*y.^5 - 5235678.*y.^6 + 1417176.*y.^7 - 1120))/80;
            case 40
                gx = (243.*y.*(9.*y - 1).*(2733750.*x.^2.*y.^2 - 20072.*y - 40144.*x - 2558790.*x.^2.*y.^3 - 3411720.*x.^3.*y.^2 + 885735.*x.^2.*y.^4 + 1574640.*x.^3.*y.^3 + 1476225.*x.^4.*y.^2 + 294888.*x.*y - 852930.*x.*y.^2 - 1279395.*x.^2.*y + 1215000.*x.*y.^3 + 2430000.*x.^3.*y - 852930.*x.*y.^4 - 2132325.*x.^4.*y + 236196.*x.*y.^5 + 708588.*x.^5.*y + 221166.*x.^2 - 568620.*x.^3 + 759375.*x.^4 - 511758.*x.^5 + 137781.*x.^6 + 73722.*y.^2 - 142155.*y.^3 + 151875.*y.^4 - 85293.*y.^5 + 19683.*y.^6 + 2240))/160;
                gy = (243.*x.*(20072.*x + 80464.*y - 14248305.*x.^2.*y.^2 + 36216720.*x.^2.*y.^3 + 18961290.*x.^3.*y.^2 - 39858075.*x.^2.*y.^4 - 32280120.*x.^3.*y.^3 - 12400290.*x.^4.*y.^2 + 15943230.*x.^2.*y.^5 + 17714700.*x.^3.*y.^4 + 10628820.*x.^4.*y.^3 + 3188646.*x.^5.*y.^2 - 656184.*x.*y + 5260383.*x.*y.^2 + 2179926.*x.^2.*y - 17782740.*x.*y.^3 - 3773790.*x.^3.*y + 29469825.*x.*y.^4 + 3586680.*x.^4.*y - 23737698.*x.*y.^5 - 1771470.*x.^5.*y + 7440174.*x.*y.^6 + 354294.*x.^6.*y - 73722.*x.^2 + 142155.*x.^3 - 151875.*x.^4 + 85293.*x.^5 - 19683.*x.^6 - 763110.*y.^2 + 3222612.*y.^3 - 7156350.*y.^4 + 8713008.*y.^5 - 5511240.*y.^6 + 1417176.*y.^7 - 2240))/160;
            case 41
                gx = (243.*y.*(80464.*x + 20072.*y - 14248305.*x.^2.*y.^2 + 18961290.*x.^2.*y.^3 + 36216720.*x.^3.*y.^2 - 12400290.*x.^2.*y.^4 - 32280120.*x.^3.*y.^3 - 39858075.*x.^4.*y.^2 + 3188646.*x.^2.*y.^5 + 10628820.*x.^3.*y.^4 + 17714700.*x.^4.*y.^3 + 15943230.*x.^5.*y.^2 - 656184.*x.*y + 2179926.*x.*y.^2 + 5260383.*x.^2.*y - 3773790.*x.*y.^3 - 17782740.*x.^3.*y + 3586680.*x.*y.^4 + 29469825.*x.^4.*y - 1771470.*x.*y.^5 - 23737698.*x.^5.*y + 354294.*x.*y.^6 + 7440174.*x.^6.*y - 763110.*x.^2 + 3222612.*x.^3 - 7156350.*x.^4 + 8713008.*x.^5 - 5511240.*x.^6 + 1417176.*x.^7 - 73722.*y.^2 + 142155.*y.^3 - 151875.*y.^4 + 85293.*y.^5 - 19683.*y.^6 - 2240))/160;
                gy = (243.*x.*(9.*x - 1).*(2733750.*x.^2.*y.^2 - 40144.*y - 20072.*x - 3411720.*x.^2.*y.^3 - 2558790.*x.^3.*y.^2 + 1476225.*x.^2.*y.^4 + 1574640.*x.^3.*y.^3 + 885735.*x.^4.*y.^2 + 294888.*x.*y - 1279395.*x.*y.^2 - 852930.*x.^2.*y + 2430000.*x.*y.^3 + 1215000.*x.^3.*y - 2132325.*x.*y.^4 - 852930.*x.^4.*y + 708588.*x.*y.^5 + 236196.*x.^5.*y + 73722.*x.^2 - 142155.*x.^3 + 151875.*x.^4 - 85293.*x.^5 + 19683.*x.^6 + 221166.*y.^2 - 568620.*y.^3 + 759375.*y.^4 - 511758.*y.^5 + 137781.*y.^6 + 2240))/160;
            case 42
                gx = -(243.*y.*(45272.*x + 7516.*y - 5912190.*x.^2.*y.^2 + 6068925.*x.^2.*y.^3 + 18414540.*x.^3.*y.^2 - 2952450.*x.^2.*y.^4 - 13384440.*x.^3.*y.^3 - 23619600.*x.^4.*y.^2 + 531441.*x.^2.*y.^5 + 3542940.*x.^3.*y.^4 + 8857350.*x.^4.*y.^3 + 10628820.*x.^5.*y.^2 - 282732.*x.*y + 695790.*x.*y.^2 + 2764854.*x.^2.*y - 843210.*x.*y.^3 - 10978740.*x.^3.*y + 503010.*x.*y.^4 + 20612475.*x.^4.*y - 118098.*x.*y.^5 - 18305190.*x.^5.*y + 6200145.*x.^6.*y - 500328.*x.^2 + 2399652.*x.^3 - 5892750.*x.^4 + 7768224.*x.^5 - 5235678.*x.^6 + 1417176.*x.^7 - 19950.*y.^2 + 26190.*y.^3 - 17010.*y.^4 + 4374.*y.^5 - 1120))/80;
                gy = -(243.*x.*(81.*x.^2 - 27.*x + 2).*(3758.*x + 7516.*y - 153090.*x.^2.*y.^2 + 87480.*x.^2.*y.^3 + 65610.*x.^3.*y.^2 - 39900.*x.*y + 117855.*x.*y.^2 + 78570.*x.^2.*y - 136080.*x.*y.^3 - 68040.*x.^3.*y + 54675.*x.*y.^4 + 21870.*x.^4.*y - 9975.*x.^2 + 13095.*x.^3 - 8505.*x.^4 + 2187.*x.^5 - 29925.*y.^2 + 52380.*y.^3 - 42525.*y.^4 + 13122.*y.^5 - 560))/80;
            case 43
                gx = (729.*y.*(4796.*x + 550.*y - 368874.*x.^2.*y.^2 + 244944.*x.^2.*y.^3 + 1371978.*x.^3.*y.^2 - 59049.*x.^2.*y.^4 - 708588.*x.^3.*y.^3 - 2066715.*x.^4.*y.^2 + 118098.*x.^3.*y.^4 + 590490.*x.^4.*y.^3 + 1062882.*x.^5.*y.^2 - 22170.*x.*y + 38025.*x.*y.^2 + 240435.*x.^2.*y - 28674.*x.*y.^3 - 1082808.*x.^3.*y + 8019.*x.*y.^4 + 2285415.*x.^4.*y - 2243862.*x.^5.*y + 826686.*x.^6.*y - 57456.*x.^2 + 302202.*x.^3 - 809190.*x.^4 + 1150362.*x.^5 - 826686.*x.^6 + 236196.*x.^7 - 1005.*y.^2 + 810.*y.^3 - 243.*y.^4 - 112))/32;
                gy = (729.*x.*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(4374.*x.^2.*y.^2 - 1100.*y - 550.*x + 4020.*x.*y - 7290.*x.*y.^2 - 4860.*x.^2.*y + 3888.*x.*y.^3 + 1944.*x.^3.*y + 1005.*x.^2 - 810.*x.^3 + 243.*x.^4 + 3015.*y.^2 - 3240.*y.^3 + 1215.*y.^4 + 112))/64;
            case 44
                gx = - (243.*y.*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/80 - (243.*x.*y.*(1890.*x - 7290.*x.^2 + 8748.*x.^3 - 150).*(191.*x + 191.*y - 432.*x.*y + 243.*x.*y.^2 + 243.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 216.*y.^2 + 81.*y.^3 - 56))/80 - (243.*x.*y.*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8).*(486.*x.*y - 432.*y - 432.*x + 243.*x.^2 + 243.*y.^2 + 191))/80;
                gy = -(243.*x.*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8).*(191.*x + 382.*y - 864.*x.*y + 729.*x.*y.^2 + 486.*x.^2.*y - 216.*x.^2 + 81.*x.^3 - 648.*y.^2 + 324.*y.^3 - 56))/80;
            case 45
                gx = (243.*y.*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160 + (243.*x.*y.*(18.*x + 18.*y - 17).*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40))/160 + (243.*x.*y.*(61965.*x.^2 - 12150.*x - 131220.*x.^3 + 98415.*x.^4 + 822).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160;
                gy = (243.*x.*(822.*x - 6075.*x.^2 + 20655.*x.^3 - 32805.*x.^4 + 19683.*x.^5 - 40).*(36.*x.*y - 34.*y - 17.*x + 9.*x.^2 + 27.*y.^2 + 8))/160;
            case 46
                gx = -(729.*y.*(9.*y - 1).*(1213785.*x.^2.*y.^2 - 3758.*y - 17596.*x - 984150.*x.^2.*y.^3 - 1924560.*x.^3.*y.^2 + 295245.*x.^2.*y.^4 + 787320.*x.^3.*y.^3 + 984150.*x.^4.*y.^2 + 107544.*x.*y - 258120.*x.*y.^2 - 656505.*x.^2.*y + 303750.*x.*y.^3 + 1550340.*x.^3.*y - 174960.*x.*y.^4 - 1585575.*x.^4.*y + 39366.*x.*y.^5 + 590490.*x.^5.*y + 131391.*x.^2 - 411480.*x.^3 + 631800.*x.^4 - 472392.*x.^5 + 137781.*x.^6 + 9975.*y.^2 - 13095.*y.^3 + 8505.*y.^4 - 2187.*y.^5 + 560))/160;
                gy = -(729.*x.*(9.*x - 1).*(1213785.*x.^2.*y.^2 - 17596.*y - 3758.*x - 1924560.*x.^2.*y.^3 - 984150.*x.^3.*y.^2 + 984150.*x.^2.*y.^4 + 787320.*x.^3.*y.^3 + 295245.*x.^4.*y.^2 + 107544.*x.*y - 656505.*x.*y.^2 - 258120.*x.^2.*y + 1550340.*x.*y.^3 + 303750.*x.^3.*y - 1585575.*x.*y.^4 - 174960.*x.^4.*y + 590490.*x.*y.^5 + 39366.*x.^5.*y + 9975.*x.^2 - 13095.*x.^3 + 8505.*x.^4 - 2187.*x.^5 + 131391.*y.^2 - 411480.*y.^3 + 631800.*y.^4 - 472392.*y.^5 + 137781.*y.^6 + 560))/160;
            case 47
                gx = (729.*y.*(9.*y - 1).*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160 + (729.*x.*y.*(9.*y - 1).*(18.*x + 18.*y - 17).*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8))/160 + (729.*x.*y.*(9.*y - 1).*(1890.*x - 7290.*x.^2 + 8748.*x.^3 - 150).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160;
                gy = (729.*x.*(945.*x.^2 - 150.*x - 2430.*x.^3 + 2187.*x.^4 + 8).*(17.*x + 178.*y - 342.*x.*y + 486.*x.*y.^2 + 162.*x.^2.*y - 9.*x.^2 - 486.*y.^2 + 324.*y.^3 - 8))/160;
            case 48
                gx = (729.*y.*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8).*(178.*x + 17.*y - 342.*x.*y + 162.*x.*y.^2 + 486.*x.^2.*y - 486.*x.^2 + 324.*x.^3 - 9.*y.^2 - 8))/160;
                gy = (729.*x.*(9.*x - 1).*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160 + (729.*x.*y.*(9.*x - 1).*(18.*x + 18.*y - 17).*(945.*y.^2 - 150.*y - 2430.*y.^3 + 2187.*y.^4 + 8))/160 + (729.*x.*y.*(9.*x - 1).*(1890.*y - 7290.*y.^2 + 8748.*y.^3 - 150).*(18.*x.*y - 17.*y - 17.*x + 9.*x.^2 + 9.*y.^2 + 8))/160;
            case 49
                gx = (243.*y.*(81.*y.^2 - 27.*y + 2).*(596.*x + 34.*y - 4374.*x.^2.*y.^2 + 8748.*x.^3.*y.^2 - 1194.*x.*y + 594.*x.*y.^2 + 10044.*x.^2.*y - 28188.*x.^3.*y + 21870.*x.^4.*y - 5625.*x.^2 + 19980.*x.^3 - 27945.*x.^4 + 13122.*x.^5 - 18.*y.^2 - 16))/32;
                gy = (243.*x.*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(2187.*x.^2.*y.^2 - 500.*y - 34.*x + 990.*x.*y - 5589.*x.*y.^2 - 486.*x.^2.*y + 5832.*x.*y.^3 + 18.*x.^2 + 3375.*y.^2 - 6480.*y.^3 + 3645.*y.^4 + 16))/32;
            case 50
                gx = (243.*y.*(33.*y - 162.*y.^2 + 243.*y.^3 - 2).*(2187.*x.^2.*y.^2 - 34.*y - 500.*x + 990.*x.*y - 486.*x.*y.^2 - 5589.*x.^2.*y + 5832.*x.^3.*y + 3375.*x.^2 - 6480.*x.^3 + 3645.*x.^4 + 18.*y.^2 + 16))/32;
                gy = (243.*x.*(81.*x.^2 - 27.*x + 2).*(34.*x + 596.*y - 4374.*x.^2.*y.^2 + 8748.*x.^2.*y.^3 - 1194.*x.*y + 10044.*x.*y.^2 + 594.*x.^2.*y - 28188.*x.*y.^3 + 21870.*x.*y.^4 - 18.*x.^2 - 5625.*y.^2 + 19980.*y.^3 - 27945.*y.^4 + 13122.*y.^5 - 16))/32;
            case 51
                gx = -(243.*y.*(33.*y - 162.*y.^2 + 243.*y.^3 - 2).*(6561.*x.^2.*y.^2 - 191.*y - 1390.*x + 4302.*x.*y - 4374.*x.*y.^2 - 12393.*x.^2.*y + 1458.*x.*y.^3 + 8748.*x.^3.*y + 5805.*x.^2 - 8100.*x.^3 + 3645.*x.^4 + 216.*y.^2 - 81.*y.^3 + 56))/32;
                gy = -(243.*x.*(9.*x - 1).*(129033.*x.^2.*y.^2 - 4460.*y - 382.*x - 367416.*x.^2.*y.^3 - 39366.*x.^3.*y.^2 + 295245.*x.^2.*y.^4 + 78732.*x.^3.*y.^3 + 14334.*x.*y - 137052.*x.*y.^2 - 15228.*x.^2.*y + 497664.*x.*y.^3 + 5346.*x.^3.*y - 721710.*x.*y.^4 + 354294.*x.*y.^5 + 432.*x.^2 - 162.*x.^3 + 47421.*y.^2 - 207360.*y.^3 + 420390.*y.^4 - 393660.*y.^5 + 137781.*y.^6 + 112))/32;
            case 52
                gx = (243.*y.*(81.*y.^2 - 27.*y + 2).*(3116.*x + 550.*y - 69984.*x.^2.*y.^2 + 26244.*x.^2.*y.^3 + 52488.*x.^3.*y.^2 - 13920.*x.*y + 22950.*x.*y.^2 + 61560.*x.^2.*y - 16524.*x.*y.^3 - 91368.*x.^3.*y + 4374.*x.*y.^4 + 43740.*x.^4.*y - 17865.*x.^2 + 39420.*x.^3 - 37665.*x.^4 + 13122.*x.^5 - 1005.*y.^2 + 810.*y.^3 - 243.*y.^4 - 112))/32;
                gy = (243.*x.*(9.*x - 1).*(449793.*x.^2.*y.^2 - 8248.*y - 1100.*x - 944784.*x.^2.*y.^3 - 275562.*x.^3.*y.^2 + 590490.*x.^2.*y.^4 + 314928.*x.^3.*y.^3 + 59049.*x.^4.*y.^2 + 37740.*x.*y - 311040.*x.*y.^2 - 63990.*x.^2.*y + 921456.*x.*y.^3 + 47628.*x.^3.*y - 1115370.*x.*y.^4 - 13122.*x.^4.*y + 472392.*x.*y.^5 + 2010.*x.^2 - 1620.*x.^3 + 486.*x.^4 + 77796.*y.^2 - 293220.*y.^3 + 518805.*y.^4 - 433026.*y.^5 + 137781.*y.^6 + 224))/32;
            case 53
                gx = (243.*y.*(9.*y - 1).*(449793.*x.^2.*y.^2 - 1100.*y - 8248.*x - 275562.*x.^2.*y.^3 - 944784.*x.^3.*y.^2 + 59049.*x.^2.*y.^4 + 314928.*x.^3.*y.^3 + 590490.*x.^4.*y.^2 + 37740.*x.*y - 63990.*x.*y.^2 - 311040.*x.^2.*y + 47628.*x.*y.^3 + 921456.*x.^3.*y - 13122.*x.*y.^4 - 1115370.*x.^4.*y + 472392.*x.^5.*y + 77796.*x.^2 - 293220.*x.^3 + 518805.*x.^4 - 433026.*x.^5 + 137781.*x.^6 + 2010.*y.^2 - 1620.*y.^3 + 486.*y.^4 + 224))/32;
                gy = (243.*x.*(81.*x.^2 - 27.*x + 2).*(550.*x + 3116.*y - 69984.*x.^2.*y.^2 + 52488.*x.^2.*y.^3 + 26244.*x.^3.*y.^2 - 13920.*x.*y + 61560.*x.*y.^2 + 22950.*x.^2.*y - 91368.*x.*y.^3 - 16524.*x.^3.*y + 43740.*x.*y.^4 + 4374.*x.^4.*y - 1005.*x.^2 + 810.*x.^3 - 243.*x.^4 - 17865.*y.^2 + 39420.*y.^3 - 37665.*y.^4 + 13122.*y.^5 - 112))/32;
            case 54
                gx = -(243.*y.*(9.*y - 1).*(129033.*x.^2.*y.^2 - 382.*y - 4460.*x - 39366.*x.^2.*y.^3 - 367416.*x.^3.*y.^2 + 78732.*x.^3.*y.^3 + 295245.*x.^4.*y.^2 + 14334.*x.*y - 15228.*x.*y.^2 - 137052.*x.^2.*y + 5346.*x.*y.^3 + 497664.*x.^3.*y - 721710.*x.^4.*y + 354294.*x.^5.*y + 47421.*x.^2 - 207360.*x.^3 + 420390.*x.^4 - 393660.*x.^5 + 137781.*x.^6 + 432.*y.^2 - 162.*y.^3 + 112))/32;
                gy = -(243.*x.*(33.*x - 162.*x.^2 + 243.*x.^3 - 2).*(6561.*x.^2.*y.^2 - 1390.*y - 191.*x + 4302.*x.*y - 12393.*x.*y.^2 - 4374.*x.^2.*y + 8748.*x.*y.^3 + 1458.*x.^3.*y + 216.*x.^2 - 81.*x.^3 + 5805.*y.^2 - 8100.*y.^3 + 3645.*y.^4 + 56))/32;
            case 55
                gx = -(27.*y.*(81.*y.^2 - 27.*y + 2).*(3788.*x + 382.*y - 72171.*x.^2.*y.^2 + 19683.*x.^2.*y.^3 + 78732.*x.^3.*y.^2 - 12042.*x.*y + 12636.*x.*y.^2 + 82863.*x.^2.*y - 4374.*x.*y.^3 - 166212.*x.^3.*y + 98415.*x.^4.*y - 30375.*x.^2 + 85860.*x.^3 - 98415.*x.^4 + 39366.*x.^5 - 432.*y.^2 + 162.*y.^3 - 112))/8;
                gy = -(27.*x.*(81.*x.^2 - 27.*x + 2).*(382.*x + 3788.*y - 72171.*x.^2.*y.^2 + 78732.*x.^2.*y.^3 + 19683.*x.^3.*y.^2 - 12042.*x.*y + 82863.*x.*y.^2 + 12636.*x.^2.*y - 166212.*x.*y.^3 - 4374.*x.^3.*y + 98415.*x.*y.^4 - 432.*x.^2 + 162.*x.^3 - 30375.*y.^2 + 85860.*y.^3 - 98415.*y.^4 + 39366.*y.^5 - 112))/8;
        end
    case 10
        switch i
            case 1
                gx = (177133.*x)/252 + (177133.*y)/252 - (26846875.*x.^2.*y.^2)/36 + (23478125.*x.^2.*y.^3)/9 + (23478125.*x.^3.*y.^2)/9 - (47265625.*x.^2.*y.^4)/9 - (189062500.*x.^3.*y.^3)/27 - (47265625.*x.^4.*y.^2)/9 + (55000000.*x.^2.*y.^5)/9 + (275000000.*x.^3.*y.^4)/27 + (275000000.*x.^4.*y.^3)/27 + (55000000.*x.^5.*y.^2)/9 - (34375000.*x.^2.*y.^6)/9 - (68750000.*x.^3.*y.^5)/9 - (85937500.*x.^4.*y.^4)/9 - (68750000.*x.^5.*y.^3)/9 - (34375000.*x.^6.*y.^2)/9 + (62500000.*x.^2.*y.^7)/63 + (62500000.*x.^3.*y.^6)/27 + (31250000.*x.^4.*y.^5)/9 + (31250000.*x.^5.*y.^4)/9 + (62500000.*x.^6.*y.^3)/27 + (62500000.*x.^7.*y.^2)/63 - (10511875.*x.*y)/756 + (42711625.*x.*y.^2)/378 + (42711625.*x.^2.*y)/378 - (26846875.*x.*y.^3)/54 - (26846875.*x.^3.*y)/54 + (23478125.*x.*y.^4)/18 + (23478125.*x.^4.*y)/18 - (18906250.*x.*y.^5)/9 - (18906250.*x.^5.*y)/9 + (55000000.*x.*y.^6)/27 + (55000000.*x.^6.*y)/27 - (68750000.*x.*y.^7)/63 - (68750000.*x.^7.*y)/63 + (15625000.*x.*y.^8)/63 + (15625000.*x.^8.*y)/63 - (10511875.*x.^2)/1512 + (42711625.*x.^3)/1134 - (26846875.*x.^4)/216 + (4695625.*x.^5)/18 - (9453125.*x.^6)/27 + (55000000.*x.^7)/189 - (8593750.*x.^8)/63 + (15625000.*x.^9)/567 - (10511875.*y.^2)/1512 + (42711625.*y.^3)/1134 - (26846875.*y.^4)/216 + (4695625.*y.^5)/18 - (9453125.*y.^6)/27 + (55000000.*y.^7)/189 - (8593750.*y.^8)/63 + (15625000.*y.^9)/567 - 7381/252;
                gy = (177133.*x)/252 + (177133.*y)/252 - (26846875.*x.^2.*y.^2)/36 + (23478125.*x.^2.*y.^3)/9 + (23478125.*x.^3.*y.^2)/9 - (47265625.*x.^2.*y.^4)/9 - (189062500.*x.^3.*y.^3)/27 - (47265625.*x.^4.*y.^2)/9 + (55000000.*x.^2.*y.^5)/9 + (275000000.*x.^3.*y.^4)/27 + (275000000.*x.^4.*y.^3)/27 + (55000000.*x.^5.*y.^2)/9 - (34375000.*x.^2.*y.^6)/9 - (68750000.*x.^3.*y.^5)/9 - (85937500.*x.^4.*y.^4)/9 - (68750000.*x.^5.*y.^3)/9 - (34375000.*x.^6.*y.^2)/9 + (62500000.*x.^2.*y.^7)/63 + (62500000.*x.^3.*y.^6)/27 + (31250000.*x.^4.*y.^5)/9 + (31250000.*x.^5.*y.^4)/9 + (62500000.*x.^6.*y.^3)/27 + (62500000.*x.^7.*y.^2)/63 - (10511875.*x.*y)/756 + (42711625.*x.*y.^2)/378 + (42711625.*x.^2.*y)/378 - (26846875.*x.*y.^3)/54 - (26846875.*x.^3.*y)/54 + (23478125.*x.*y.^4)/18 + (23478125.*x.^4.*y)/18 - (18906250.*x.*y.^5)/9 - (18906250.*x.^5.*y)/9 + (55000000.*x.*y.^6)/27 + (55000000.*x.^6.*y)/27 - (68750000.*x.*y.^7)/63 - (68750000.*x.^7.*y)/63 + (15625000.*x.*y.^8)/63 + (15625000.*x.^8.*y)/63 - (10511875.*x.^2)/1512 + (42711625.*x.^3)/1134 - (26846875.*x.^4)/216 + (4695625.*x.^5)/18 - (9453125.*x.^6)/27 + (55000000.*x.^7)/189 - (8593750.*x.^8)/63 + (15625000.*x.^9)/567 - (10511875.*y.^2)/1512 + (42711625.*y.^3)/1134 - (26846875.*y.^4)/216 + (4695625.*y.^5)/18 - (9453125.*y.^6)/27 + (55000000.*y.^7)/189 - (8593750.*y.^8)/63 + (15625000.*y.^9)/567 - 7381/252;
            case 2
                gx = (7129.*x)/126 - (162875.*x.^2)/168 + (4523000.*x.^3)/567 - (296875.*x.^4)/8 + (1883125.*x.^5)/18 - (546875.*x.^6)/3 + (36250000.*x.^7)/189 - (781250.*x.^8)/7 + (15625000.*x.^9)/567 - 1;
                gy = 0;
            case 3
                gx = 0;
                gy = (7129.*y)/126 - (162875.*y.^2)/168 + (4523000.*y.^3)/567 - (296875.*y.^4)/8 + (1883125.*y.^5)/18 - (546875.*y.^6)/3 + (36250000.*y.^7)/189 - (781250.*y.^8)/7 + (15625000.*y.^9)/567 - 1;
            case 4
                gx = (25.*y.*(147655.*x.^2 - 9132.*x - 1121400.*x.^3 + 4676875.*x.^4 - 11340000.*x.^5 + 15925000.*x.^6 - 12000000.*x.^7 + 3750000.*x.^8 + 168))/378;
                gy = (25.*x.*(147655.*x.^2 - 13698.*x - 841050.*x.^3 + 2806125.*x.^4 - 5670000.*x.^5 + 6825000.*x.^6 - 4500000.*x.^7 + 1250000.*x.^8 + 504))/1134;
            case 5
                gx = (25.*y.*(10.*y - 1).*(3267.*x - 49245.*x.^2 + 338450.*x.^3 - 1225000.*x.^4 + 2415000.*x.^5 - 2450000.*x.^6 + 1000000.*x.^7 - 63))/252;
                gy = (25.*x.*(20.*y - 1).*(3267.*x - 32830.*x.^2 + 169225.*x.^3 - 490000.*x.^4 + 805000.*x.^5 - 700000.*x.^6 + 250000.*x.^7 - 126))/504;
            case 6
                gx = (100.*y.*(50.*y.^2 - 15.*y + 1).*(6090.*x.^2 - 441.*x - 36750.*x.^3 + 109375.*x.^4 - 157500.*x.^5 + 87500.*x.^6 + 9))/189;
                gy = (50.*x.*(150.*y.^2 - 30.*y + 1).*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18))/189;
            case 7
                gx = (25.*y.*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(274.*x - 3375.*x.^2 + 17000.*x.^3 - 37500.*x.^4 + 30000.*x.^5 - 6))/108;
                gy = (25.*x.*(110.*y - 900.*y.^2 + 2000.*y.^3 - 3).*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6))/108;
            case 8
                gx = (y.*(2625.*x.^2 - 250.*x - 10000.*x.^3 + 12500.*x.^4 + 6).*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6))/9;
                gy = (x.*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(2625.*y.^2 - 250.*y - 10000.*y.^3 + 12500.*y.^4 + 6))/9;
            case 9
                gx = (25.*y.*(110.*x - 900.*x.^2 + 2000.*x.^3 - 3).*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6))/108;
                gy = (25.*x.*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(274.*y - 3375.*y.^2 + 17000.*y.^3 - 37500.*y.^4 + 30000.*y.^5 - 6))/108;
            case 10
                gx = (50.*y.*(150.*x.^2 - 30.*x + 1).*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18))/189;
                gy = (100.*x.*(50.*x.^2 - 15.*x + 1).*(6090.*y.^2 - 441.*y - 36750.*y.^3 + 109375.*y.^4 - 157500.*y.^5 + 87500.*y.^6 + 9))/189;
            case 11
                gx = (25.*y.*(20.*x - 1).*(3267.*y - 32830.*y.^2 + 169225.*y.^3 - 490000.*y.^4 + 805000.*y.^5 - 700000.*y.^6 + 250000.*y.^7 - 126))/504;
                gy = (25.*x.*(10.*x - 1).*(3267.*y - 49245.*y.^2 + 338450.*y.^3 - 1225000.*y.^4 + 2415000.*y.^5 - 2450000.*y.^6 + 1000000.*y.^7 - 63))/252;
            case 12
                gx = (25.*y.*(147655.*y.^2 - 13698.*y - 841050.*y.^3 + 2806125.*y.^4 - 5670000.*y.^5 + 6825000.*y.^6 - 4500000.*y.^7 + 1250000.*y.^8 + 504))/1134;
                gy = (25.*x.*(147655.*y.^2 - 9132.*y - 1121400.*y.^3 + 4676875.*y.^4 - 11340000.*y.^5 + 15925000.*y.^6 - 12000000.*y.^7 + 3750000.*y.^8 + 168))/378;
            case 13
                gx = -(25.*y.*(147655.*y.^2 - 13698.*y - 841050.*y.^3 + 2806125.*y.^4 - 5670000.*y.^5 + 6825000.*y.^6 - 4500000.*y.^7 + 1250000.*y.^8 + 504))/1134;
                gy = (38050.*x.*y)/63 - (13150.*y)/21 - (100.*x)/9 - (3691375.*x.*y.^2)/378 + (222500.*x.*y.^3)/3 - (16703125.*x.*y.^4)/54 + 750000.*x.*y.^5 - (28437500.*x.*y.^6)/27 + (50000000.*x.*y.^7)/63 - (15625000.*x.*y.^8)/63 + (4033825.*y.^2)/378 - (49435250.*y.^3)/567 + (21709375.*y.^4)/54 - (10090625.*y.^5)/9 + (52062500.*y.^6)/27 - (377500000.*y.^7)/189 + (71875000.*y.^8)/63 - (156250000.*y.^9)/567 + 100/9;
            case 14
                gx = (25.*y.*(20.*x + 20.*y - 19).*(3267.*y - 32830.*y.^2 + 169225.*y.^3 - 490000.*y.^4 + 805000.*y.^5 - 700000.*y.^6 + 250000.*y.^7 - 126))/504;
                gy = (25.*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9).*(3267.*y - 32830.*y.^2 + 169225.*y.^3 - 490000.*y.^4 + 805000.*y.^5 - 700000.*y.^6 + 250000.*y.^7 - 126))/504 + (25.*y.*(20.*x + 20.*y - 19).*(3267.*y - 32830.*y.^2 + 169225.*y.^3 - 490000.*y.^4 + 805000.*y.^5 - 700000.*y.^6 + 250000.*y.^7 - 126))/504 + (25.*y.*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9).*(507675.*y.^2 - 65660.*y - 1960000.*y.^3 + 4025000.*y.^4 - 4200000.*y.^5 + 1750000.*y.^6 + 3267))/504;
            case 15
                gx = -(50.*y.*(300.*x.*y - 270.*y - 270.*x + 150.*x.^2 + 150.*y.^2 + 121).*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18))/189;
                gy = - (50.*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/189 - (50.*y.*(8120.*y - 55125.*y.^2 + 175000.*y.^3 - 262500.*y.^4 + 150000.*y.^5 - 441).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/189 - (50.*y.*(300.*x.*y - 270.*y - 270.*x + 150.*x.^2 + 150.*y.^2 + 121).*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18))/189;
            case 16
                gx = (25.*y.*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(4310.*x + 4310.*y - 10200.*x.*y + 6000.*x.*y.^2 + 6000.*x.^2.*y - 5100.*x.^2 + 2000.*x.^3 - 5100.*y.^2 + 2000.*y.^3 - 1207))/108;
                gy = (25.*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(3000.*x.^2.*y.^2 - 1207.*y - 1207.*x + 4310.*x.*y - 5100.*x.*y.^2 - 5100.*x.^2.*y + 2000.*x.*y.^3 + 2000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/108 + (25.*y.*(12750.*y.^2 - 2250.*y - 30000.*y.^3 + 25000.*y.^4 + 137).*(3000.*x.^2.*y.^2 - 1207.*y - 1207.*x + 4310.*x.*y - 5100.*x.*y.^2 - 5100.*x.^2.*y + 2000.*x.*y.^3 + 2000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/108 + (25.*y.*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(4310.*x + 4310.*y - 10200.*x.*y + 6000.*x.*y.^2 + 6000.*x.^2.*y - 5100.*x.^2 + 2000.*x.^3 - 5100.*y.^2 + 2000.*y.^3 - 1207))/108;
            case 17
                gx = -(y.*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(75000.*x.^2.*y.^2 - 25000.*y - 25000.*x + 95250.*x.*y - 120000.*x.*y.^2 - 120000.*x.^2.*y + 50000.*x.*y.^3 + 50000.*x.^3.*y + 47625.*x.^2 - 40000.*x.^3 + 12500.*x.^4 + 47625.*y.^2 - 40000.*y.^3 + 12500.*y.^4 + 4881))/9;
                gy = (17250625.*x.^2.*y.^2)/3 - 27508.*y - 3254.*x - (322287500.*x.^2.*y.^3)/9 - 6346875.*x.^3.*y.^2 + (1029687500.*x.^2.*y.^4)/9 + (103750000.*x.^3.*y.^3)/3 + 3437500.*x.^4.*y.^2 - (581875000.*x.^2.*y.^5)/3 - (807812500.*x.^3.*y.^4)/9 - (143750000.*x.^4.*y.^3)/9 - (2187500.*x.^5.*y.^2)/3 + (1487500000.*x.^2.*y.^6)/9 + (325000000.*x.^3.*y.^5)/3 + 31250000.*x.^4.*y.^4 + (25000000.*x.^5.*y.^3)/9 - (500000000.*x.^2.*y.^7)/9 - (437500000.*x.^3.*y.^6)/9 - (62500000.*x.^4.*y.^5)/3 - (31250000.*x.^5.*y.^4)/9 + (506750.*x.*y)/3 - (7681625.*x.*y.^2)/3 - (3696500.*x.^2.*y)/9 + (161082500.*x.*y.^3)/9 + (4448750.*x.^3.*y)/9 - 67471875.*x.*y.^4 - (2650000.*x.^4.*y)/9 + (436250000.*x.*y.^5)/3 + (625000.*x.^5.*y)/9 - (1610000000.*x.*y.^6)/9 + (350000000.*x.*y.^7)/3 - 31250000.*x.*y.^8 + (25000.*x.^2)/3 - (31750.*x.^3)/3 + (20000.*x.^4)/3 - (5000.*x.^5)/3 + 448875.*y.^2 - (31274500.*y.^3)/9 + (135371875.*y.^4)/9 - (117216250.*y.^5)/3 + (560000000.*y.^6)/9 - (535000000.*y.^7)/9 + 31250000.*y.^8 - (62500000.*y.^9)/9 + 504;
            case 18
                gx = (25.*y.*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(44524.*x + 44524.*y - 675000.*x.^2.*y.^2 + 300000.*x.^2.*y.^3 + 300000.*x.^3.*y.^2 - 245250.*x.*y + 501000.*x.*y.^2 + 501000.*x.^2.*y - 450000.*x.*y.^3 - 450000.*x.^3.*y + 150000.*x.*y.^4 + 150000.*x.^4.*y - 122625.*x.^2 + 167000.*x.^3 - 112500.*x.^4 + 30000.*x.^5 - 122625.*y.^2 + 167000.*y.^3 - 112500.*y.^4 + 30000.*y.^5 - 6393))/108;
                gy = (53275.*x)/12 + (168775.*y)/6 - (118120625.*x.^2.*y.^2)/12 + (1559275000.*x.^2.*y.^3)/27 + (138265625.*x.^3.*y.^2)/9 - (517578125.*x.^2.*y.^4)/3 - (693437500.*x.^3.*y.^3)/9 - (39453125.*x.^4.*y.^2)/3 + 273437500.*x.^2.*y.^5 + (4890625000.*x.^3.*y.^4)/27 + (1468750000.*x.^4.*y.^3)/27 + (17500000.*x.^5.*y.^2)/3 - 218750000.*x.^2.*y.^6 - (593750000.*x.^3.*y.^5)/3 - (273437500.*x.^4.*y.^4)/3 - 18750000.*x.^5.*y.^3 - (3125000.*x.^6.*y.^2)/3 + (625000000.*x.^2.*y.^7)/9 + (2187500000.*x.^3.*y.^6)/27 + (156250000.*x.^4.*y.^5)/3 + (156250000.*x.^5.*y.^4)/9 + (62500000.*x.^6.*y.^3)/27 - (4043225.*x.*y)/18 + (118364875.*x.*y.^2)/36 + (39807125.*x.^2.*y)/54 - (198325625.*x.*y.^3)/9 - (22909375.*x.^3.*y)/18 + (2142875000.*x.*y.^4)/27 + (32921875.*x.^4.*y)/27 - (490375000.*x.*y.^5)/3 - (1843750.*x.^5.*y)/3 + (5201875000.*x.*y.^6)/27 + (3437500.*x.^6.*y)/27 - (362500000.*x.*y.^7)/3 + 31250000.*x.*y.^8 - (278275.*x.^2)/18 + (340625.*x.^3)/12 - (521875.*x.^4)/18 + 15625.*x.^5 - (31250.*x.^6)/9 - (1792225.*y.^2)/4 + (91073375.*y.^3)/27 - (510353125.*y.^4)/36 + (107321875.*y.^5)/3 - (498968750.*y.^6)/9 + (155000000.*y.^7)/3 - 26562500.*y.^8 + (156250000.*y.^9)/27 - 525;
            case 19
                gx = -(50.*y.*(50.*y.^2 - 15.*y + 1).*(7612500.*x.^2.*y.^2 - 152978.*y - 152978.*x - 7350000.*x.^2.*y.^3 - 7350000.*x.^3.*y.^2 + 2625000.*x.^2.*y.^4 + 3500000.*x.^3.*y.^3 + 2625000.*x.^4.*y.^2 + 1158360.*x.*y - 3454500.*x.*y.^2 - 3454500.*x.^2.*y + 5075000.*x.*y.^3 + 5075000.*x.^3.*y - 3675000.*x.*y.^4 - 3675000.*x.^4.*y + 1050000.*x.*y.^5 + 1050000.*x.^5.*y + 579180.*x.^2 - 1151500.*x.^3 + 1268750.*x.^4 - 735000.*x.^5 + 175000.*x.^6 + 579180.*y.^2 - 1151500.*y.^3 + 1268750.*y.^4 - 735000.*y.^5 + 175000.*y.^6 + 16566))/189;
                gy = (101710000.*x.^2.*y.^2)/9 - (1308200.*y)/63 - (276100.*x)/63 - (1640150000.*x.^2.*y.^3)/27 - (70150000.*x.^3.*y.^2)/3 + (1503125000.*x.^2.*y.^4)/9 + (2802500000.*x.^3.*y.^3)/27 + (251875000.*x.^4.*y.^2)/9 - (2213750000.*x.^2.*y.^5)/9 - (5875000000.*x.^3.*y.^4)/27 - (875000000.*x.^4.*y.^3)/9 - (173125000.*x.^5.*y.^2)/9 + (1662500000.*x.^2.*y.^6)/9 + (1937500000.*x.^3.*y.^5)/9 + (1250000000.*x.^4.*y.^4)/9 + (425000000.*x.^5.*y.^3)/9 + (62500000.*x.^6.*y.^2)/9 - (500000000.*x.^2.*y.^7)/9 - (2187500000.*x.^3.*y.^6)/27 - (625000000.*x.^4.*y.^5)/9 - (312500000.*x.^5.*y.^4)/9 - (250000000.*x.^6.*y.^3)/27 - (62500000.*x.^7.*y.^2)/63 + (40146800.*x.*y)/189 - (20567500.*x.*y.^2)/7 - 913500.*x.^2.*y + (499660000.*x.*y.^3)/27 + (57820000.*x.^3.*y)/27 - (1696437500.*x.*y.^4)/27 - (79812500.*x.^4.*y)/27 + (1104875000.*x.*y.^5)/9 + (21625000.*x.^5.*y)/9 - (3731875000.*x.*y.^6)/27 - (28750000.*x.^6.*y)/27 + (250000000.*x.*y.^7)/3 + (12500000.*x.^7.*y)/63 - (62500000.*x.*y.^8)/3 + (546350.*x.^2)/27 - (1379000.*x.^3)/27 + (2056250.*x.^4)/27 - (1812500.*x.^5)/27 + (875000.*x.^6)/27 - (1250000.*x.^7)/189 + (20028950.*y.^2)/63 - (433739000.*y.^3)/189 + (83431250.*y.^4)/9 - (67737500.*y.^5)/3 + (305375000.*y.^6)/9 - (1940000000.*y.^7)/63 + (325000000.*y.^8)/21 - (625000000.*y.^9)/189 + 400;
            case 20
                gx = (25.*y.*(10.*y - 1).*(790254.*x + 790254.*y - 109200000.*x.^2.*y.^2 + 174300000.*x.^2.*y.^3 + 174300000.*x.^3.*y.^2 - 136500000.*x.^2.*y.^4 - 182000000.*x.^3.*y.^3 - 136500000.*x.^4.*y.^2 + 42000000.*x.^2.*y.^5 + 70000000.*x.^3.*y.^4 + 70000000.*x.^4.*y.^3 + 42000000.*x.^5.*y.^2 - 8064420.*x.*y + 33530700.*x.*y.^2 + 33530700.*x.^2.*y - 72800000.*x.*y.^3 - 72800000.*x.^3.*y + 87150000.*x.*y.^4 + 87150000.*x.^4.*y - 54600000.*x.*y.^5 - 54600000.*x.^5.*y + 14000000.*x.*y.^6 + 14000000.*x.^6.*y - 4032210.*x.^2 + 11176900.*x.^3 - 18200000.*x.^4 + 17430000.*x.^5 - 9100000.*x.^6 + 2000000.*x.^7 - 4032210.*y.^2 + 11176900.*y.^3 - 18200000.*y.^4 + 17430000.*y.^5 - 9100000.*y.^6 + 2000000.*y.^7 - 64818))/504;
                gy = (90025.*x)/28 + (153025.*y)/14 - (33980625.*x.^2.*y.^2)/4 + (364381250.*x.^2.*y.^3)/9 + (66146875.*x.^3.*y.^2)/3 - (909765625.*x.^2.*y.^4)/9 - 83750000.*x.^3.*y.^3 - (100703125.*x.^4.*y.^2)/3 + 137812500.*x.^2.*y.^5 + (1398437500.*x.^3.*y.^4)/9 + (859375000.*x.^4.*y.^3)/9 + 30000000.*x.^5.*y.^2 - (875000000.*x.^2.*y.^6)/9 - (418750000.*x.^3.*y.^5)/3 - 117187500.*x.^4.*y.^4 - (512500000.*x.^5.*y.^3)/9 - (43750000.*x.^6.*y.^2)/3 + (250000000.*x.^2.*y.^7)/9 + (437500000.*x.^3.*y.^6)/9 + (156250000.*x.^4.*y.^5)/3 + (312500000.*x.^5.*y.^4)/9 + (125000000.*x.^6.*y.^3)/9 + (62500000.*x.^7.*y.^2)/21 - (1997825.*x.*y)/14 + (49728125.*x.*y.^2)/28 + (16632250.*x.^2.*y)/21 - (91962500.*x.*y.^3)/9 - (21980000.*x.^3.*y)/9 + 32234375.*x.*y.^4 + (27465625.*x.^4.*y)/6 - (178062500.*x.*y.^5)/3 - (48062500.*x.^5.*y)/9 + (573125000.*x.*y.^6)/9 + (34062500.*x.^6.*y)/9 - (775000000.*x.*y.^7)/21 - (31250000.*x.^7.*y)/21 + (62500000.*x.*y.^8)/7 + (15625000.*x.^8.*y)/63 - (1097575.*x.^2)/56 + (2400125.*x.^3)/36 - (9979375.*x.^4)/72 + (1625000.*x.^5)/9 - (1296875.*x.^6)/9 + (4062500.*x.^7)/63 - (781250.*x.^8)/63 - (8694225.*y.^2)/56 + (66191750.*y.^3)/63 - (289909375.*y.^4)/72 + (56396875.*y.^5)/6 - (122828125.*y.^6)/9 + (758750000.*y.^7)/63 - (41406250.*y.^8)/7 + (78125000.*y.^9)/63 - 225;
            case 21
                gx = -(25.*y.*(528333750.*x.^2.*y.^2 - 1438434.*y - 1438434.*x - 1266300000.*x.^2.*y.^3 - 1266300000.*x.^3.*y.^2 + 1661625000.*x.^2.*y.^4 + 2215500000.*x.^3.*y.^3 + 1661625000.*x.^4.*y.^2 - 1134000000.*x.^2.*y.^5 - 1890000000.*x.^3.*y.^4 - 1890000000.*x.^4.*y.^3 - 1134000000.*x.^5.*y.^2 + 315000000.*x.^2.*y.^6 + 630000000.*x.^3.*y.^5 + 787500000.*x.^4.*y.^4 + 630000000.*x.^5.*y.^3 + 315000000.*x.^6.*y.^2 + 19918230.*x.*y - 114174900.*x.*y.^2 - 114174900.*x.^2.*y + 352222500.*x.*y.^3 + 352222500.*x.^3.*y - 633150000.*x.*y.^4 - 633150000.*x.^4.*y + 664650000.*x.*y.^5 + 664650000.*x.^5.*y - 378000000.*x.*y.^6 - 378000000.*x.^6.*y + 90000000.*x.*y.^7 + 90000000.*x.^7.*y + 9959115.*x.^2 - 38058300.*x.^3 + 88055625.*x.^4 - 126630000.*x.^5 + 110775000.*x.^6 - 54000000.*x.^7 + 11250000.*x.^8 + 9959115.*y.^2 - 38058300.*y.^3 + 88055625.*y.^4 - 126630000.*y.^5 + 110775000.*y.^6 - 54000000.*y.^7 + 11250000.*y.^8 + 87498))/1134;
                gy = 3775625.*x.^2.*y.^2 - (243050.*y)/63 - (121525.*x)/63 - (419312500.*x.^2.*y.^3)/27 - (104828125.*x.^3.*y.^2)/9 + (104687500.*x.^2.*y.^4)/3 + (335000000.*x.^3.*y.^3)/9 + 20937500.*x.^4.*y.^2 - (131875000.*x.^2.*y.^5)/3 - (1648437500.*x.^3.*y.^4)/27 - (1318750000.*x.^4.*y.^3)/27 - (65937500.*x.^5.*y.^2)/3 + (87500000.*x.^2.*y.^6)/3 + 50000000.*x.^3.*y.^5 + (156250000.*x.^4.*y.^4)/3 + (100000000.*x.^5.*y.^3)/3 + 12500000.*x.^6.*y.^2 - (500000000.*x.^2.*y.^7)/63 - (437500000.*x.^3.*y.^6)/27 - (62500000.*x.^4.*y.^5)/3 - (156250000.*x.^5.*y.^4)/9 - (250000000.*x.^6.*y.^3)/27 - (62500000.*x.^7.*y.^2)/21 + (3995650.*x.*y)/63 - (82992625.*x.*y.^2)/126 - (82992625.*x.^2.*y)/189 + (30205000.*x.*y.^3)/9 + (15102500.*x.^3.*y)/9 - (524140625.*x.*y.^4)/54 - (104828125.*x.^4.*y)/27 + 16750000.*x.*y.^5 + (16750000.*x.^5.*y)/3 - (461562500.*x.*y.^6)/27 - (131875000.*x.^6.*y)/27 + (200000000.*x.*y.^7)/21 + (50000000.*x.^7.*y)/21 - (15625000.*x.*y.^8)/7 - (31250000.*x.^8.*y)/63 + (1997825.*x.^2)/126 - (82992625.*x.^3)/1134 + (3775625.*x.^4)/18 - (20965625.*x.^5)/54 + (4187500.*x.^6)/9 - (65937500.*x.^7)/189 + (3125000.*x.^8)/21 - (15625000.*x.^9)/567 + (1997825.*y.^2)/42 - (165985250.*y.^3)/567 + (18878125.*y.^4)/18 - (20965625.*y.^5)/9 + (29312500.*y.^6)/9 - (527500000.*y.^7)/189 + (9375000.*y.^8)/7 - (156250000.*y.^9)/567 + 100;
            case 22
                gx = 3775625.*x.^2.*y.^2 - (121525.*y)/63 - (243050.*x)/63 - (104828125.*x.^2.*y.^3)/9 - (419312500.*x.^3.*y.^2)/27 + 20937500.*x.^2.*y.^4 + (335000000.*x.^3.*y.^3)/9 + (104687500.*x.^4.*y.^2)/3 - (65937500.*x.^2.*y.^5)/3 - (1318750000.*x.^3.*y.^4)/27 - (1648437500.*x.^4.*y.^3)/27 - (131875000.*x.^5.*y.^2)/3 + 12500000.*x.^2.*y.^6 + (100000000.*x.^3.*y.^5)/3 + (156250000.*x.^4.*y.^4)/3 + 50000000.*x.^5.*y.^3 + (87500000.*x.^6.*y.^2)/3 - (62500000.*x.^2.*y.^7)/21 - (250000000.*x.^3.*y.^6)/27 - (156250000.*x.^4.*y.^5)/9 - (62500000.*x.^5.*y.^4)/3 - (437500000.*x.^6.*y.^3)/27 - (500000000.*x.^7.*y.^2)/63 + (3995650.*x.*y)/63 - (82992625.*x.*y.^2)/189 - (82992625.*x.^2.*y)/126 + (15102500.*x.*y.^3)/9 + (30205000.*x.^3.*y)/9 - (104828125.*x.*y.^4)/27 - (524140625.*x.^4.*y)/54 + (16750000.*x.*y.^5)/3 + 16750000.*x.^5.*y - (131875000.*x.*y.^6)/27 - (461562500.*x.^6.*y)/27 + (50000000.*x.*y.^7)/21 + (200000000.*x.^7.*y)/21 - (31250000.*x.*y.^8)/63 - (15625000.*x.^8.*y)/7 + (1997825.*x.^2)/42 - (165985250.*x.^3)/567 + (18878125.*x.^4)/18 - (20965625.*x.^5)/9 + (29312500.*x.^6)/9 - (527500000.*x.^7)/189 + (9375000.*x.^8)/7 - (156250000.*x.^9)/567 + (1997825.*y.^2)/126 - (82992625.*y.^3)/1134 + (3775625.*y.^4)/18 - (20965625.*y.^5)/54 + (4187500.*y.^6)/9 - (65937500.*y.^7)/189 + (3125000.*y.^8)/21 - (15625000.*y.^9)/567 + 100;
                gy = -(25.*x.*(528333750.*x.^2.*y.^2 - 1438434.*y - 1438434.*x - 1266300000.*x.^2.*y.^3 - 1266300000.*x.^3.*y.^2 + 1661625000.*x.^2.*y.^4 + 2215500000.*x.^3.*y.^3 + 1661625000.*x.^4.*y.^2 - 1134000000.*x.^2.*y.^5 - 1890000000.*x.^3.*y.^4 - 1890000000.*x.^4.*y.^3 - 1134000000.*x.^5.*y.^2 + 315000000.*x.^2.*y.^6 + 630000000.*x.^3.*y.^5 + 787500000.*x.^4.*y.^4 + 630000000.*x.^5.*y.^3 + 315000000.*x.^6.*y.^2 + 19918230.*x.*y - 114174900.*x.*y.^2 - 114174900.*x.^2.*y + 352222500.*x.*y.^3 + 352222500.*x.^3.*y - 633150000.*x.*y.^4 - 633150000.*x.^4.*y + 664650000.*x.*y.^5 + 664650000.*x.^5.*y - 378000000.*x.*y.^6 - 378000000.*x.^6.*y + 90000000.*x.*y.^7 + 90000000.*x.^7.*y + 9959115.*x.^2 - 38058300.*x.^3 + 88055625.*x.^4 - 126630000.*x.^5 + 110775000.*x.^6 - 54000000.*x.^7 + 11250000.*x.^8 + 9959115.*y.^2 - 38058300.*y.^3 + 88055625.*y.^4 - 126630000.*y.^5 + 110775000.*y.^6 - 54000000.*y.^7 + 11250000.*y.^8 + 87498))/1134;
            case 23
                gx = (153025.*x)/14 + (90025.*y)/28 - (33980625.*x.^2.*y.^2)/4 + (66146875.*x.^2.*y.^3)/3 + (364381250.*x.^3.*y.^2)/9 - (100703125.*x.^2.*y.^4)/3 - 83750000.*x.^3.*y.^3 - (909765625.*x.^4.*y.^2)/9 + 30000000.*x.^2.*y.^5 + (859375000.*x.^3.*y.^4)/9 + (1398437500.*x.^4.*y.^3)/9 + 137812500.*x.^5.*y.^2 - (43750000.*x.^2.*y.^6)/3 - (512500000.*x.^3.*y.^5)/9 - 117187500.*x.^4.*y.^4 - (418750000.*x.^5.*y.^3)/3 - (875000000.*x.^6.*y.^2)/9 + (62500000.*x.^2.*y.^7)/21 + (125000000.*x.^3.*y.^6)/9 + (312500000.*x.^4.*y.^5)/9 + (156250000.*x.^5.*y.^4)/3 + (437500000.*x.^6.*y.^3)/9 + (250000000.*x.^7.*y.^2)/9 - (1997825.*x.*y)/14 + (16632250.*x.*y.^2)/21 + (49728125.*x.^2.*y)/28 - (21980000.*x.*y.^3)/9 - (91962500.*x.^3.*y)/9 + (27465625.*x.*y.^4)/6 + 32234375.*x.^4.*y - (48062500.*x.*y.^5)/9 - (178062500.*x.^5.*y)/3 + (34062500.*x.*y.^6)/9 + (573125000.*x.^6.*y)/9 - (31250000.*x.*y.^7)/21 - (775000000.*x.^7.*y)/21 + (15625000.*x.*y.^8)/63 + (62500000.*x.^8.*y)/7 - (8694225.*x.^2)/56 + (66191750.*x.^3)/63 - (289909375.*x.^4)/72 + (56396875.*x.^5)/6 - (122828125.*x.^6)/9 + (758750000.*x.^7)/63 - (41406250.*x.^8)/7 + (78125000.*x.^9)/63 - (1097575.*y.^2)/56 + (2400125.*y.^3)/36 - (9979375.*y.^4)/72 + (1625000.*y.^5)/9 - (1296875.*y.^6)/9 + (4062500.*y.^7)/63 - (781250.*y.^8)/63 - 225;
                gy = (25.*x.*(10.*x - 1).*(790254.*x + 790254.*y - 109200000.*x.^2.*y.^2 + 174300000.*x.^2.*y.^3 + 174300000.*x.^3.*y.^2 - 136500000.*x.^2.*y.^4 - 182000000.*x.^3.*y.^3 - 136500000.*x.^4.*y.^2 + 42000000.*x.^2.*y.^5 + 70000000.*x.^3.*y.^4 + 70000000.*x.^4.*y.^3 + 42000000.*x.^5.*y.^2 - 8064420.*x.*y + 33530700.*x.*y.^2 + 33530700.*x.^2.*y - 72800000.*x.*y.^3 - 72800000.*x.^3.*y + 87150000.*x.*y.^4 + 87150000.*x.^4.*y - 54600000.*x.*y.^5 - 54600000.*x.^5.*y + 14000000.*x.*y.^6 + 14000000.*x.^6.*y - 4032210.*x.^2 + 11176900.*x.^3 - 18200000.*x.^4 + 17430000.*x.^5 - 9100000.*x.^6 + 2000000.*x.^7 - 4032210.*y.^2 + 11176900.*y.^3 - 18200000.*y.^4 + 17430000.*y.^5 - 9100000.*y.^6 + 2000000.*y.^7 - 64818))/504;
            case 24
                gx = (101710000.*x.^2.*y.^2)/9 - (276100.*y)/63 - (1308200.*x)/63 - (70150000.*x.^2.*y.^3)/3 - (1640150000.*x.^3.*y.^2)/27 + (251875000.*x.^2.*y.^4)/9 + (2802500000.*x.^3.*y.^3)/27 + (1503125000.*x.^4.*y.^2)/9 - (173125000.*x.^2.*y.^5)/9 - (875000000.*x.^3.*y.^4)/9 - (5875000000.*x.^4.*y.^3)/27 - (2213750000.*x.^5.*y.^2)/9 + (62500000.*x.^2.*y.^6)/9 + (425000000.*x.^3.*y.^5)/9 + (1250000000.*x.^4.*y.^4)/9 + (1937500000.*x.^5.*y.^3)/9 + (1662500000.*x.^6.*y.^2)/9 - (62500000.*x.^2.*y.^7)/63 - (250000000.*x.^3.*y.^6)/27 - (312500000.*x.^4.*y.^5)/9 - (625000000.*x.^5.*y.^4)/9 - (2187500000.*x.^6.*y.^3)/27 - (500000000.*x.^7.*y.^2)/9 + (40146800.*x.*y)/189 - 913500.*x.*y.^2 - (20567500.*x.^2.*y)/7 + (57820000.*x.*y.^3)/27 + (499660000.*x.^3.*y)/27 - (79812500.*x.*y.^4)/27 - (1696437500.*x.^4.*y)/27 + (21625000.*x.*y.^5)/9 + (1104875000.*x.^5.*y)/9 - (28750000.*x.*y.^6)/27 - (3731875000.*x.^6.*y)/27 + (12500000.*x.*y.^7)/63 + (250000000.*x.^7.*y)/3 - (62500000.*x.^8.*y)/3 + (20028950.*x.^2)/63 - (433739000.*x.^3)/189 + (83431250.*x.^4)/9 - (67737500.*x.^5)/3 + (305375000.*x.^6)/9 - (1940000000.*x.^7)/63 + (325000000.*x.^8)/21 - (625000000.*x.^9)/189 + (546350.*y.^2)/27 - (1379000.*y.^3)/27 + (2056250.*y.^4)/27 - (1812500.*y.^5)/27 + (875000.*y.^6)/27 - (1250000.*y.^7)/189 + 400;
                gy = -(50.*x.*(50.*x.^2 - 15.*x + 1).*(7612500.*x.^2.*y.^2 - 152978.*y - 152978.*x - 7350000.*x.^2.*y.^3 - 7350000.*x.^3.*y.^2 + 2625000.*x.^2.*y.^4 + 3500000.*x.^3.*y.^3 + 2625000.*x.^4.*y.^2 + 1158360.*x.*y - 3454500.*x.*y.^2 - 3454500.*x.^2.*y + 5075000.*x.*y.^3 + 5075000.*x.^3.*y - 3675000.*x.*y.^4 - 3675000.*x.^4.*y + 1050000.*x.*y.^5 + 1050000.*x.^5.*y + 579180.*x.^2 - 1151500.*x.^3 + 1268750.*x.^4 - 735000.*x.^5 + 175000.*x.^6 + 579180.*y.^2 - 1151500.*y.^3 + 1268750.*y.^4 - 735000.*y.^5 + 175000.*y.^6 + 16566))/189;
            case 25
                gx = (168775.*x)/6 + (53275.*y)/12 - (118120625.*x.^2.*y.^2)/12 + (138265625.*x.^2.*y.^3)/9 + (1559275000.*x.^3.*y.^2)/27 - (39453125.*x.^2.*y.^4)/3 - (693437500.*x.^3.*y.^3)/9 - (517578125.*x.^4.*y.^2)/3 + (17500000.*x.^2.*y.^5)/3 + (1468750000.*x.^3.*y.^4)/27 + (4890625000.*x.^4.*y.^3)/27 + 273437500.*x.^5.*y.^2 - (3125000.*x.^2.*y.^6)/3 - 18750000.*x.^3.*y.^5 - (273437500.*x.^4.*y.^4)/3 - (593750000.*x.^5.*y.^3)/3 - 218750000.*x.^6.*y.^2 + (62500000.*x.^3.*y.^6)/27 + (156250000.*x.^4.*y.^5)/9 + (156250000.*x.^5.*y.^4)/3 + (2187500000.*x.^6.*y.^3)/27 + (625000000.*x.^7.*y.^2)/9 - (4043225.*x.*y)/18 + (39807125.*x.*y.^2)/54 + (118364875.*x.^2.*y)/36 - (22909375.*x.*y.^3)/18 - (198325625.*x.^3.*y)/9 + (32921875.*x.*y.^4)/27 + (2142875000.*x.^4.*y)/27 - (1843750.*x.*y.^5)/3 - (490375000.*x.^5.*y)/3 + (3437500.*x.*y.^6)/27 + (5201875000.*x.^6.*y)/27 - (362500000.*x.^7.*y)/3 + 31250000.*x.^8.*y - (1792225.*x.^2)/4 + (91073375.*x.^3)/27 - (510353125.*x.^4)/36 + (107321875.*x.^5)/3 - (498968750.*x.^6)/9 + (155000000.*x.^7)/3 - 26562500.*x.^8 + (156250000.*x.^9)/27 - (278275.*y.^2)/18 + (340625.*y.^3)/12 - (521875.*y.^4)/18 + 15625.*y.^5 - (31250.*y.^6)/9 - 525;
                gy = (25.*x.*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(44524.*x + 44524.*y - 675000.*x.^2.*y.^2 + 300000.*x.^2.*y.^3 + 300000.*x.^3.*y.^2 - 245250.*x.*y + 501000.*x.*y.^2 + 501000.*x.^2.*y - 450000.*x.*y.^3 - 450000.*x.^3.*y + 150000.*x.*y.^4 + 150000.*x.^4.*y - 122625.*x.^2 + 167000.*x.^3 - 112500.*x.^4 + 30000.*x.^5 - 122625.*y.^2 + 167000.*y.^3 - 112500.*y.^4 + 30000.*y.^5 - 6393))/108;
            case 26
                gx = (17250625.*x.^2.*y.^2)/3 - 3254.*y - 27508.*x - 6346875.*x.^2.*y.^3 - (322287500.*x.^3.*y.^2)/9 + 3437500.*x.^2.*y.^4 + (103750000.*x.^3.*y.^3)/3 + (1029687500.*x.^4.*y.^2)/9 - (2187500.*x.^2.*y.^5)/3 - (143750000.*x.^3.*y.^4)/9 - (807812500.*x.^4.*y.^3)/9 - (581875000.*x.^5.*y.^2)/3 + (25000000.*x.^3.*y.^5)/9 + 31250000.*x.^4.*y.^4 + (325000000.*x.^5.*y.^3)/3 + (1487500000.*x.^6.*y.^2)/9 - (31250000.*x.^4.*y.^5)/9 - (62500000.*x.^5.*y.^4)/3 - (437500000.*x.^6.*y.^3)/9 - (500000000.*x.^7.*y.^2)/9 + (506750.*x.*y)/3 - (3696500.*x.*y.^2)/9 - (7681625.*x.^2.*y)/3 + (4448750.*x.*y.^3)/9 + (161082500.*x.^3.*y)/9 - (2650000.*x.*y.^4)/9 - 67471875.*x.^4.*y + (625000.*x.*y.^5)/9 + (436250000.*x.^5.*y)/3 - (1610000000.*x.^6.*y)/9 + (350000000.*x.^7.*y)/3 - 31250000.*x.^8.*y + 448875.*x.^2 - (31274500.*x.^3)/9 + (135371875.*x.^4)/9 - (117216250.*x.^5)/3 + (560000000.*x.^6)/9 - (535000000.*x.^7)/9 + 31250000.*x.^8 - (62500000.*x.^9)/9 + (25000.*y.^2)/3 - (31750.*y.^3)/3 + (20000.*y.^4)/3 - (5000.*y.^5)/3 + 504;
                gy = -(x.*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(75000.*x.^2.*y.^2 - 25000.*y - 25000.*x + 95250.*x.*y - 120000.*x.*y.^2 - 120000.*x.^2.*y + 50000.*x.*y.^3 + 50000.*x.^3.*y + 47625.*x.^2 - 40000.*x.^3 + 12500.*x.^4 + 47625.*y.^2 - 40000.*y.^3 + 12500.*y.^4 + 4881))/9;
            case 27
                gx = (25.*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(3000.*x.^2.*y.^2 - 1207.*y - 1207.*x + 4310.*x.*y - 5100.*x.*y.^2 - 5100.*x.^2.*y + 2000.*x.*y.^3 + 2000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/108 + (25.*x.*(12750.*x.^2 - 2250.*x - 30000.*x.^3 + 25000.*x.^4 + 137).*(3000.*x.^2.*y.^2 - 1207.*y - 1207.*x + 4310.*x.*y - 5100.*x.*y.^2 - 5100.*x.^2.*y + 2000.*x.*y.^3 + 2000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/108 + (25.*x.*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(4310.*x + 4310.*y - 10200.*x.*y + 6000.*x.*y.^2 + 6000.*x.^2.*y - 5100.*x.^2 + 2000.*x.^3 - 5100.*y.^2 + 2000.*y.^3 - 1207))/108;
                gy = (25.*x.*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(4310.*x + 4310.*y - 10200.*x.*y + 6000.*x.*y.^2 + 6000.*x.^2.*y - 5100.*x.^2 + 2000.*x.^3 - 5100.*y.^2 + 2000.*y.^3 - 1207))/108;
            case 28
                gx = - (50.*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/189 - (50.*x.*(8120.*x - 55125.*x.^2 + 175000.*x.^3 - 262500.*x.^4 + 150000.*x.^5 - 441).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/189 - (50.*x.*(300.*x.*y - 270.*y - 270.*x + 150.*x.^2 + 150.*y.^2 + 121).*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18))/189;
                gy = -(50.*x.*(300.*x.*y - 270.*y - 270.*x + 150.*x.^2 + 150.*y.^2 + 121).*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18))/189;
            case 29
                gx = (25.*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9).*(3267.*x - 32830.*x.^2 + 169225.*x.^3 - 490000.*x.^4 + 805000.*x.^5 - 700000.*x.^6 + 250000.*x.^7 - 126))/504 + (25.*x.*(20.*x + 20.*y - 19).*(3267.*x - 32830.*x.^2 + 169225.*x.^3 - 490000.*x.^4 + 805000.*x.^5 - 700000.*x.^6 + 250000.*x.^7 - 126))/504 + (25.*x.*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9).*(507675.*x.^2 - 65660.*x - 1960000.*x.^3 + 4025000.*x.^4 - 4200000.*x.^5 + 1750000.*x.^6 + 3267))/504;
                gy = (25.*x.*(20.*x + 20.*y - 19).*(3267.*x - 32830.*x.^2 + 169225.*x.^3 - 490000.*x.^4 + 805000.*x.^5 - 700000.*x.^6 + 250000.*x.^7 - 126))/504;
            case 30
                gx = (38050.*x.*y)/63 - (100.*y)/9 - (13150.*x)/21 - (3691375.*x.^2.*y)/378 + (222500.*x.^3.*y)/3 - (16703125.*x.^4.*y)/54 + 750000.*x.^5.*y - (28437500.*x.^6.*y)/27 + (50000000.*x.^7.*y)/63 - (15625000.*x.^8.*y)/63 + (4033825.*x.^2)/378 - (49435250.*x.^3)/567 + (21709375.*x.^4)/54 - (10090625.*x.^5)/9 + (52062500.*x.^6)/27 - (377500000.*x.^7)/189 + (71875000.*x.^8)/63 - (156250000.*x.^9)/567 + 100/9;
                gy = -(25.*x.*(147655.*x.^2 - 13698.*x - 841050.*x.^3 + 2806125.*x.^4 - 5670000.*x.^5 + 6825000.*x.^6 - 4500000.*x.^7 + 1250000.*x.^8 + 504))/1134;
            case 31
                gx = (125.*y.*(50296050.*x.^2.*y.^2 - 64818.*y - 129636.*x - 109200000.*x.^2.*y.^3 - 145600000.*x.^3.*y.^2 + 130725000.*x.^2.*y.^4 + 232400000.*x.^3.*y.^3 + 217875000.*x.^4.*y.^2 - 81900000.*x.^2.*y.^5 - 182000000.*x.^3.*y.^4 - 227500000.*x.^4.*y.^3 - 163800000.*x.^5.*y.^2 + 21000000.*x.^2.*y.^6 + 56000000.*x.^3.*y.^5 + 87500000.*x.^4.*y.^4 + 84000000.*x.^5.*y.^3 + 49000000.*x.^6.*y.^2 + 1580508.*x.*y - 8064420.*x.*y.^2 - 12096630.*x.^2.*y + 22353800.*x.*y.^3 + 44707600.*x.^3.*y - 36400000.*x.*y.^4 - 91000000.*x.^4.*y + 34860000.*x.*y.^5 + 104580000.*x.^5.*y - 18200000.*x.*y.^6 - 63700000.*x.^6.*y + 4000000.*x.*y.^7 + 16000000.*x.^7.*y + 1185381.*x.^2 - 5376280.*x.^3 + 13971125.*x.^4 - 21840000.*x.^5 + 20335000.*x.^6 - 10400000.*x.^7 + 2250000.*x.^8 + 395127.*y.^2 - 1344070.*y.^3 + 2794225.*y.^4 - 3640000.*y.^5 + 2905000.*y.^6 - 1300000.*y.^7 + 250000.*y.^8 + 4536))/126;
                gy = (125.*x.*(50296050.*x.^2.*y.^2 - 129636.*y - 64818.*x - 145600000.*x.^2.*y.^3 - 109200000.*x.^3.*y.^2 + 217875000.*x.^2.*y.^4 + 232400000.*x.^3.*y.^3 + 130725000.*x.^4.*y.^2 - 163800000.*x.^2.*y.^5 - 227500000.*x.^3.*y.^4 - 182000000.*x.^4.*y.^3 - 81900000.*x.^5.*y.^2 + 49000000.*x.^2.*y.^6 + 84000000.*x.^3.*y.^5 + 87500000.*x.^4.*y.^4 + 56000000.*x.^5.*y.^3 + 21000000.*x.^6.*y.^2 + 1580508.*x.*y - 12096630.*x.*y.^2 - 8064420.*x.^2.*y + 44707600.*x.*y.^3 + 22353800.*x.^3.*y - 91000000.*x.*y.^4 - 36400000.*x.^4.*y + 104580000.*x.*y.^5 + 34860000.*x.^5.*y - 63700000.*x.*y.^6 - 18200000.*x.^6.*y + 16000000.*x.*y.^7 + 4000000.*x.^7.*y + 395127.*x.^2 - 1344070.*x.^3 + 2794225.*x.^4 - 3640000.*x.^5 + 2905000.*x.^6 - 1300000.*x.^7 + 250000.*x.^8 + 1185381.*y.^2 - 5376280.*y.^3 + 13971125.*y.^4 - 21840000.*y.^5 + 20335000.*y.^6 - 10400000.*y.^7 + 2250000.*y.^8 + 4536))/126;
            case 32
                gx = -(125.*y.*(6534.*x.*y - 126.*y - 6786.*x - 98490.*x.^2.*y + 676900.*x.^3.*y - 2450000.*x.^4.*y + 4830000.*x.^5.*y - 4900000.*x.^6.*y + 2000000.*x.^7.*y + 108291.*x.^2 - 808220.*x.^3 + 3296125.*x.^4 - 7770000.*x.^5 + 10535000.*x.^6 - 7600000.*x.^7 + 2250000.*x.^8 + 126))/126;
                gy = -(125.*x.*(x + 2.*y - 1).*(3267.*x - 32830.*x.^2 + 169225.*x.^3 - 490000.*x.^4 + 805000.*x.^5 - 700000.*x.^6 + 250000.*x.^7 - 126))/126;
            case 33
                gx = -(125.*y.*(2.*x + y - 1).*(3267.*y - 32830.*y.^2 + 169225.*y.^3 - 490000.*y.^4 + 805000.*y.^5 - 700000.*y.^6 + 250000.*y.^7 - 126))/126;
                gy = -(125.*x.*(6534.*x.*y - 6786.*y - 126.*x - 98490.*x.*y.^2 + 676900.*x.*y.^3 - 2450000.*x.*y.^4 + 4830000.*x.*y.^5 - 4900000.*x.*y.^6 + 2000000.*x.*y.^7 + 108291.*y.^2 - 808220.*y.^3 + 3296125.*y.^4 - 7770000.*y.^5 + 10535000.*y.^6 - 7600000.*y.^7 + 2250000.*y.^8 + 126))/126;
            case 34
                gx = -(250.*y.*(10.*y - 1).*(918.*x + 18.*y - 882.*x.*y + 12180.*x.^2.*y - 73500.*x.^3.*y + 218750.*x.^4.*y - 315000.*x.^5.*y + 175000.*x.^6.*y - 13503.*x.^2 + 89740.*x.^3 - 310625.*x.^4 + 577500.*x.^5 - 542500.*x.^6 + 200000.*x.^7 - 18))/63;
                gy = -(250.*x.*(20.*x.*y - 22.*y - x + 30.*y.^2 + 1).*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18))/63;
            case 35
                gx = -(250.*y.*(50.*y.^2 - 15.*y + 1).*(274.*x.*y - 6.*y - 286.*x - 3375.*x.^2.*y + 17000.*x.^3.*y - 37500.*x.^4.*y + 30000.*x.^5.*y + 3786.*x.^2 - 21500.*x.^3 + 58750.*x.^4 - 75000.*x.^5 + 35000.*x.^6 + 6))/27;
                gy = -(250.*x.*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(x + 32.*y - 30.*x.*y + 150.*x.*y.^2 - 195.*y.^2 + 200.*y.^3 - 1))/27;
            case 36
                gx = -(25.*y.*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(262.*x + 6.*y - 250.*x.*y + 2625.*x.^2.*y - 10000.*x.^3.*y + 12500.*x.^4.*y - 3000.*x.^2 + 13500.*x.^3 - 25000.*x.^4 + 15000.*x.^5 - 6))/9;
                gy = -(25.*x.*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(110.*x.*y - 116.*y - 3.*x - 900.*x.*y.^2 + 2000.*x.*y.^3 + 1065.*y.^2 - 3200.*y.^3 + 2500.*y.^4 + 3))/9;
            case 37
                gx = -(25.*y.*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(110.*x.*y - 3.*y - 116.*x - 900.*x.^2.*y + 2000.*x.^3.*y + 1065.*x.^2 - 3200.*x.^3 + 2500.*x.^4 + 3))/9;
                gy = -(25.*x.*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(6.*x + 262.*y - 250.*x.*y + 2625.*x.*y.^2 - 10000.*x.*y.^3 + 12500.*x.*y.^4 - 3000.*y.^2 + 13500.*y.^3 - 25000.*y.^4 + 15000.*y.^5 - 6))/9;
            case 38
                gx = -(250.*y.*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(32.*x + y - 30.*x.*y + 150.*x.^2.*y - 195.*x.^2 + 200.*x.^3 - 1))/27;
                gy = -(250.*x.*(50.*x.^2 - 15.*x + 1).*(274.*x.*y - 286.*y - 6.*x - 3375.*x.*y.^2 + 17000.*x.*y.^3 - 37500.*x.*y.^4 + 30000.*x.*y.^5 + 3786.*y.^2 - 21500.*y.^3 + 58750.*y.^4 - 75000.*y.^5 + 35000.*y.^6 + 6))/27;
            case 39
                gx = -(250.*y.*(20.*x.*y - y - 22.*x + 30.*x.^2 + 1).*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18))/63;
                gy = -(250.*x.*(10.*x - 1).*(18.*x + 918.*y - 882.*x.*y + 12180.*x.*y.^2 - 73500.*x.*y.^3 + 218750.*x.*y.^4 - 315000.*x.*y.^5 + 175000.*x.*y.^6 - 13503.*y.^2 + 89740.*y.^3 - 310625.*y.^4 + 577500.*y.^5 - 542500.*y.^6 + 200000.*y.^7 - 18))/63;
            case 40
                gx = (250.*y.*(40.*x.*y - 19.*y - 38.*x + 30.*x.^2 + 10.*y.^2 + 9).*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18))/63;
                gy = (250.*x.*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9).*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18))/63 + (250.*x.*y.*(8120.*y - 55125.*y.^2 + 175000.*y.^3 - 262500.*y.^4 + 150000.*y.^5 - 441).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/63 + (250.*x.*y.*(20.*x + 20.*y - 19).*(4060.*y.^2 - 441.*y - 18375.*y.^3 + 43750.*y.^4 - 52500.*y.^5 + 25000.*y.^6 + 18))/63;
            case 41
                gx = -(250.*y.*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(242.*x + 121.*y - 540.*x.*y + 300.*x.*y.^2 + 450.*x.^2.*y - 405.*x.^2 + 200.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/27;
                gy = - (250.*x.*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/27 - (250.*x.*y.*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(300.*x.*y - 270.*y - 270.*x + 150.*x.^2 + 150.*y.^2 + 121))/27 - (250.*x.*y.*(12750.*y.^2 - 2250.*y - 30000.*y.^3 + 25000.*y.^4 + 137).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/27;
            case 42
                gx = (25.*y.*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(9000.*x.^2.*y.^2 - 1207.*y - 2414.*x + 8620.*x.*y - 10200.*x.*y.^2 - 15300.*x.^2.*y + 4000.*x.*y.^3 + 8000.*x.^3.*y + 6465.*x.^2 - 6800.*x.^3 + 2500.*x.^4 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/9;
                gy = (25.*x.*(7623375.*x.^2.*y.^2 - 77484.*y - 7242.*x - 40900000.*x.^2.*y.^3 - 5212500.*x.^3.*y.^2 + 103812500.*x.^2.*y.^4 + 24000000.*x.^3.*y.^3 + 1312500.*x.^4.*y.^2 - 121500000.*x.^2.*y.^5 - 46250000.*x.^3.*y.^4 - 5000000.*x.^4.*y.^3 + 52500000.*x.^2.*y.^6 + 30000000.*x.^3.*y.^5 + 6250000.*x.^4.*y.^4 + 353470.*x.*y - 4876425.*x.*y.^2 - 599950.*x.^2.*y + 29753000.*x.*y.^3 + 449000.*x.^3.*y - 92525000.*x.*y.^4 - 125000.*x.^4.*y + 151650000.*x.*y.^5 - 124250000.*x.*y.^6 + 40000000.*x.*y.^7 + 12930.*x.^2 - 10200.*x.^3 + 3000.*x.^4 + 1152915.*y.^2 - 7862800.*y.^3 + 28743125.*y.^4 - 59730000.*y.^5 + 70525000.*y.^6 - 44000000.*y.^7 + 11250000.*y.^8 + 1512))/9;
            case 43
                gx = -(25.*y.*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(9762.*x + 4881.*y - 180000.*x.^2.*y.^2 + 75000.*x.^2.*y.^3 + 100000.*x.^3.*y.^2 - 50000.*x.*y + 95250.*x.*y.^2 + 142875.*x.^2.*y - 80000.*x.*y.^3 - 160000.*x.^3.*y + 25000.*x.*y.^4 + 62500.*x.^4.*y - 37500.*x.^2 + 63500.*x.^3 - 50000.*x.^4 + 15000.*x.^5 - 12500.*y.^2 + 15875.*y.^3 - 10000.*y.^4 + 2500.*y.^5 - 756))/9;
                gy = -(25.*x.*(19648125.*x.^2.*y.^2 - 112446.*y - 14643.*x - 95650000.*x.^2.*y.^3 - 21112500.*x.^3.*y.^2 + 215937500.*x.^2.*y.^4 + 85250000.*x.^3.*y.^3 + 11062500.*x.^4.*y.^2 - 225000000.*x.^2.*y.^5 - 137500000.*x.^3.*y.^4 - 35000000.*x.^4.*y.^3 - 2250000.*x.^5.*y.^2 + 87500000.*x.^2.*y.^6 + 75000000.*x.^3.*y.^5 + 31250000.*x.^4.*y.^4 + 5000000.*x.^5.*y.^3 + 686910.*x.*y - 8946525.*x.*y.^2 - 1660750.*x.^2.*y + 50719500.*x.*y.^3 + 1986250.*x.^3.*y - 145125000.*x.*y.^4 - 1175000.*x.^4.*y + 219000000.*x.*y.^5 + 275000.*x.^5.*y - 166250000.*x.*y.^6 + 50000000.*x.*y.^7 + 37500.*x.^2 - 47625.*x.^3 + 30000.*x.^4 - 7500.*x.^5 + 1598265.*y.^2 - 10309700.*y.^3 + 35468125.*y.^4 - 69420000.*y.^5 + 77525000.*y.^6 - 46000000.*y.^7 + 11250000.*y.^8 + 2268))/9;
            case 44
                gx = (250.*y.*(50.*y.^2 - 15.*y + 1).*(751500.*x.^2.*y.^2 - 6393.*y - 12786.*x - 675000.*x.^2.*y.^3 - 900000.*x.^3.*y.^2 + 225000.*x.^2.*y.^4 + 400000.*x.^3.*y.^3 + 375000.*x.^4.*y.^2 + 89048.*x.*y - 245250.*x.*y.^2 - 367875.*x.^2.*y + 334000.*x.*y.^3 + 668000.*x.^3.*y - 225000.*x.*y.^4 - 562500.*x.^4.*y + 60000.*x.*y.^5 + 180000.*x.^5.*y + 66786.*x.^2 - 163500.*x.^3 + 208750.*x.^4 - 135000.*x.^5 + 35000.*x.^6 + 22262.*y.^2 - 40875.*y.^3 + 41750.*y.^4 - 22500.*y.^5 + 5000.*y.^6 + 756))/27;
                gy = (250.*x.*(9608925.*x.^2.*y.^2 - 35466.*y - 6393.*x - 40455000.*x.^2.*y.^3 - 14321250.*x.^3.*y.^2 + 79875000.*x.^2.*y.^4 + 47300000.*x.^3.*y.^3 + 11550000.*x.^4.*y.^2 - 74250000.*x.^2.*y.^5 - 63750000.*x.^3.*y.^4 - 27000000.*x.^4.*y.^3 - 4725000.*x.^5.*y.^2 + 26250000.*x.^2.*y.^6 + 30000000.*x.^3.*y.^5 + 18750000.*x.^4.*y.^4 + 6000000.*x.^5.*y.^3 + 750000.*x.^6.*y.^2 + 280838.*x.*y - 3330405.*x.*y.^2 - 913110.*x.^2.*y + 16930300.*x.*y.^3 + 1560250.*x.^3.*y - 43743750.*x.*y.^4 - 1477500.*x.^4.*y + 60405000.*x.*y.^5 + 735000.*x.^5.*y - 42525000.*x.*y.^6 - 150000.*x.^6.*y + 12000000.*x.*y.^7 + 22262.*x.^2 - 40875.*x.^3 + 41750.*x.^4 - 22500.*x.^5 + 5000.*x.^6 + 467871.*y.^2 - 2777820.*y.^3 + 8839875.*y.^4 - 16155000.*y.^5 + 17010000.*y.^6 - 9600000.*y.^7 + 2250000.*y.^8 + 756))/27;
            case 45
                gx = -(250.*y.*(10.*y - 1).*(33132.*x + 16566.*y - 5181750.*x.^2.*y.^2 + 7612500.*x.^2.*y.^3 + 10150000.*x.^3.*y.^2 - 5512500.*x.^2.*y.^4 - 9800000.*x.^3.*y.^3 - 9187500.*x.^4.*y.^2 + 1575000.*x.^2.*y.^5 + 3500000.*x.^3.*y.^4 + 4375000.*x.^4.*y.^3 + 3150000.*x.^5.*y.^2 - 305956.*x.*y + 1158360.*x.*y.^2 + 1737540.*x.^2.*y - 2303000.*x.*y.^3 - 4606000.*x.^3.*y + 2537500.*x.*y.^4 + 6343750.*x.^4.*y - 1470000.*x.*y.^5 - 4410000.*x.^5.*y + 350000.*x.*y.^6 + 1225000.*x.^6.*y - 229467.*x.^2 + 772240.*x.^3 - 1439375.*x.^4 + 1522500.*x.^5 - 857500.*x.^6 + 200000.*x.^7 - 76489.*y.^2 + 193060.*y.^3 - 287875.*y.^4 + 253750.*y.^5 - 122500.*y.^6 + 25000.*y.^7 - 1512))/63;
                gy = -(250.*x.*(22557150.*x.^2.*y.^2 - 63372.*y - 16566.*x - 79240000.*x.^2.*y.^3 - 42157500.*x.^3.*y.^2 + 136062500.*x.^2.*y.^4 + 111300000.*x.^3.*y.^3 + 43575000.*x.^4.*y.^2 - 113400000.*x.^2.*y.^5 - 126875000.*x.^3.*y.^4 - 77000000.*x.^4.*y.^3 - 23625000.*x.^5.*y.^2 + 36750000.*x.^2.*y.^6 + 52500000.*x.^3.*y.^5 + 43750000.*x.^4.*y.^4 + 21000000.*x.^5.*y.^3 + 5250000.*x.^6.*y.^2 + 637276.*x.*y - 6326880.*x.*y.^2 - 2688140.*x.^2.*y + 27773200.*x.*y.^3 + 6164200.*x.^3.*y - 63918750.*x.*y.^4 - 8295000.*x.^4.*y + 80535000.*x.*y.^5 + 6545000.*x.^5.*y - 52675000.*x.*y.^6 - 2800000.*x.^6.*y + 14000000.*x.*y.^7 + 500000.*x.^7.*y + 76489.*x.^2 - 193060.*x.^3 + 287875.*x.^4 - 253750.*x.^5 + 122500.*x.^6 - 25000.*x.^7 + 726447.*y.^2 - 3831800.*y.^3 + 11092375.*y.^4 - 18795000.*y.^5 + 18620000.*y.^6 - 10000000.*y.^7 + 2250000.*y.^8 + 1512))/63;
            case 46
                gx = -(250.*y.*(22557150.*x.^2.*y.^2 - 16566.*y - 63372.*x - 42157500.*x.^2.*y.^3 - 79240000.*x.^3.*y.^2 + 43575000.*x.^2.*y.^4 + 111300000.*x.^3.*y.^3 + 136062500.*x.^4.*y.^2 - 23625000.*x.^2.*y.^5 - 77000000.*x.^3.*y.^4 - 126875000.*x.^4.*y.^3 - 113400000.*x.^5.*y.^2 + 5250000.*x.^2.*y.^6 + 21000000.*x.^3.*y.^5 + 43750000.*x.^4.*y.^4 + 52500000.*x.^5.*y.^3 + 36750000.*x.^6.*y.^2 + 637276.*x.*y - 2688140.*x.*y.^2 - 6326880.*x.^2.*y + 6164200.*x.*y.^3 + 27773200.*x.^3.*y - 8295000.*x.*y.^4 - 63918750.*x.^4.*y + 6545000.*x.*y.^5 + 80535000.*x.^5.*y - 2800000.*x.*y.^6 - 52675000.*x.^6.*y + 500000.*x.*y.^7 + 14000000.*x.^7.*y + 726447.*x.^2 - 3831800.*x.^3 + 11092375.*x.^4 - 18795000.*x.^5 + 18620000.*x.^6 - 10000000.*x.^7 + 2250000.*x.^8 + 76489.*y.^2 - 193060.*y.^3 + 287875.*y.^4 - 253750.*y.^5 + 122500.*y.^6 - 25000.*y.^7 + 1512))/63;
                gy = -(250.*x.*(10.*x - 1).*(16566.*x + 33132.*y - 5181750.*x.^2.*y.^2 + 10150000.*x.^2.*y.^3 + 7612500.*x.^3.*y.^2 - 9187500.*x.^2.*y.^4 - 9800000.*x.^3.*y.^3 - 5512500.*x.^4.*y.^2 + 3150000.*x.^2.*y.^5 + 4375000.*x.^3.*y.^4 + 3500000.*x.^4.*y.^3 + 1575000.*x.^5.*y.^2 - 305956.*x.*y + 1737540.*x.*y.^2 + 1158360.*x.^2.*y - 4606000.*x.*y.^3 - 2303000.*x.^3.*y + 6343750.*x.*y.^4 + 2537500.*x.^4.*y - 4410000.*x.*y.^5 - 1470000.*x.^5.*y + 1225000.*x.*y.^6 + 350000.*x.^6.*y - 76489.*x.^2 + 193060.*x.^3 - 287875.*x.^4 + 253750.*x.^5 - 122500.*x.^6 + 25000.*x.^7 - 229467.*y.^2 + 772240.*y.^3 - 1439375.*y.^4 + 1522500.*y.^5 - 857500.*y.^6 + 200000.*y.^7 - 1512))/63;
            case 47
                gx = (250.*y.*(9608925.*x.^2.*y.^2 - 6393.*y - 35466.*x - 14321250.*x.^2.*y.^3 - 40455000.*x.^3.*y.^2 + 11550000.*x.^2.*y.^4 + 47300000.*x.^3.*y.^3 + 79875000.*x.^4.*y.^2 - 4725000.*x.^2.*y.^5 - 27000000.*x.^3.*y.^4 - 63750000.*x.^4.*y.^3 - 74250000.*x.^5.*y.^2 + 750000.*x.^2.*y.^6 + 6000000.*x.^3.*y.^5 + 18750000.*x.^4.*y.^4 + 30000000.*x.^5.*y.^3 + 26250000.*x.^6.*y.^2 + 280838.*x.*y - 913110.*x.*y.^2 - 3330405.*x.^2.*y + 1560250.*x.*y.^3 + 16930300.*x.^3.*y - 1477500.*x.*y.^4 - 43743750.*x.^4.*y + 735000.*x.*y.^5 + 60405000.*x.^5.*y - 150000.*x.*y.^6 - 42525000.*x.^6.*y + 12000000.*x.^7.*y + 467871.*x.^2 - 2777820.*x.^3 + 8839875.*x.^4 - 16155000.*x.^5 + 17010000.*x.^6 - 9600000.*x.^7 + 2250000.*x.^8 + 22262.*y.^2 - 40875.*y.^3 + 41750.*y.^4 - 22500.*y.^5 + 5000.*y.^6 + 756))/27;
                gy = (250.*x.*(50.*x.^2 - 15.*x + 1).*(751500.*x.^2.*y.^2 - 12786.*y - 6393.*x - 900000.*x.^2.*y.^3 - 675000.*x.^3.*y.^2 + 375000.*x.^2.*y.^4 + 400000.*x.^3.*y.^3 + 225000.*x.^4.*y.^2 + 89048.*x.*y - 367875.*x.*y.^2 - 245250.*x.^2.*y + 668000.*x.*y.^3 + 334000.*x.^3.*y - 562500.*x.*y.^4 - 225000.*x.^4.*y + 180000.*x.*y.^5 + 60000.*x.^5.*y + 22262.*x.^2 - 40875.*x.^3 + 41750.*x.^4 - 22500.*x.^5 + 5000.*x.^6 + 66786.*y.^2 - 163500.*y.^3 + 208750.*y.^4 - 135000.*y.^5 + 35000.*y.^6 + 756))/27;
            case 48
                gx = -(25.*y.*(19648125.*x.^2.*y.^2 - 14643.*y - 112446.*x - 21112500.*x.^2.*y.^3 - 95650000.*x.^3.*y.^2 + 11062500.*x.^2.*y.^4 + 85250000.*x.^3.*y.^3 + 215937500.*x.^4.*y.^2 - 2250000.*x.^2.*y.^5 - 35000000.*x.^3.*y.^4 - 137500000.*x.^4.*y.^3 - 225000000.*x.^5.*y.^2 + 5000000.*x.^3.*y.^5 + 31250000.*x.^4.*y.^4 + 75000000.*x.^5.*y.^3 + 87500000.*x.^6.*y.^2 + 686910.*x.*y - 1660750.*x.*y.^2 - 8946525.*x.^2.*y + 1986250.*x.*y.^3 + 50719500.*x.^3.*y - 1175000.*x.*y.^4 - 145125000.*x.^4.*y + 275000.*x.*y.^5 + 219000000.*x.^5.*y - 166250000.*x.^6.*y + 50000000.*x.^7.*y + 1598265.*x.^2 - 10309700.*x.^3 + 35468125.*x.^4 - 69420000.*x.^5 + 77525000.*x.^6 - 46000000.*x.^7 + 11250000.*x.^8 + 37500.*y.^2 - 47625.*y.^3 + 30000.*y.^4 - 7500.*y.^5 + 2268))/9;
                gy = -(25.*x.*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(4881.*x + 9762.*y - 180000.*x.^2.*y.^2 + 100000.*x.^2.*y.^3 + 75000.*x.^3.*y.^2 - 50000.*x.*y + 142875.*x.*y.^2 + 95250.*x.^2.*y - 160000.*x.*y.^3 - 80000.*x.^3.*y + 62500.*x.*y.^4 + 25000.*x.^4.*y - 12500.*x.^2 + 15875.*x.^3 - 10000.*x.^4 + 2500.*x.^5 - 37500.*y.^2 + 63500.*y.^3 - 50000.*y.^4 + 15000.*y.^5 - 756))/9;
            case 49
                gx = (25.*y.*(7623375.*x.^2.*y.^2 - 7242.*y - 77484.*x - 5212500.*x.^2.*y.^3 - 40900000.*x.^3.*y.^2 + 1312500.*x.^2.*y.^4 + 24000000.*x.^3.*y.^3 + 103812500.*x.^4.*y.^2 - 5000000.*x.^3.*y.^4 - 46250000.*x.^4.*y.^3 - 121500000.*x.^5.*y.^2 + 6250000.*x.^4.*y.^4 + 30000000.*x.^5.*y.^3 + 52500000.*x.^6.*y.^2 + 353470.*x.*y - 599950.*x.*y.^2 - 4876425.*x.^2.*y + 449000.*x.*y.^3 + 29753000.*x.^3.*y - 125000.*x.*y.^4 - 92525000.*x.^4.*y + 151650000.*x.^5.*y - 124250000.*x.^6.*y + 40000000.*x.^7.*y + 1152915.*x.^2 - 7862800.*x.^3 + 28743125.*x.^4 - 59730000.*x.^5 + 70525000.*x.^6 - 44000000.*x.^7 + 11250000.*x.^8 + 12930.*y.^2 - 10200.*y.^3 + 3000.*y.^4 + 1512))/9;
                gy = (25.*x.*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(9000.*x.^2.*y.^2 - 2414.*y - 1207.*x + 8620.*x.*y - 15300.*x.*y.^2 - 10200.*x.^2.*y + 8000.*x.*y.^3 + 4000.*x.^3.*y + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 6465.*y.^2 - 6800.*y.^3 + 2500.*y.^4 + 252))/9;
            case 50
                gx = - (250.*y.*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/27 - (250.*x.*y.*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(300.*x.*y - 270.*y - 270.*x + 150.*x.^2 + 150.*y.^2 + 121))/27 - (250.*x.*y.*(12750.*x.^2 - 2250.*x - 30000.*x.^3 + 25000.*x.^4 + 137).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/27;
                gy = -(250.*x.*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(121.*x + 242.*y - 540.*x.*y + 450.*x.*y.^2 + 300.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 405.*y.^2 + 200.*y.^3 - 36))/27;
            case 51
                gx = (250.*y.*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9).*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18))/63 + (250.*x.*y.*(8120.*x - 55125.*x.^2 + 175000.*x.^3 - 262500.*x.^4 + 150000.*x.^5 - 441).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/63 + (250.*x.*y.*(20.*x + 20.*y - 19).*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18))/63;
                gy = (250.*x.*(40.*x.*y - 38.*y - 19.*x + 10.*x.^2 + 30.*y.^2 + 9).*(4060.*x.^2 - 441.*x - 18375.*x.^3 + 43750.*x.^4 - 52500.*x.^5 + 25000.*x.^6 + 18))/63;
            case 52
                gx = (125.*y.*(10.*y - 1).*(27906.*x + 6393.*y - 4430250.*x.^2.*y.^2 + 5685000.*x.^2.*y.^3 + 10920000.*x.^3.*y.^2 - 3600000.*x.^2.*y.^4 - 9400000.*x.^3.*y.^3 - 11625000.*x.^4.*y.^2 + 900000.*x.^2.*y.^5 + 3000000.*x.^3.*y.^4 + 5000000.*x.^4.*y.^3 + 4500000.*x.^5.*y.^2 - 216908.*x.*y + 690490.*x.*y.^2 + 1703595.*x.^2.*y - 1151500.*x.*y.^3 - 5573000.*x.^3.*y + 1060000.*x.*y.^4 + 8912500.*x.^4.*y - 510000.*x.*y.^5 - 6930000.*x.^5.*y + 100000.*x.*y.^6 + 2100000.*x.^6.*y - 258576.*x.^2 + 1053980.*x.^3 - 2252500.*x.^4 + 2640000.*x.^5 - 1610000.*x.^6 + 400000.*x.^7 - 22262.*y.^2 + 40875.*y.^3 - 41750.*y.^4 + 22500.*y.^5 - 5000.*y.^6 - 756))/18;
                gy = (125.*x.*(10.*x - 1).*(6393.*x + 27906.*y - 4430250.*x.^2.*y.^2 + 10920000.*x.^2.*y.^3 + 5685000.*x.^3.*y.^2 - 11625000.*x.^2.*y.^4 - 9400000.*x.^3.*y.^3 - 3600000.*x.^4.*y.^2 + 4500000.*x.^2.*y.^5 + 5000000.*x.^3.*y.^4 + 3000000.*x.^4.*y.^3 + 900000.*x.^5.*y.^2 - 216908.*x.*y + 1703595.*x.*y.^2 + 690490.*x.^2.*y - 5573000.*x.*y.^3 - 1151500.*x.^3.*y + 8912500.*x.*y.^4 + 1060000.*x.^4.*y - 6930000.*x.*y.^5 - 510000.*x.^5.*y + 2100000.*x.*y.^6 + 100000.*x.^6.*y - 22262.*x.^2 + 40875.*x.^3 - 41750.*x.^4 + 22500.*x.^5 - 5000.*x.^6 - 258576.*y.^2 + 1053980.*y.^3 - 2252500.*y.^4 + 2640000.*y.^5 - 1610000.*y.^6 + 400000.*y.^7 - 756))/18;
            case 53
                gx = (125.*y.*(10.*y - 1).*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/18 + (125.*x.*y.*(10.*y - 1).*(12750.*x.^2 - 2250.*x - 30000.*x.^3 + 25000.*x.^4 + 137).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/18 + (125.*x.*y.*(10.*y - 1).*(20.*x + 20.*y - 19).*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6))/18;
                gy = (125.*x.*(137.*x - 1125.*x.^2 + 4250.*x.^3 - 7500.*x.^4 + 5000.*x.^5 - 6).*(19.*x + 218.*y - 420.*x.*y + 600.*x.*y.^2 + 200.*x.^2.*y - 10.*x.^2 - 600.*y.^2 + 400.*y.^3 - 9))/18;
            case 54
                gx = (125.*y.*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(218.*x + 19.*y - 420.*x.*y + 200.*x.*y.^2 + 600.*x.^2.*y - 600.*x.^2 + 400.*x.^3 - 10.*y.^2 - 9))/18;
                gy = (125.*x.*(10.*x - 1).*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/18 + (125.*x.*y.*(10.*x - 1).*(12750.*y.^2 - 2250.*y - 30000.*y.^3 + 25000.*y.^4 + 137).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/18 + (125.*x.*y.*(10.*x - 1).*(20.*x + 20.*y - 19).*(137.*y - 1125.*y.^2 + 4250.*y.^3 - 7500.*y.^4 + 5000.*y.^5 - 6))/18;
            case 55
                gx = (50.*y.*(50.*y.^2 - 15.*y + 1).*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/9 + (50.*x.*y.*(50.*y.^2 - 15.*y + 1).*(1750.*x - 7500.*x.^2 + 10000.*x.^3 - 125).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/9 + (50.*x.*y.*(20.*x + 20.*y - 19).*(50.*y.^2 - 15.*y + 1).*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6))/9;
                gy = (50.*x.*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(1500.*x.^2.*y.^2 - 308.*y - 19.*x + 610.*x.*y - 3750.*x.*y.^2 - 300.*x.^2.*y + 4000.*x.*y.^3 + 10.*x.^2 + 2235.*y.^2 - 4400.*y.^3 + 2500.*y.^4 + 9))/9;
            case 56
                gx = (125.*y.*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(1104.*x + 57.*y - 9000.*x.^2.*y.^2 + 20000.*x.^3.*y.^2 - 2210.*x.*y + 1100.*x.*y.^2 + 20400.*x.^2.*y - 62000.*x.^3.*y + 50000.*x.^4.*y - 11325.*x.^2 + 43000.*x.^3 - 62500.*x.^4 + 30000.*x.^5 - 30.*y.^2 - 27))/36;
                gy = (125.*x.*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(57.*x + 1104.*y - 9000.*x.^2.*y.^2 + 20000.*x.^2.*y.^3 - 2210.*x.*y + 20400.*x.*y.^2 + 1100.*x.^2.*y - 62000.*x.*y.^3 + 50000.*x.*y.^4 - 30.*x.^2 - 11325.*y.^2 + 43000.*y.^3 - 62500.*y.^4 + 30000.*y.^5 - 27))/36;
            case 57
                gx = (50.*y.*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(1500.*x.^2.*y.^2 - 19.*y - 308.*x + 610.*x.*y - 300.*x.*y.^2 - 3750.*x.^2.*y + 4000.*x.^3.*y + 2235.*x.^2 - 4400.*x.^3 + 2500.*x.^4 + 10.*y.^2 + 9))/9;
                gy = (50.*x.*(50.*x.^2 - 15.*x + 1).*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/9 + (50.*x.*y.*(50.*x.^2 - 15.*x + 1).*(1750.*y - 7500.*y.^2 + 10000.*y.^3 - 125).*(20.*x.*y - 19.*y - 19.*x + 10.*x.^2 + 10.*y.^2 + 9))/9 + (50.*x.*y.*(20.*x + 20.*y - 19).*(50.*x.^2 - 15.*x + 1).*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6))/9;
            case 58
                gx = -(50.*y.*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(4500.*x.^2.*y.^2 - 121.*y - 962.*x + 2960.*x.*y - 3000.*x.*y.^2 - 8550.*x.^2.*y + 1000.*x.*y.^3 + 6000.*x.^3.*y + 4035.*x.^2 - 5600.*x.^3 + 2500.*x.^4 + 135.*y.^2 - 50.*y.^3 + 36))/9;
                gy = - (50.*x.*(10.*x - 1).*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/9 - (50.*x.*y.*(10.*x - 1).*(875.*y.^2 - 125.*y - 2500.*y.^3 + 2500.*y.^4 + 6).*(300.*x.*y - 270.*y - 270.*x + 150.*x.^2 + 150.*y.^2 + 121))/9 - (50.*x.*y.*(10.*x - 1).*(1750.*y - 7500.*y.^2 + 10000.*y.^3 - 125).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/9;
            case 59
                gx = (125.*y.*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(7454.*x + 1207.*y - 162000.*x.^2.*y.^2 + 60000.*x.^2.*y.^3 + 120000.*x.^3.*y.^2 - 32760.*x.*y + 53300.*x.*y.^2 + 144600.*x.^2.*y - 38000.*x.*y.^3 - 212000.*x.^3.*y + 10000.*x.*y.^4 + 100000.*x.^4.*y - 42675.*x.^2 + 93000.*x.^3 - 87500.*x.^4 + 30000.*x.^5 - 2155.*y.^2 + 1700.*y.^3 - 500.*y.^4 - 252))/36;
                gy = (125.*x.*(10.*x - 1).*(3621.*x + 34962.*y - 2808000.*x.^2.*y.^2 + 11090000.*x.^2.*y.^3 + 1860000.*x.^3.*y.^2 - 17250000.*x.^2.*y.^4 - 5800000.*x.^3.*y.^3 - 450000.*x.^4.*y.^2 + 9000000.*x.^2.*y.^5 + 5000000.*x.^3.*y.^4 + 1000000.*x.^4.*y.^3 - 158630.*x.*y + 1843350.*x.*y.^2 + 267650.*x.^2.*y - 8732000.*x.*y.^3 - 199000.*x.^3.*y + 18975000.*x.*y.^4 + 55000.*x.^4.*y - 18900000.*x.*y.^5 + 7000000.*x.*y.^6 - 6465.*x.^2 + 5100.*x.^3 - 1500.*x.^4 - 445350.*y.^2 + 2446900.*y.^3 - 6725000.*y.^4 + 9690000.*y.^5 - 7000000.*y.^6 + 2000000.*y.^7 - 756))/36;
            case 60
                gx = -(50.*y.*(50.*y.^2 - 15.*y + 1).*(1608750.*x.^2.*y.^2 - 4881.*y - 24882.*x - 1275000.*x.^2.*y.^3 - 2500000.*x.^3.*y.^2 + 375000.*x.^2.*y.^4 + 1000000.*x.^3.*y.^3 + 1250000.*x.^4.*y.^2 + 147620.*x.*y - 345250.*x.*y.^2 - 892875.*x.^2.*y + 397500.*x.*y.^3 + 2065000.*x.^3.*y - 225000.*x.*y.^4 - 2062500.*x.^4.*y + 50000.*x.*y.^5 + 750000.*x.^5.*y + 183930.*x.^2 - 563500.*x.^3 + 843750.*x.^4 - 615000.*x.^5 + 175000.*x.^6 + 12500.*y.^2 - 15875.*y.^3 + 10000.*y.^4 - 2500.*y.^5 + 756))/9;
                gy = -(50.*x.*(10.*x - 1).*(4881.*x + 32442.*y - 4198125.*x.^2.*y.^2 + 13225000.*x.^2.*y.^3 + 4256250.*x.^3.*y.^2 - 16875000.*x.^2.*y.^4 - 9500000.*x.^3.*y.^3 - 2062500.*x.^4.*y.^2 + 7500000.*x.^2.*y.^5 + 6250000.*x.^3.*y.^4 + 2500000.*x.^4.*y.^3 + 375000.*x.^5.*y.^2 - 196430.*x.*y + 2000025.*x.*y.^2 + 470250.*x.^2.*y - 8017500.*x.*y.^3 - 556250.*x.^3.*y + 14968750.*x.*y.^4 + 325000.*x.^4.*y - 13125000.*x.*y.^5 - 75000.*x.^5.*y + 4375000.*x.*y.^6 - 12500.*x.^2 + 15875.*x.^3 - 10000.*x.^4 + 2500.*x.^5 - 370545.*y.^2 + 1789700.*y.^3 - 4365625.*y.^4 + 5677500.*y.^5 - 3762500.*y.^6 + 1000000.*y.^7 - 756))/9;
            case 61
                gx = -(50.*y.*(10.*y - 1).*(32442.*x + 4881.*y - 4198125.*x.^2.*y.^2 + 4256250.*x.^2.*y.^3 + 13225000.*x.^3.*y.^2 - 2062500.*x.^2.*y.^4 - 9500000.*x.^3.*y.^3 - 16875000.*x.^4.*y.^2 + 375000.*x.^2.*y.^5 + 2500000.*x.^3.*y.^4 + 6250000.*x.^4.*y.^3 + 7500000.*x.^5.*y.^2 - 196430.*x.*y + 470250.*x.*y.^2 + 2000025.*x.^2.*y - 556250.*x.*y.^3 - 8017500.*x.^3.*y + 325000.*x.*y.^4 + 14968750.*x.^4.*y - 75000.*x.*y.^5 - 13125000.*x.^5.*y + 4375000.*x.^6.*y - 370545.*x.^2 + 1789700.*x.^3 - 4365625.*x.^4 + 5677500.*x.^5 - 3762500.*x.^6 + 1000000.*x.^7 - 12500.*y.^2 + 15875.*y.^3 - 10000.*y.^4 + 2500.*y.^5 - 756))/9;
                gy = -(50.*x.*(50.*x.^2 - 15.*x + 1).*(1608750.*x.^2.*y.^2 - 24882.*y - 4881.*x - 2500000.*x.^2.*y.^3 - 1275000.*x.^3.*y.^2 + 1250000.*x.^2.*y.^4 + 1000000.*x.^3.*y.^3 + 375000.*x.^4.*y.^2 + 147620.*x.*y - 892875.*x.*y.^2 - 345250.*x.^2.*y + 2065000.*x.*y.^3 + 397500.*x.^3.*y - 2062500.*x.*y.^4 - 225000.*x.^4.*y + 750000.*x.*y.^5 + 50000.*x.^5.*y + 12500.*x.^2 - 15875.*x.^3 + 10000.*x.^4 - 2500.*x.^5 + 183930.*y.^2 - 563500.*y.^3 + 843750.*y.^4 - 615000.*y.^5 + 175000.*y.^6 + 756))/9;
            case 62
                gx = (125.*y.*(10.*y - 1).*(34962.*x + 3621.*y - 2808000.*x.^2.*y.^2 + 1860000.*x.^2.*y.^3 + 11090000.*x.^3.*y.^2 - 450000.*x.^2.*y.^4 - 5800000.*x.^3.*y.^3 - 17250000.*x.^4.*y.^2 + 1000000.*x.^3.*y.^4 + 5000000.*x.^4.*y.^3 + 9000000.*x.^5.*y.^2 - 158630.*x.*y + 267650.*x.*y.^2 + 1843350.*x.^2.*y - 199000.*x.*y.^3 - 8732000.*x.^3.*y + 55000.*x.*y.^4 + 18975000.*x.^4.*y - 18900000.*x.^5.*y + 7000000.*x.^6.*y - 445350.*x.^2 + 2446900.*x.^3 - 6725000.*x.^4 + 9690000.*x.^5 - 7000000.*x.^6 + 2000000.*x.^7 - 6465.*y.^2 + 5100.*y.^3 - 1500.*y.^4 - 756))/36;
                gy = (125.*x.*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(1207.*x + 7454.*y - 162000.*x.^2.*y.^2 + 120000.*x.^2.*y.^3 + 60000.*x.^3.*y.^2 - 32760.*x.*y + 144600.*x.*y.^2 + 53300.*x.^2.*y - 212000.*x.*y.^3 - 38000.*x.^3.*y + 100000.*x.*y.^4 + 10000.*x.^4.*y - 2155.*x.^2 + 1700.*x.^3 - 500.*x.^4 - 42675.*y.^2 + 93000.*y.^3 - 87500.*y.^4 + 30000.*y.^5 - 252))/36;
            case 63
                gx = - (50.*y.*(10.*y - 1).*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/9 - (50.*x.*y.*(10.*y - 1).*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(300.*x.*y - 270.*y - 270.*x + 150.*x.^2 + 150.*y.^2 + 121))/9 - (50.*x.*y.*(10.*y - 1).*(1750.*x - 7500.*x.^2 + 10000.*x.^3 - 125).*(121.*x + 121.*y - 270.*x.*y + 150.*x.*y.^2 + 150.*x.^2.*y - 135.*x.^2 + 50.*x.^3 - 135.*y.^2 + 50.*y.^3 - 36))/9;
                gy = -(50.*x.*(875.*x.^2 - 125.*x - 2500.*x.^3 + 2500.*x.^4 + 6).*(4500.*x.^2.*y.^2 - 962.*y - 121.*x + 2960.*x.*y - 8550.*x.*y.^2 - 3000.*x.^2.*y + 6000.*x.*y.^3 + 1000.*x.^3.*y + 135.*x.^2 - 50.*x.^3 + 4035.*y.^2 - 5600.*y.^3 + 2500.*y.^4 + 36))/9;
            case 64
                gx = (250.*y.*(50.*y.^2 - 15.*y + 1).*(561750.*x.^2.*y.^2 - 1207.*y - 9974.*x - 345000.*x.^2.*y.^3 - 1200000.*x.^3.*y.^2 + 75000.*x.^2.*y.^4 + 400000.*x.^3.*y.^3 + 750000.*x.^4.*y.^2 + 44830.*x.*y - 74850.*x.*y.^2 - 390300.*x.^2.*y + 55000.*x.*y.^3 + 1176000.*x.^3.*y - 15000.*x.*y.^4 - 1425000.*x.^4.*y + 600000.*x.^5.*y + 98580.*x.^2 - 377500.*x.^3 + 668750.*x.^4 - 555000.*x.^5 + 175000.*x.^6 + 2155.*y.^2 - 1700.*y.^3 + 500.*y.^4 + 252))/27;
                gy = (250.*x.*(50.*x.^2 - 15.*x + 1).*(561750.*x.^2.*y.^2 - 9974.*y - 1207.*x - 1200000.*x.^2.*y.^3 - 345000.*x.^3.*y.^2 + 750000.*x.^2.*y.^4 + 400000.*x.^3.*y.^3 + 75000.*x.^4.*y.^2 + 44830.*x.*y - 390300.*x.*y.^2 - 74850.*x.^2.*y + 1176000.*x.*y.^3 + 55000.*x.^3.*y - 1425000.*x.*y.^4 - 15000.*x.^4.*y + 600000.*x.*y.^5 + 2155.*x.^2 - 1700.*x.^3 + 500.*x.^4 + 98580.*y.^2 - 377500.*y.^3 + 668750.*y.^4 - 555000.*y.^5 + 175000.*y.^6 + 252))/27;
            case 65
                gx = -(250.*y.*(50.*y.^2 - 15.*y + 1).*(146250.*x.^2.*y.^2 - 363.*y - 4686.*x - 45000.*x.^2.*y.^3 - 450000.*x.^3.*y.^2 + 100000.*x.^3.*y.^3 + 375000.*x.^4.*y.^2 + 14930.*x.*y - 15750.*x.*y.^2 - 154800.*x.^2.*y + 5500.*x.*y.^3 + 599000.*x.^3.*y - 900000.*x.^4.*y + 450000.*x.^5.*y + 53580.*x.^2 - 247500.*x.^3 + 518750.*x.^4 - 495000.*x.^5 + 175000.*x.^6 + 405.*y.^2 - 150.*y.^3 + 108))/27;
                gy = -(250.*x.*(55.*x - 300.*x.^2 + 500.*x.^3 - 3).*(121.*x + 1322.*y - 27000.*x.^2.*y.^2 + 30000.*x.^2.*y.^3 + 7500.*x.^3.*y.^2 - 4170.*x.*y + 30750.*x.*y.^2 + 4350.*x.^2.*y - 63000.*x.*y.^3 - 1500.*x.^3.*y + 37500.*x.*y.^4 - 135.*x.^2 + 50.*x.^3 - 11250.*y.^2 + 32500.*y.^3 - 37500.*y.^4 + 15000.*y.^5 - 36))/27;
            case 66
                gx = -(250.*y.*(55.*y - 300.*y.^2 + 500.*y.^3 - 3).*(1322.*x + 121.*y - 27000.*x.^2.*y.^2 + 7500.*x.^2.*y.^3 + 30000.*x.^3.*y.^2 - 4170.*x.*y + 4350.*x.*y.^2 + 30750.*x.^2.*y - 1500.*x.*y.^3 - 63000.*x.^3.*y + 37500.*x.^4.*y - 11250.*x.^2 + 32500.*x.^3 - 37500.*x.^4 + 15000.*x.^5 - 135.*y.^2 + 50.*y.^3 - 36))/27;
                gy = -(250.*x.*(50.*x.^2 - 15.*x + 1).*(146250.*x.^2.*y.^2 - 4686.*y - 363.*x - 450000.*x.^2.*y.^3 - 45000.*x.^3.*y.^2 + 375000.*x.^2.*y.^4 + 100000.*x.^3.*y.^3 + 14930.*x.*y - 154800.*x.*y.^2 - 15750.*x.^2.*y + 599000.*x.*y.^3 + 5500.*x.^3.*y - 900000.*x.*y.^4 + 450000.*x.*y.^5 + 405.*x.^2 - 150.*x.^3 + 53580.*y.^2 - 247500.*y.^3 + 518750.*y.^4 - 495000.*y.^5 + 175000.*y.^6 + 108))/27;
        end
end

end

function [NQ,xq,yq,wq]=quadratura(fdq)

switch fdq
    case {'degree=1','degree=2','degree=0'} % N = 3:
        q = [...
            0.1666666666667   0.6666666666667  0.6666666666667
            0.6666666666667   0.1666666666667  0.6666666666667
            0.1666666666667   0.1666666666667  0.6666666666667];
    case {'degree=3','degree=4'} % N = 6:
        q = [...
            0.0915762135098   0.0915762135098  0.2199034873106
            0.8168475729805   0.0915762135098  0.2199034873106
            0.0915762135098   0.8168475729805  0.2199034873106
            0.1081030181681   0.4459484909160  0.4467631793560
            0.4459484909160   0.1081030181681  0.4467631793560
            0.4459484909160   0.4459484909160  0.4467631793560];
    case {'degree=5','degree=6'} % N = 10:
        q = [...
            0.0000000000000   1.0000000000000  0.0262712099504
            1.0000000000000   0.0000000000000  0.0262716612068
            0.0000000000000   0.0000000000000  0.0274163947600
            0.2673273531185   0.6728199218710  0.2348383865823
            0.6728175529461   0.2673288599482  0.2348412238268
            0.0649236350054   0.6716530111494  0.2480251793114
            0.6716498539042   0.0649251690029  0.2480304922521
            0.0654032456800   0.2693789366453  0.2518604605529
            0.2693767069140   0.0654054874919  0.2518660533658
            0.3386738503896   0.3386799893027  0.4505789381914];
    case 'degree=7' % N = 15:
        q = [...
            1.0000000000000   0.0000000000000  0.0102558174092
            0.0000000000000   0.0000000000000  0.0102558174092
            0.0000000000000   1.0000000000000  0.0102558174092
            0.7839656651012   0.0421382841642  0.1116047046647
            0.1738960507345   0.7839656651012  0.1116047046647
            0.1738960507345   0.0421382841642  0.1116047046647
            0.0421382841642   0.1738960507345  0.1116047046647
            0.7839656651012   0.1738960507345  0.1116047046647
            0.0421382841642   0.7839656651012  0.1116047046647
            0.4743880861752   0.4743880861752  0.1679775595335
            0.4743880861752   0.0512238276497  0.1679775595335
            0.0512238276497   0.4743880861752  0.1679775595335
            0.2385615300181   0.5228769399639  0.2652238803946
            0.5228769399639   0.2385615300181  0.2652238803946
            0.2385615300181   0.2385615300181  0.2652238803946];
    case {'degree=8','degree=9'} % N = 21:
        q = [
            0.0451890097844   0.0451890097844  0.0519871420646
            0.0451890097844   0.9096219804312  0.0519871420646
            0.9096219804312   0.0451890097844  0.0519871420646
            0.7475124727339   0.0304243617288  0.0707034101784
            0.2220631655373   0.0304243617288  0.0707034101784
            0.7475124727339   0.2220631655373  0.0707034101784
            0.2220631655373   0.7475124727339  0.0707034101784
            0.0304243617288   0.7475124727339  0.0707034101784
            0.0304243617288   0.2220631655373  0.0707034101784
            0.1369912012649   0.2182900709714  0.0909390760952
            0.6447187277637   0.2182900709714  0.0909390760952
            0.1369912012649   0.6447187277637  0.0909390760952
            0.2182900709714   0.6447187277637  0.0909390760952
            0.2182900709714   0.1369912012649  0.0909390760952
            0.6447187277637   0.1369912012649  0.0909390760952
            0.0369603304334   0.4815198347833  0.1032344051380
            0.4815198347833   0.0369603304334  0.1032344051380
            0.4815198347833   0.4815198347833  0.1032344051380
            0.4036039798179   0.1927920403641  0.1881601469167
            0.4036039798179   0.4036039798179  0.1881601469167
            0.1927920403641   0.4036039798179  0.1881601469167];
    case {'degree=10','degree=11'} % N = 28:
        q = [...
            0.0000000000000   0.9451704450174  0.0114082494033
            0.9451704450173   0.0000000000000  0.0114082494033
            0.9289002405719   0.0685505797224  0.0132691285720
            0.0685505797224   0.9289002405717  0.0132691285720
            0.0243268355615   0.0243268355616  0.0155865773350
            0.1279662835335   0.0277838749488  0.0408274780428
            0.0277838749488   0.1279662835337  0.0408274780429
            0.0287083428360   0.7498347588657  0.0579849665116
            0.7498347588656   0.0287083428360  0.0579849665116
            0.7228007909707   0.2497602062385  0.0601385247663
            0.2497602062386   0.7228007909707  0.0601385247663
            0.0865562992839   0.8325513856997  0.0625273888433
            0.8325513856998   0.0865562992839  0.0625273888433
            0.3061619157672   0.0303526617491  0.0639684321504
            0.0303526617491   0.3061619157675  0.0639684321504
            0.4868610595047   0.4868610595047  0.0661325872161
            0.6657904293017   0.1765456154219  0.0668503236820
            0.1765456154221   0.6657904293016  0.0668503236821
            0.0293121007360   0.5295657488669  0.0686904305977
            0.5295657488667   0.0293121007360  0.0686904305977
            0.1444673824391   0.1444673824391  0.1002717543859
            0.3299740111411   0.5361815729050  0.1143136784099
            0.5361815729052   0.3299740111409  0.1143136784099
            0.5511507516862   0.1437790861923  0.1223648146752
            0.1437790861923   0.5511507516862  0.1223648146752
            0.3348066587327   0.1529619437161  0.1394422334178
            0.1529619437161   0.3348066587327  0.1394422334178
            0.3430183498147   0.3430183498147  0.1744377829182];
    case {'degree=12','degree=13'} % N = 36:
        q = [...
            0.0242935351590   0.9493059293846  0.0166240998757
            0.0265193427722   0.0242695130640  0.0166811699778
            0.9492126023551   0.0265067966437  0.0166830569067
            0.0033775763749   0.4767316412363  0.0175680870083
            0.4757672298101   0.5198921829102  0.0184474661845
            0.5190783193471   0.0055912706202  0.0197942410188
            0.8616839745321   0.0133996048618  0.0203540395855
            0.1249209759926   0.8613054321334  0.0206852863940
            0.0138565453861   0.1247733717358  0.0208271366086
            0.0211887064222   0.8438438351223  0.0317819778279
            0.8432296787219   0.1354563645830  0.0320472035241
            0.1354231797865   0.0213482820656  0.0320607681146
            0.3088853510679   0.0221919663014  0.0430765959183
            0.6685057595169   0.3089012879389  0.0438473415339
            0.0226545012557   0.6691709943321  0.0439209672733
            0.2808515408772   0.6924718155106  0.0479951923691
            0.6922446749051   0.0268723345026  0.0483806260733
            0.0268617447119   0.2810093973222  0.0484867423375
            0.1141778485470   0.7973581413586  0.0556964488024
            0.7974807922061   0.0879806508791  0.0561026364356
            0.0892807293894   0.1145020561128  0.0565190123693
            0.1052487892455   0.6686904119922  0.0689289890670
            0.6663022280740   0.2275051631832  0.0717213336089
            0.2307803737547   0.1054572561221  0.0727453920976
            0.1705059157540   0.5174064398658  0.0788807336737
            0.5086593973043   0.3170523855209  0.0810114345512
            0.3141823862281   0.1810706361659  0.0825725299055
            0.4617460817864   0.4678594539804  0.0842044567330
            0.0693087496081   0.4622856042085  0.0843585533305
            0.4651955259268   0.0724357805669  0.0851969868488
            0.2578625857893   0.6131395039177  0.0902845328052
            0.6112627766779   0.1300360834609  0.0914283143485
            0.1305182135934   0.2581713828884  0.0916279065409
            0.4281437991828   0.2362005969817  0.1025573374896
            0.3356995783730   0.4311026308588  0.1033159661413
            0.2305424298836   0.3456013949376  0.1035854367193];
    case 'degree=14' % N = 45:
        q = [...
            0.0000000000000   1.0000000000000  0.0010616711990
            1.0000000000000   0.0000000000000  0.0010616711990
            0.0000000000000   0.0000000000000  0.0010616711990
            0.0573330873026   0.0151382269814  0.0131460236101
            0.0573330873026   0.9275286857160  0.0131460236101
            0.9275286857160   0.0573330873026  0.0131460236101
            0.0151382269814   0.0573330873026  0.0131460236101
            0.9275286857160   0.0151382269814  0.0131460236101
            0.0151382269814   0.9275286857160  0.0131460236101
            0.8159625040711   0.1659719969565  0.0242881926949
            0.8159625040711   0.0180654989724  0.0242881926949
            0.1659719969565   0.8159625040711  0.0242881926949
            0.0180654989724   0.8159625040711  0.0242881926949
            0.1659719969565   0.0180654989724  0.0242881926949
            0.0180654989724   0.1659719969565  0.0242881926949
            0.3165475556378   0.0186886898773  0.0316799866332
            0.6647637544849   0.0186886898773  0.0316799866332
            0.0186886898773   0.6647637544849  0.0316799866332
            0.0186886898773   0.3165475556378  0.0316799866332
            0.3165475556378   0.6647637544849  0.0316799866332
            0.6647637544849   0.3165475556378  0.0316799866332
            0.0192662192492   0.4903668903754  0.0349317947036
            0.4903668903754   0.0192662192492  0.0349317947036
            0.4903668903754   0.4903668903754  0.0349317947036
            0.0875134669581   0.8249730660837  0.0383664533945
            0.0875134669581   0.0875134669581  0.0383664533945
            0.8249730660837   0.0875134669581  0.0383664533945
            0.0935526036219   0.2079865423167  0.0578369491210
            0.0935526036219   0.6984608540613  0.0578369491210
            0.2079865423167   0.0935526036219  0.0578369491210
            0.6984608540613   0.0935526036219  0.0578369491210
            0.6984608540613   0.2079865423167  0.0578369491210
            0.2079865423167   0.6984608540613  0.0578369491210
            0.0974892983467   0.5380088595149  0.0725821687394
            0.3645018421383   0.0974892983467  0.0725821687394
            0.5380088595149   0.0974892983467  0.0725821687394
            0.5380088595149   0.3645018421383  0.0725821687394
            0.3645018421383   0.5380088595149  0.0725821687394
            0.0974892983467   0.3645018421383  0.0725821687394
            0.2217145894873   0.5565708210253  0.0897856524107
            0.5565708210253   0.2217145894873  0.0897856524107
            0.2217145894873   0.2217145894873  0.0897856524107
            0.3860471669296   0.2279056661408  0.1034544533617
            0.2279056661408   0.3860471669296  0.1034544533617
            0.3860471669296   0.3860471669296  0.1034544533617];
    case {'degree=15','degree=16'} % N = 55:
        q = [...
            1.0000000000000   0.0000000000000  0.0006202599851
            0.0000000000000   1.0000000000000  0.0006315174712
            0.0000000000000   0.0000000000000  0.0007086601559
            0.9398863583577   0.0049848744634  0.0055163716168
            0.0543806683058   0.9386405618617  0.0062692407656
            0.0093940049164   0.0526424462697  0.0078531408826
            0.0164345086362   0.9469035517351  0.0094551483864
            0.9469487269862   0.0363373677167  0.0097824511271
            0.0426604005768   0.0151224541799  0.0099861643489
            0.0122269495439   0.8693773510664  0.0137553818816
            0.8673696521047   0.1204917285774  0.0140979178040
            0.8456744021389   0.0157763967870  0.0149646864337
            0.1395759632103   0.8448120870375  0.0156097503612
            0.1317821743231   0.0135009605584  0.0157683693348
            0.0157955126300   0.1455274938536  0.0175794546383
            0.7365462884436   0.0155697540908  0.0204113840270
            0.0139688430330   0.7379836894450  0.0209562878616
            0.2547895186039   0.7297615689771  0.0210713412998
            0.7316386522555   0.2543076683315  0.0217646760202
            0.0157253728951   0.2696239795791  0.0222288408699
            0.2662302843647   0.0144783956308  0.0224186693682
            0.8673504065214   0.0591679410400  0.0230122616993
            0.0741493666957   0.8634782575061  0.0236813902500
            0.0159285948360   0.4191238955238  0.0257464643368
            0.0156061028068   0.5809222921146  0.0257956801608
            0.5910094817484   0.0159251452651  0.0258072327610
            0.4034771496889   0.5806700368104  0.0260343232059
            0.5694745628526   0.4149495146302  0.0265768141609
            0.0678493700650   0.0761218678591  0.0265784761831
            0.4265968590272   0.0157509692312  0.0267532329238
            0.0670982507890   0.7741898312421  0.0375787806641
            0.7528310231480   0.0819119495639  0.0383065894195
            0.7753727783557   0.1577128457292  0.0384849695025
            0.1689073157787   0.7503943099742  0.0389619825852
            0.1687335832919   0.0708311507268  0.0394604111547
            0.0821244708436   0.1762996626771  0.0412364778098
            0.6288705363345   0.0807744953317  0.0512872438483
            0.0811413015266   0.3054373589776  0.0516405641935
            0.2969112065080   0.6227485988871  0.0518230042269
            0.0767542314171   0.6247247149546  0.0528527988181
            0.6223022333845   0.3011485821166  0.0538505573027
            0.3103786288051   0.0779098365079  0.0541895329319
            0.0819218215187   0.4603633038351  0.0584737146444
            0.4717022665013   0.0821554006797  0.0592863168363
            0.4546603415250   0.4637565033890  0.0594358276749
            0.1701091339237   0.6422277808188  0.0631800255863
            0.6406004329487   0.1898293537256  0.0632926845153
            0.1912267583717   0.1739955685343  0.0640707361772
            0.1885315767070   0.4798914070406  0.0812040595918
            0.4772929957691   0.3348356598119  0.0814437513530
            0.3126974621760   0.4957972197259  0.0814679201241
            0.4961225945946   0.1927553668904  0.0815050548084
            0.1928805312867   0.3161015807261  0.0815164664939
            0.3360041453816   0.1894892801290  0.0816931059623
            0.3337280550848   0.3343571021811  0.0923218334531];
    case {'degree=17','degree=18'} % N = 66:
        q = [...
            0.0116731059668   0.9812565951289  0.0025165756986
            0.9810030858388   0.0071462504863  0.0025273452007
            0.0106966317092   0.0115153933376  0.0033269295333
            0.9382476983551   0.0495570591341  0.0081503492125
            0.0126627518417   0.9370123620615  0.0086135525742
            0.0598109409984   0.0121364578922  0.0087786746179
            0.0137363297927   0.0612783625597  0.0097099585562
            0.9229527959405   0.0141128270602  0.0102466211915
            0.0633107354993   0.9220197291727  0.0108397688341
            0.0117265100335   0.1500520475229  0.0129385390176
            0.1554720587323   0.8325147121589  0.0136339823583
            0.8343293888982   0.0125228158759  0.0138477328147
            0.8501638031957   0.1371997508736  0.0139421540105
            0.0128816350522   0.8477627063479  0.0144121399968
            0.1510801608959   0.0136526924039  0.0153703455534
            0.0101917879217   0.5770438618345  0.0162489802253
            0.2813372399303   0.7066853759623  0.0169718304280
            0.7124374628501   0.0124569780990  0.0170088532421
            0.2763025250863   0.0121741311386  0.0170953520675
            0.0109658368561   0.4194306712466  0.0173888854559
            0.4289110517884   0.5599616067469  0.0174543962439
            0.4215420555115   0.0116475994785  0.0178406757287
            0.5711258590444   0.0118218313989  0.0178446863879
            0.5826868270511   0.4057889581177  0.0179046337552
            0.0130567806713   0.2725023750868  0.0181259756201
            0.0130760400964   0.7224712523233  0.0184784838882
            0.7263437062407   0.2602984019251  0.0185793564371
            0.0687230068637   0.0631417277210  0.0203217151777
            0.8652302101529   0.0720611837338  0.0213771661809
            0.0648599071037   0.8590433543910  0.0231916854098
            0.1483494943362   0.7888788352240  0.0274426710859
            0.0624359898396   0.1493935499354  0.0290301922340
            0.7871369011735   0.0656382042757  0.0294522738505
            0.0519104921610   0.5255635695605  0.0299436251629
            0.1543129927444   0.0716383926917  0.0307026948119
            0.2617842745603   0.0621479485288  0.0325263365863
            0.7667257872813   0.1658211554831  0.0327884208506
            0.2582103676627   0.6800119766139  0.0331234675192
            0.0679065925147   0.7571515437782  0.0346167526875
            0.5293578274804   0.4121503841107  0.0347081373976
            0.0666036150484   0.2612513087886  0.0347372049404
            0.0585675461899   0.3902236114535  0.0348528762454
            0.0644535360411   0.6373626559761  0.0348601561186
            0.6748138429151   0.0637583342061  0.0355471569975
            0.3914602310369   0.5503238090563  0.0360182996383
            0.6487701492307   0.2836728360263  0.0362926285843
            0.3946498220408   0.0605175522554  0.0381897702083
            0.5390137151933   0.0611990176936  0.0392252800118
            0.1627895082785   0.6861322141035  0.0482710125888
            0.6812436322641   0.1567968345899  0.0489912121566
            0.1542832878020   0.1667512624020  0.0497220833872
            0.2522727750445   0.2504803933395  0.0507065736986
            0.2547981532407   0.4994090649043  0.0509771994043
            0.1485580549194   0.5756023096087  0.0521360063667
            0.2930239606436   0.5656897354162  0.0523460874925
            0.2808991272310   0.1437921574248  0.0524440683552
            0.4820989592971   0.2518557535865  0.0527459644823
            0.5641878245444   0.1462966743153  0.0529449063728
            0.1307699644344   0.4489577586117  0.0542395594501
            0.1479692221948   0.3001174386829  0.0543470203419
            0.5638684222946   0.2813772089298  0.0547100548639
            0.4361157428790   0.4252053446420  0.0557288345913
            0.3603263935285   0.2599190004889  0.0577734264233
            0.4224188334674   0.1453238443303  0.0585393781623
            0.3719001833052   0.3780122703567  0.0609039250680
            0.2413645006928   0.3847563284940  0.0637273964449];
    case {'degree=19','degree=20'} % N = 78:
        q = [...
            0.0089411337112   0.0086983293702  0.0021744545399
            0.9792622629807   0.0102644133744  0.0028987135265
            0.0105475382112   0.9785514202515  0.0030846029337
            0.0023777061947   0.0636551098604  0.0034401633104
            0.0630425115795   0.0041506347509  0.0041898472012
            0.9308422496730   0.0048053482263  0.0044738051498
            0.0629076555490   0.9316790069481  0.0047054420814
            0.9315962246381   0.0626264881801  0.0048867935750
            0.0061951689415   0.9293587058564  0.0051927643369
            0.0287125819237   0.0310202122997  0.0074073058981
            0.9293844478305   0.0342152968219  0.0079755410301
            0.0375457566621   0.9257868884669  0.0083550522910
            0.0086895739064   0.1584971251510  0.0096166660864
            0.1547597053965   0.8363606657688  0.0096318257850
            0.8331025294185   0.0089257244824  0.0098577460758
            0.8374231073526   0.1529167304078  0.0102657880301
            0.1559362505234   0.0094966240058  0.0103188103111
            0.0098599642095   0.8342211493596  0.0106291001630
            0.4055873733289   0.0074389302008  0.0106881306895
            0.5964727898618   0.3956330809311  0.0106969021010
            0.0080747800416   0.4031319425903  0.0109026461714
            0.0075073977721   0.5851609594681  0.0109899783575
            0.3936764519237   0.5974896592899  0.0113423055229
            0.5846530726212   0.0087250464968  0.0120535642930
            0.4870804112120   0.0202129229912  0.0139619193821
            0.2683512811785   0.7202340088668  0.0141147991536
            0.7223956288748   0.2662399366456  0.0141930347046
            0.2716826742357   0.0112882698808  0.0144212676268
            0.0112580842046   0.7169695963325  0.0144704346855
            0.0115034734370   0.2740067110166  0.0144949769872
            0.7140525900564   0.0113511560497  0.0145386775694
            0.4902871053112   0.4936491841468  0.0145964190926
            0.0201423425209   0.4832573459601  0.0147314578466
            0.0361107464859   0.0935679501582  0.0167463963304
            0.8607998819851   0.0397379067075  0.0168955500458
            0.1005891526001   0.8586343419352  0.0169422662884
            0.0918740717058   0.0395513001973  0.0173070172095
            0.8604888296191   0.0966224057079  0.0174524546493
            0.0439842178673   0.8561886349107  0.0177217222159
            0.2011017606735   0.7449115835626  0.0282824024023
            0.7449993726263   0.0536865638166  0.0284996712488
            0.0532186641310   0.1963754275935  0.0285005646539
            0.7453984647401   0.1982065805550  0.0300647223478
            0.1957289932876   0.0555713833156  0.0302031277082
            0.1092532057988   0.6100036182413  0.0303987136077
            0.0567625702001   0.7409121894959  0.0305668796074
            0.0483837933475   0.6075135660978  0.0306067413002
            0.1080612809760   0.1122081510437  0.0309330068201
            0.6185605900991   0.2698753703035  0.0309773820835
            0.7721296013497   0.1114117395333  0.0313146250545
            0.6115734801133   0.3389367677931  0.0313573493392
            0.3381326103376   0.0494693938787  0.0314320469287
            0.1173084128254   0.7696451309795  0.0315182143894
            0.2674551260596   0.1115718808154  0.0324248137985
            0.6542100160026   0.1906548314700  0.0347512152386
            0.0538297481158   0.3358616826849  0.0350393454927
            0.1848840324117   0.1551831523851  0.0350717420310
            0.3376267104744   0.6081402596294  0.0352129215334
            0.6067102034499   0.0542632795598  0.0352615504981
            0.4612614085496   0.0688176670722  0.0366403220343
            0.1525465365671   0.6510240845749  0.0367733107670
            0.0700582543543   0.4661904392742  0.0371675662937
            0.4704201379032   0.4634826455353  0.0373371571606
            0.1216461693746   0.2381494875516  0.0403973346588
            0.6371404052702   0.1238399384513  0.0413580040638
            0.2379904515119   0.6370216452326  0.0421957791870
            0.1483929857177   0.4894188577780  0.0495451004037
            0.3598069571550   0.1452880866253  0.0500419261141
            0.4941441055095   0.3610216383818  0.0505794587115
            0.1440630687981   0.3513508341887  0.0520037210188
            0.5019764440004   0.1435491663293  0.0521533567886
            0.3555423834298   0.5016491599502  0.0524899152358
            0.2443439540771   0.2406052129104  0.0599159762516
            0.2437064989342   0.5109017277055  0.0599609997426
            0.5122200807321   0.2452737973543  0.0599915272129
            0.2526038315178   0.3700319555094  0.0634133183449
            0.3759895652851   0.2505406611631  0.0635311861108
            0.3729077987144   0.3753750277549  0.0637206605672];
    case 'degree=21' % N = 91:
        q = [...
            0.0035524391922   0.0035524391922  0.0006704436439
            0.0035524391922   0.9928951216156  0.0006704436439
            0.9928951216156   0.0035524391922  0.0006704436439
            0.9553548273730   0.0087898929093  0.0045472608074
            0.0358552797177   0.0087898929093  0.0045472608074
            0.9553548273730   0.0358552797177  0.0045472608074
            0.0087898929093   0.0358552797177  0.0045472608074
            0.0087898929093   0.9553548273730  0.0045472608074
            0.0358552797177   0.9553548273730  0.0045472608074
            0.8865264879047   0.1082329745017  0.0052077585320
            0.8865264879047   0.0052405375935  0.0052077585320
            0.0052405375935   0.1082329745017  0.0052077585320
            0.0052405375935   0.8865264879047  0.0052077585320
            0.1082329745017   0.8865264879047  0.0052077585320
            0.1082329745017   0.0052405375935  0.0052077585320
            0.0466397432150   0.9067205135700  0.0065435432887
            0.0466397432150   0.0466397432150  0.0065435432887
            0.9067205135700   0.0466397432150  0.0065435432887
            0.2075720456946   0.0082759241284  0.0092737841533
            0.2075720456946   0.7841520301770  0.0092737841533
            0.7841520301770   0.2075720456946  0.0092737841533
            0.0082759241284   0.7841520301770  0.0092737841533
            0.0082759241284   0.2075720456946  0.0092737841533
            0.7841520301770   0.0082759241284  0.0092737841533
            0.0858119489725   0.0314836947701  0.0095937782623
            0.8827043562574   0.0314836947701  0.0095937782623
            0.0314836947701   0.0858119489725  0.0095937782623
            0.0858119489725   0.8827043562574  0.0095937782623
            0.8827043562574   0.0858119489725  0.0095937782623
            0.0314836947701   0.8827043562574  0.0095937782623
            0.6688778233826   0.0095150760625  0.0114247809167
            0.0095150760625   0.3216071005550  0.0114247809167
            0.0095150760625   0.6688778233826  0.0114247809167
            0.6688778233826   0.3216071005550  0.0114247809167
            0.3216071005550   0.6688778233826  0.0114247809167
            0.3216071005550   0.0095150760625  0.0114247809167
            0.4379999543113   0.0099859785681  0.0117216964174
            0.0099859785681   0.5520140671206  0.0117216964174
            0.4379999543113   0.5520140671206  0.0117216964174
            0.0099859785681   0.4379999543113  0.0117216964174
            0.5520140671206   0.4379999543113  0.0117216964174
            0.5520140671206   0.0099859785681  0.0117216964174
            0.7974931072148   0.0405093994119  0.0188197155232
            0.0405093994119   0.1619974933734  0.0188197155232
            0.0405093994119   0.7974931072148  0.0188197155232
            0.1619974933734   0.7974931072148  0.0188197155232
            0.7974931072148   0.1619974933734  0.0188197155232
            0.1619974933734   0.0405093994119  0.0188197155232
            0.3864215551955   0.3864215551955  0.0235260980271
            0.3864215551955   0.2271568896090  0.0235260980271
            0.2271568896090   0.3864215551955  0.0235260980271
            0.8090129379329   0.0954935310336  0.0235571466151
            0.0954935310336   0.8090129379329  0.0235571466151
            0.0954935310336   0.0954935310336  0.0235571466151
            0.2745425238718   0.0479840480721  0.0268246207430
            0.0479840480721   0.6774734280561  0.0268246207430
            0.6774734280561   0.0479840480721  0.0268246207430
            0.6774734280561   0.2745425238718  0.0268246207430
            0.2745425238718   0.6774734280561  0.0268246207430
            0.0479840480721   0.2745425238718  0.0268246207430
            0.4053472446667   0.5429849622344  0.0314289776779
            0.0516677930989   0.4053472446667  0.0314289776779
            0.4053472446667   0.0516677930989  0.0314289776779
            0.5429849622344   0.0516677930989  0.0314289776779
            0.0516677930989   0.5429849622344  0.0314289776779
            0.5429849622344   0.4053472446667  0.0314289776779
            0.1877738615539   0.1068148267588  0.0337196192159
            0.7054113116872   0.1877738615539  0.0337196192159
            0.7054113116872   0.1068148267588  0.0337196192159
            0.1068148267588   0.7054113116872  0.0337196192159
            0.1877738615539   0.7054113116872  0.0337196192159
            0.1068148267588   0.1877738615539  0.0337196192159
            0.1195059712009   0.3057122990643  0.0427745294213
            0.1195059712009   0.5747817297348  0.0427745294213
            0.5747817297348   0.1195059712009  0.0427745294213
            0.5747817297348   0.3057122990643  0.0427745294213
            0.3057122990643   0.5747817297348  0.0427745294213
            0.3057122990643   0.1195059712009  0.0427745294213
            0.5981245743363   0.2009377128319  0.0441138932737
            0.2009377128319   0.5981245743363  0.0441138932737
            0.2009377128319   0.2009377128319  0.0441138932737
            0.2160775200005   0.3121360256673  0.0461469594684
            0.3121360256673   0.2160775200005  0.0461469594684
            0.2160775200005   0.4717864543321  0.0461469594684
            0.3121360256673   0.4717864543321  0.0461469594684
            0.4717864543321   0.3121360256673  0.0461469594684
            0.4717864543321   0.2160775200005  0.0461469594684
            0.4376579903849   0.4376579903849  0.0469152468624
            0.4376579903849   0.1246840192303  0.0469152468624
            0.1246840192303   0.4376579903849  0.0469152468624
            0.3333333333333   0.3333333333333  0.0551199980347];
    case {'degree=22','degree=23'} % N = 105:
        q = [...
            0.0087809303836   0.9903676436772  0.0006438298261
            0.9903675314220   0.0087809216232  0.0006438413076
            0.0027029276450   0.0335914404439  0.0010134735710
            0.0335909214524   0.0027028946710  0.0010134752576
            0.0091675068606   0.0091676353051  0.0019679929935
            0.9675568182558   0.0084737176656  0.0033467313784
            0.0084737200688   0.9675569435345  0.0033467339208
            0.0078781948792   0.0676784943862  0.0042873323375
            0.0676785477700   0.0078781659291  0.0042873459885
            0.9470266955047   0.0442974541187  0.0043003801372
            0.0442974755680   0.9470266676487  0.0043003849098
            0.9144243214882   0.0081735455132  0.0056934629205
            0.0081735424459   0.9144244234031  0.0056934640134
            0.2497452292741   0.3833232434720  0.0061643868015
            0.3833232646055   0.2497451268005  0.0061644756418
            0.8876850353557   0.1035328809446  0.0062014513591
            0.1035329228297   0.8876849931840  0.0062014531952
            0.0077255923618   0.1403190991974  0.0069636330294
            0.1403192425107   0.0077255934624  0.0069636331842
            0.8104591009652   0.1809642523926  0.0075066257720
            0.1809643003717   0.8104590515334  0.0075066264565
            0.8330767948684   0.0083010939677  0.0079074768339
            0.0083010907126   0.8330768545392  0.0079074772485
            0.0348407706147   0.0348406969482  0.0080353344623
            0.2740287679608   0.7173981847948  0.0087963441074
            0.7173982224778   0.2740287304386  0.0087963448112
            0.2394976858234   0.0081859182262  0.0091304195716
            0.0081859185845   0.2394975566677  0.0091304213611
            0.0068836152075   0.4843740892687  0.0092821748751
            0.4843741485699   0.0068836232949  0.0092821815662
            0.4960767772741   0.4960767529507  0.0094499806178
            0.6112936776245   0.3804323691239  0.0094627468484
            0.3804323980345   0.6112936466533  0.0094627485294
            0.7303890713524   0.0083987179701  0.0095555772285
            0.0083987168639   0.7303890895407  0.0095555792843
            0.6128525675612   0.0075475979695  0.0096138842488
            0.0075475961037   0.6128525484582  0.0096138846826
            0.0079525316513   0.3559773826721  0.0099991524212
            0.3559774870460   0.0079525358502  0.0099991551850
            0.9110236977966   0.0437233665345  0.0100301319277
            0.0437233605166   0.9110236807446  0.0100301346636
            0.0388480061835   0.0967030908282  0.0124936676185
            0.0967032117936   0.0388479942386  0.0124936726125
            0.0873226911312   0.0873226620391  0.0140197309137
            0.0421445202084   0.8485617789108  0.0143336216896
            0.8485617974961   0.0421445420915  0.0143336272125
            0.8477921333864   0.1067435942472  0.0153604142740
            0.1067435889398   0.8477921328146  0.0153604183425
            0.1833966521991   0.0416340521608  0.0184523825614
            0.0416340541167   0.1833965196930  0.0184523863146
            0.7611632251560   0.1941599202852  0.0195833983573
            0.1941599254144   0.7611632153938  0.0195834019994
            0.7579378747173   0.0439826608586  0.0197632751342
            0.0439826512395   0.7579378242308  0.0197632766677
            0.0369760535918   0.5363186076436  0.0198806391019
            0.5363187134342   0.0369760780935  0.0198806485776
            0.1001256948921   0.7912267093545  0.0207181838484
            0.7912266693524   0.1001257554673  0.0207181934893
            0.0379866714177   0.4157413128558  0.0208943071440
            0.4157414028965   0.0379867061535  0.0208943251956
            0.6507106491463   0.0420141226713  0.0214864573885
            0.0420141133438   0.6507105645084  0.0214864586007
            0.0425548444254   0.2920626023484  0.0222218133036
            0.2920627107240   0.0425548546753  0.0222218160203
            0.5389729538180   0.4193031469005  0.0223345305455
            0.4193031828489   0.5389729093610  0.0223345378739
            0.6549472009700   0.3007352636162  0.0224758924946
            0.3007352790917   0.6549471812731  0.0224758980440
            0.3752400771585   0.3453980130752  0.0229701395845
            0.3453980282786   0.3752400695673  0.0229703394438
            0.0994532168761   0.1598308695187  0.0232798376102
            0.1598309359585   0.0994531960132  0.0232798427506
            0.1797326661667   0.7124585430924  0.0269483199647
            0.7124584461943   0.1797327722240  0.0269483307107
            0.1066065678636   0.7001701784175  0.0280438758010
            0.7001701904096   0.1066065855677  0.0280438764607
            0.0993303629801   0.6065647984796  0.0287526270172
            0.6065648052521   0.0993303896769  0.0287526387271
            0.1023223542704   0.2533381579528  0.0298980829063
            0.2533382324938   0.1023223826189  0.0298980922759
            0.6166226715217   0.2769502060575  0.0309004358516
            0.2769500693109   0.6166227900624  0.0309004385956
            0.0904184571873   0.4981522637001  0.0314031017088
            0.4981522767248   0.0904185045149  0.0314031073955
            0.0928231860168   0.3738418516908  0.0319191553024
            0.3738418699229   0.0928232584790  0.0319191668378
            0.2521678840407   0.2521680925697  0.0321429924062
            0.5087500218708   0.3905580544330  0.0330395601388
            0.3905579116731   0.5087501437661  0.0330395631829
            0.1706141469096   0.5266738039554  0.0356169095589
            0.5266737761312   0.1706142257537  0.0356169276054
            0.3487581527629   0.2588055084886  0.0365741189998
            0.2588053596017   0.3487583491703  0.0365741515204
            0.1696614558053   0.3013522183964  0.0365977646990
            0.3013521806875   0.1696615963219  0.0365978053889
            0.2580202409759   0.4584741774478  0.0369945680114
            0.4584740860198   0.2580203819011  0.0369945775059
            0.1848898683498   0.1848898704551  0.0374053623787
            0.6130740338465   0.1921611994069  0.0375550258317
            0.1921611750994   0.6130740398389  0.0375550312530
            0.4180541160599   0.1650613336416  0.0388887693486
            0.1650612642036   0.4180541199244  0.0388887708342
            0.5159205739625   0.2982719005229  0.0392705643548
            0.2982718935750   0.5159205534362  0.0392705802517
            0.4098894602340   0.4098894317792  0.0398766879831];
    case {'degree=24','degree=25'} % N = 120:
        q = [...
            0.0082881595033   0.9848202768869  0.0014873417859
            0.4618422030241   0.5381577969759  0.0014889035262
            0.0071066441239   0.0080842361390  0.0015005944380
            0.9847613141699   0.0070015755134  0.0015059208313
            0.5374447869049   0.4625552130951  0.0015318868715
            0.0000000000000   0.4887676880140  0.0023032634487
            0.4914131929361   0.0000000000000  0.0023649067042
            0.0070345937020   0.9574158053697  0.0028751143611
            0.9564734714228   0.0364655449485  0.0029862488735
            0.0370198792045   0.0070908577166  0.0030384162737
            0.1024124542747   0.8936125594937  0.0032092459688
            0.5928065811509   0.0049451705600  0.0037029598435
            0.0050948422371   0.0996676659189  0.0037407186035
            0.0081562023689   0.0415561148784  0.0038452543223
            0.0424936107568   0.9494865260352  0.0038670778668
            0.9495543500844   0.0081794507292  0.0039192555178
            0.8932787471239   0.0053224326262  0.0039573282688
            0.0069317612927   0.9065401020433  0.0044032251724
            0.9035839030665   0.0894771171077  0.0045907108173
            0.0905665738209   0.0070525342005  0.0047023669435
            0.0083929332787   0.6663179931111  0.0050014843818
            0.6261245686071   0.0092197583153  0.0052387830156
            0.0062801592979   0.8335207460527  0.0054422104092
            0.8272539257367   0.1665134939330  0.0056931248912
            0.0062005875353   0.7424693255229  0.0059107422989
            0.1676900311185   0.0065717743528  0.0059687967687
            0.7199353069567   0.0064354534962  0.0067262190287
            0.2749740090237   0.7185296120719  0.0068307848624
            0.0079257582005   0.1766411374714  0.0069531259112
            0.0069981220752   0.2704767254004  0.0072460270642
            0.8125248773263   0.0082299533210  0.0072728189613
            0.0073536969970   0.5934167875453  0.0073008930847
            0.7283665935411   0.2648817553752  0.0073604666776
            0.1800642304565   0.8115848976682  0.0074119923255
            0.2658102467762   0.0068553525429  0.0074892214336
            0.0070892364520   0.3757632659744  0.0078604067260
            0.3774054302043   0.6148573533757  0.0078621726423
            0.0369649608668   0.9210792302893  0.0080506361066
            0.9203194109805   0.0426025082114  0.0081442860473
            0.0425477806431   0.0372689941794  0.0081478804152
            0.6191278394983   0.3724055713809  0.0092444146612
            0.3762697209178   0.0081436422011  0.0094674635165
            0.0956111149690   0.8771098372601  0.0097132210137
            0.0302473410377   0.0943858903393  0.0099753581151
            0.8739905691754   0.0313198990883  0.0103367803673
            0.8604133734958   0.1049019782046  0.0112263277166
            0.0347307852352   0.8609856462886  0.0114309118745
            0.1043606608343   0.0357152881004  0.0115550567487
            0.7797622824754   0.1872318199265  0.0135575856957
            0.0185865164256   0.4834397678794  0.0135984962900
            0.0324585286618   0.7783474916042  0.0137754813837
            0.8371293901157   0.0804060570156  0.0137961015942
            0.0836602075315   0.8421414817051  0.0138408839904
            0.0784070242501   0.0849927089145  0.0140634019977
            0.4929238648458   0.4892855914710  0.0140991451009
            0.1870637584073   0.0345210858281  0.0142004111991
            0.4892636967025   0.0190774755077  0.0144518424517
            0.0401982618372   0.1691143187109  0.0150245979639
            0.7894259278865   0.0412206731484  0.0152817804122
            0.1686260456429   0.7894860640585  0.0155550724169
            0.3750901913174   0.5895318272013  0.0164570886000
            0.0356362876880   0.3681256217699  0.0165275759573
            0.5887548164804   0.0359968962541  0.0166847554451
            0.0373308082182   0.6790704673533  0.0167409312985
            0.2820769993374   0.0373639992361  0.0168674663361
            0.6819277603320   0.2803330345725  0.0168882230165
            0.0374938324382   0.2634016180014  0.0172087112691
            0.6984079204127   0.0364154673322  0.0174681068264
            0.2654390894079   0.6980717436193  0.0176663899614
            0.1429848440800   0.7612254618453  0.0182967621475
            0.7623554007647   0.0943741220275  0.0183576852459
            0.0934222022749   0.1479799836832  0.0186392569521
            0.5759004479923   0.3821329641698  0.0189781060590
            0.3822427332525   0.0426716362301  0.0191847922578
            0.0411414081675   0.5718082874432  0.0194080442044
            0.0802462538379   0.7702204382042  0.0194720072193
            0.7625229819410   0.1559420577362  0.0200855080495
            0.1524941445131   0.0842965421322  0.0201673909332
            0.0622159195833   0.4538181318873  0.0221742162761
            0.1109539036076   0.4586014071171  0.0229702440508
            0.4575627212057   0.4795313560210  0.0233465117399
            0.4322865136374   0.1230591237472  0.0234883135338
            0.5865002850241   0.0834119779793  0.0240682099018
            0.0869359250818   0.6755677013351  0.0240910792953
            0.0929594906936   0.2326500892727  0.0245677049481
            0.6661932141454   0.2448294007406  0.0246536315719
            0.4780306362227   0.0661749044835  0.0246756530052
            0.4372215294577   0.4442145585244  0.0249704602710
            0.6779224504669   0.0929096534577  0.0250026544082
            0.2423431255660   0.0889793655129  0.0250490869426
            0.2288925420305   0.6780053081672  0.0250936250125
            0.3315065049959   0.5847381559741  0.0251482076226
            0.3424200526607   0.5139245722736  0.0255010290447
            0.0862630046475   0.3340976249234  0.0256544511979
            0.5113188946635   0.1380154720554  0.0257974750630
            0.1538977841001   0.6788062619562  0.0270007753993
            0.6779951348472   0.1663358925269  0.0274431536844
            0.1664600469411   0.1582214504849  0.0277072401488
            0.0950910318888   0.5666590332543  0.0278284415364
            0.3436048136712   0.0978960873457  0.0287207381105
            0.5560417025366   0.3468917820947  0.0288826834956
            0.1452404029513   0.3599534491052  0.0293302729759
            0.1619685156238   0.5810131373330  0.0318902879557
            0.5800164844262   0.2560674640672  0.0319083660286
            0.2450201223288   0.5881469552102  0.0320938960329
            0.2557621891794   0.1652244065047  0.0321618608780
            0.2205239985511   0.3496507466106  0.0322424127534
            0.4940183111285   0.2549448448453  0.0327072446421
            0.2531570689798   0.2543369115017  0.0329946316695
            0.5846891116357   0.1666603916479  0.0331828096025
            0.1660333602278   0.2523240191705  0.0334857162651
            0.2505426292461   0.4959007627528  0.0335468472792
            0.3519336802182   0.1805380367800  0.0337049042988
            0.3502668835419   0.4358582329881  0.0340361462767
            0.4400892485512   0.2120576104941  0.0342465235323
            0.4680855471546   0.3552681570774  0.0345528817251
            0.1770237763947   0.4670352922266  0.0356782875703
            0.3900920779501   0.3323152819300  0.0364656225016
            0.2805847774120   0.3898041176680  0.0365172708706
            0.3361523347440   0.2778500044356  0.0371924811018];
end

xq = q(:,1);
yq = q(:,2);
wq = q(:,3)/4;
NQ = length(xq);

end

function gdof=getgdof(k,ip,mesh,type)
%
% vertici
%
bgdof = mesh.polygon(ip).vertices;
%
% k-1 per lato
%
if k>1
    for ie=[2 3 1]
        %
        % numero globale dell'edge
        %
        ieg = mesh.polygon(ip).edges(ie);
        %
        % aggiorno i gradi di libertà globali
        %
        gdofie = (((ieg-1)*(k-1)+1):((ieg-1)*(k-1)+(k-1)));
        %
        % devo stabilire l'orientamento dell'edge
        %
        if mesh.polygon(ip).vertices(ie) == mesh.edge(ieg).v1
            % stesso orientamento
            bgdof = [bgdof mesh.NV+gdofie(1:end)];
        elseif mesh.polygon(ip).vertices(ie) == mesh.edge(ieg).v2
            % orientamento contrario
            bgdof = [bgdof mesh.NV+gdofie(end:-1:1)];
        else
            fprintf('%s: wrong edge',mfilename)
        end
    end
    %
    % Npkm3 per triangolo
    %
    % dimensione dello spazio di polinomi di grado k-3 in due variabili
    %
    Npkm3 = (k-2)*(k-1)/2;
    %
    % dof interni
    %
    igdof = [mesh.NV+mesh.NE*(k-1)+(((ip-1)*Npkm3+1):((ip-1)*Npkm3+Npkm3))];
    %
else
    igdof = [];
end
%
if nargin == 3
    %
    % defauly: all dofs
    %
    gdof = [bgdof igdof];
    %
elseif nargin == 4
    %
    switch type
        case 'boundary'
            gdof = bgdof;
        case 'internal'
            gdof = igdof;
        otherwise
            fprintf('%s: bad type argument: %s',mfilename,type)
            return
    end
else
    fprintf('%s: bad number of arguments',mfilename)
end

end


function VEM=splitdof(VEMOLD,dof)

VEM = VEMOLD;

p = VEM.p;

Np = dimP(p);
Npm2 = dimP(p-2);

VEM.dof = dof;
%
% decompose VEM.dof
%
if 1<=VEM.dof && VEM.dof<=VEM.NV
    %
    % vertex dof
    %
    VEM.iv = VEM.dof;
    %
    VEM.ie = 0;
    VEM.ip = 0;
    %
    VEM.ialpha = 0;
    %
elseif VEM.NV+1<=VEM.dof && VEM.dof<=VEM.NV+(p-1)*VEM.NV
    %
    % edge dof
    %
    VEM.iv = 0;
    %
    VEM.ie = floor((VEM.dof-VEM.NV-1)/(p-1))+1;
    VEM.ip = VEM.dof - VEM.NV - (VEM.ie-1)*(p-1);
    %
    VEM.ialpha = 0;
    %
    % check
    %
    VEM.mydof = VEM.NV + (VEM.ie-1)*(p-1) + VEM.ip;
    %
    if not(VEM.mydof == VEM.dof)
        error('bad edge dof calculation');
    end
    %
elseif p*VEM.NV+1<=VEM.dof && VEM.dof<=p*VEM.NV+Np
    %
    % internal dof
    %
    VEM.iv = 0;
    %
    VEM.ie = 0;
    VEM.ip = 0;
    %
    VEM.ialpha = VEM.dof - p*VEM.NV;
    %
    if VEM.dof>p*VEM.NV+Npm2
        warning('using enhanced VEM!')
    end
else
    error('bad dof')
end

end

function z = VEM_m(alpha,x,y,VEM)

% polinomi che definiscono i gradi di libertà
%
% alpha = multiindex (alpha(1),alpha(2))

z = ((x-VEM.xb)/VEM.h).^alpha(1) .* ((y-VEM.yb)/VEM.h).^alpha(2);

end

function z=g(VEM,marker,x,y)

p = VEM.p;
ip = VEM.ip;

NV = VEM.NV;

iv = VEM.iv;
ie = VEM.ie;
ialpha = VEM.ialpha;

xv = VEM.xv;
yv = VEM.yv;

if ialpha>0
    %
    % internal dof - nothing to do
    %
    z = 0;
    %
    return
    %
end
%
% is a boundary dof
%
% could be a vertex dof or an edge dof
%
if iv > 0
    %
    % is a vertex dof
    %
    ivp = iv+1;
    if ivp==NV+1
        ivp = 1;
    end
    ivm = iv-1;
    if ivm==0
        ivm = NV;
    end
    %
    if marker == 100*iv
        %
        % right vertex
        %
        z = 1;
        %
    else
        %
        if marker==ivm
            %
            % compute distance from iv
            %
            d = norm([x-xv(iv), y-yv(iv)]);
            %
            % lunghezza edge
            %
            le = norm([xv(iv)-xv(ivm), yv(iv)-yv(ivm)]);
            %
            s = d/le;
            %
            z = polyval(polyfit(VEM.tgl',[1 zeros(1,p)],p),2*s-1);
            %
        elseif marker==iv 
            %
            % compute distance from iv
            %
            d = norm([x-xv(iv), y-yv(iv)]);
            %
            % lunghezza edge
            %
            le = norm([xv(ivp)-xv(iv), yv(ivp)-yv(iv)]);
            %
            s = d/le;
            %
            z = polyval(polyfit(VEM.tgl',[1 zeros(1,p)],p),2*s-1);
            %
        else
            z = 0;
        end
    end
    %
    return
    %
end
%
% can only be an edge dof
%
if ie>0
    %
    % is an edge dof
    %
    if marker == ie
        %
        % we are on the right edge
        %
        iep = ie+1;
        if iep==NV+1
            iep = 1;
        end
        %
        % compute distance from ie
        %
        d = norm([x-xv(ie), y-yv(ie)]);
        %
        % lunghezza edge
        %
        le = norm([xv(iep)-xv(ie), yv(iep)-yv(ie)]);
        %
        s = d/le;
        %        
        z = polyval(polyfit(VEM.tgl',[zeros(1,ip) 1 zeros(1,p-ip)],p),2*s-1);
        %
    else
        %
        % different edge
        %
        z = 0;
        %
    end
end
                
end

function drawuh(uh,k,mesh)

%title(sprintf('p%d solution\n',k))

%Xv = NaN*ones(maxNV+(k-1)*maxNE+1,mesh.NP);
%Yv = NaN*ones(maxNV+(k-1)*maxNE+1,mesh.NP);
%Zv = NaN*ones(maxNV+(k-1)*maxNE+1,mesh.NP);
%Wv = NaN*ones(maxNV+(k-1)*maxNE+1,mesh.NP);

Xp = [];
Yp = [];
Zp = [];

for ip=1:mesh.NP

    X = zeros(3*k,1);
    Y = zeros(3*k,1);
    Z = zeros(3*k,1);
    
    % valori nel vertici
    
    for iv=1:3
        X(k*(iv-1)+1) = mesh.vertex(mesh.polygon(ip).vertices(iv)).x;
        Y(k*(iv-1)+1) = mesh.vertex(mesh.polygon(ip).vertices(iv)).y;
        Z(k*(iv-1)+1) = uh(mesh.polygon(ip).vertices(iv));
    end
    
    for iv=1:3
        
        ivp  = mod(iv,3)+1;
        ivpp = mod(ivp,3)+1;

        x1 = mesh.vertex(mesh.polygon(ip).vertices(ivp)).x;
        y1 = mesh.vertex(mesh.polygon(ip).vertices(ivp)).y;

        x2 = mesh.vertex(mesh.polygon(ip).vertices(ivpp)).x;
        y2 = mesh.vertex(mesh.polygon(ip).vertices(ivpp)).y;
                
        ieg = mesh.polygon(ip).edges(ivp);

        for ik=1:k-1
            
            xik = (1-ik/k)*x1+(ik/k)*x2;
            yik = (1-ik/k)*y1+(ik/k)*y2;
            
            X(k*(ivp-1)+1+ik) = xik;
            Y(k*(ivp-1)+1+ik) = yik;
            
            if mesh.polygon(ip).vertices(ivp) == mesh.edge(ieg).v1
                Z(k*(ivp-1)+1+ik)           = uh(mesh.NV+(ieg-1)*(k-1)+ik);
            else
                Z(k*(ivp-1)+1+((k-1)-ik+1)) = uh(mesh.NV+(ieg-1)*(k-1)+ik);
            end
            
        end
    end
    
    Xp = [Xp X];
    Yp = [Yp Y];
    Zp = [Zp Z];

end
    
patch(Xp,Yp,Zp,'w','LineWidth',0.1,'EdgeALpha',0.5)
%patch(Xp,Yp,Zp,Zp,'EdgeALpha',0.5)
colormap('jet')
%patch(Xp,Yp,Zp,Zp)

Xe = [];
Ye = [];
Ze = [];

for ie=1:mesh.NE
    
    X = zeros(k+1,1);
    Y = zeros(k+1,1);
    Z = zeros(k+1,1);
    
    if mesh.edge(ie).marker>0
        
        v1 = mesh.edge(ie).v1;
        v2 = mesh.edge(ie).v2;
        x1 = mesh.vertex(v1).x;
        y1 = mesh.vertex(v1).y;
        x2 = mesh.vertex(v2).x;
        y2 = mesh.vertex(v2).y;
        
        X(1) = x1;
        Y(1) = y1;
        X(k+1) = x2;
        Y(k+1) = y2;
        Z(1) = uh(v1);
        Z(k+1) = uh (v2);
        
        for ik=1:k-1
            
            xik = (1-ik/k)*x1+(ik/k)*x2;
            yik = (1-ik/k)*y1+(ik/k)*y2;
            
            X(1+ik) = xik;
            Y(1+ik) = yik;
            
            Z(1+ik) = uh(mesh.NV+(ie-1)*(k-1)+ik);
            
        end
        
        Xe = [Xe X];
        Ye = [Ye Y];
        Ze = [Ze Z];
        
    end
    
end

line(Xe,Ye,Ze,'Color','k','LineWidth',1)
    
end