function mesh=square2mesh(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7)
%
% square2mesh(arg0)
%      tries to load ./meshes/arg0.mat
%   
% square2mesh(arg0,arg1)
%      domain = unit square
%      switch arg0
%          case {'s','t','square','triangle'}
%                structured mesh, arg1=xparts and yparts
%          case 'u'
%                unstructured triangular mesh, arg1=area
%          case 'v'
%                voronoi mesh, arg1=polygons, 20 lloyd iterations
%      end
%
% square2mesh(arg0,arg1,arg2)
%      domain = unit square
%      switch arg0
%          case {'s','t','square','triangle'}
%                structured mesh, arg1=xparts, arg2=yparts
%                with 't' squares are cut left-to-right
%                with 't-' square are cut right-to-left
%          case 'u'
%                NA
%          case 'v'
%                voronoi mesh, arg1=polygons, arg2=Lloyd iteraitions
%      end
%
% square2mesh(arg0,arg1,arg2,arg3)
%      domain = unit square
%      switch arg0
%          case 'd-s'
%                structured deformed quad mesh, arg1=xparts, arg2=yparts, arg3=deformation
%          case 'dd-s'
%                structured deformed2 quad mesh, arg1=xparts, arg2=yparts, arg3=deformation
%      end
%
% square2mesh(arg0,arg1,arg2,arg3,arg4,arg5)
%      domain = rectangle (arg1,arg2)-(arg3,arg4)
%      switch type
%          case {'s','t','square','triangle'}
%                structured mesh, arg5=xparts and yparts
%          case 'u'
%                unstructured triangular mesh, arg5=area
%          case 'v'
%                voronoi mesh, arg5=polygons, 20 Lloyd iterations
%      end
%
%
% square2mesh(arg0,arg1,arg2,arg3,arg4,arg5,arg6)
%      domain = rectangle (arg1,arg2)-(arg3,arg4)
%      switch type
%          case {'s','t','square','triangle'}
%                structured mesh, arg5=xparts, arg6=yparts
%          case 'u'
%                unstructured mesh, arg5=area, arg6=dummy
%          case 'v'
%                voronoi mesh, arg5=polygons, arg6=Lloyd iterations
%      end
%
% square2mesh(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7)
%      domain = rectangle (arg1,arg2)-(arg3,arg4)
%      switch type
%          case {'d-s'}
%                structured deformed quad mesh, arg5=xparts, arg6=yparts, arg7=deformation
%          case {'dd-s'}
%                structured deformed2 quad mesh, arg5=xparts, arg6=yparts, arg7=deformation
%      end
%
defaultLloyd = 20;
%
VERBOSE = false;
%
% reset random generator to replicate the result
%
%stream =RandStream.getGlobalStream;
%reset(stream);

%
switch nargin
    
    case 1
        
        % load meshes/arg0.mat
        
        meshname = arg0;
        
        meshfilename = ['./meshes/' meshname];
        fprintf('%s: info: loading mesh %s\n',mfilename,meshfilename)
        load([meshfilename '.mat'],'mesh');
        
    otherwise
        
        type = arg0;
        
        switch type
            case {'square','s'}
                type = 's';
            case {'deformed-square','d-s'}
                type = 'd-s';
            case {'dd-s'}
                type = 'dd-s';
            case {'triangle','t','t+'}
                type = 't';
            case 't-'
                type = 't-';
            case 't*'
                type = 't*';
            case {'babuaziz','b'}
                type = 'b';
            case {'unstructured','u'}
                type = 'u';
            case {'voronoi','v'}
                type = 'v';
            otherwise
                fprintf('%s: wrong type %s\n\n',mfilename,type);
                return
        end
        
        switch nargin
            case 2
                x1 = 0.0;
                y1 = 0.0;
                x2 = 1.0;
                y2 = 1.0;
                switch type
                    case {'t','t+','t-','t*','s','b'}
                        Nx = arg1;
                        Ny = Nx;
                    case {'u'}
                        h = arg1;
                    case {'v'}
                        NP = arg1;
                        niter = defaultLloyd;
                end
            case 3
                x1 = 0.0;
                y1 = 0.0;
                x2 = 1.0;
                y2 = 1.0;
                switch type
                    case {'t','t+','t-','t*','s','b'}
                        Nx = arg1;
                        Ny = arg2;
                    case {'u'}
                        h = arg1;
                    case {'v'}
                        NP = arg1;
                        niter = arg2;
                end
            case 4
                x1 = 0.0;
                y1 = 0.0;
                x2 = 1.0;
                y2 = 1.0;
                %
                Nx = arg1;
                Ny = arg2;
                d = arg3;
            case 6
                x1 = arg1;
                y1 = arg2;
                x2 = arg3;
                y2 = arg4;
                switch type
                    case {'t','t+','t-','t*','s','b'}
                        Nx = arg5;
                        Ny = Nx;
                    case {'u'}
                        h = arg5;
                    case {'v'}
                        NP = arg5;
                        niter = defaultLloyd;
                end
            case 7
                x1 = arg1;
                y1 = arg2;
                x2 = arg3;
                y2 = arg4;
                switch type
                    case {'t','t+','t-','t*','s','b'}
                        Nx = arg5;
                        Ny = arg6;
                    case {'u'}
                        h = arg5;
                    case {'v'}
                        NP = arg5;
                        niter = arg6;
                end
            case 8
                x1 = arg1;
                y1 = arg2;
                x2 = arg3;
                y2 = arg4;
                %
                Nx = arg5;
                Ny = arg6;
                d = arg7;
            otherwise
                fprintf('%s: # arguments can only be 2, 3, 4, 6, 7, or 8 - not %d\n\n',mfilename,nargin);
                mesh = [];
                return
        end
        
        switch type
            case {'t','t+'}
                [xn,yn,nodi] = uniformtriangles(x1,y1,x2,y2,Nx,Ny,+1);
            case 't-'
                [xn,yn,nodi] = uniformtriangles(x1,y1,x2,y2,Nx,Ny,-1);
            case 't*'
                [xn,yn,nodi] = uniformtriangles(x1,y1,x2,y2,Nx,Ny,0);
            case 'b'
                [xn,yn,nodi] = uniformtriangles_babuaziz(x1,y1,x2,y2,Nx,Ny);
            case 's'
                [xn,yn,nodi] = uniformsquares(x1,y1,x2,y2,Nx,Ny);
            case 'd-s'
                [xn,yn,nodi] = deformedsquares(x1,y1,x2,y2,Nx,Ny,d);
            case 'dd-s'
                [xn,yn,nodi] = deformedsquares2(x1,y1,x2,y2,Nx,Ny,d);
            case 'u'
                [xn,yn,nodi] = unstructuredtriangles(x1,y1,x2,y2,h);
            case 'v'
                [xn,yn,nodi] = voronoi(x1,y1,x2,y2,NP,niter);
        end
        
        mesh = conn2mesh(xn,yn,nodi);
        
        switch type
            case 'u'
                mesh.name = ['square2mesh-area=' num2str(h)];
            case {'t','t+','t-','t*','s'}
                mesh.name = ['square2mesh-' num2str(Nx) 'x' num2str(Ny)];
            case {'d-s'}
                mesh.name = ['square2mesh-deformed-' num2str(Nx) 'x' num2str(Ny)];
            case {'dd-s'}
                mesh.name = ['square2mesh-deformed2-' num2str(Nx) 'x' num2str(Ny)];
            case {'v'}
                mesh.name = ['square2mesh-polygons=' num2str(NP)];
        end
        
end

end

function [xn,yn,nodi] = uniformsquares(x1,y1,x2,y2,Nx,Ny)
%
[xn,yn,nodi] = deformedsquares(x1,y1,x2,y2,Nx,Ny,0);
%
end

function [xn,yn,nodi] = deformedsquares(x1,y1,x2,y2,Nx,Ny,d)
%
hx = (x2-x1)/Nx;
hy = (y2-y1)/Ny;
%
nnod = (Nx+1)*(Ny+1);
%
nele = Nx*Ny;
%
xn = zeros(nnod,1);
yn = zeros(nnod,1);
%
% nodi
%
nodi = zeros(nele,4);
%
for i=1:Nx+1
    if mod(i,2)==1
        si = -1;
    else
        si = +1;
    end
    %
    for j=1:Ny+1
        if mod(j,2)==1
            sj = +1;
        else
            sj = -1;
        end
        %
        inod = (j-1)*(Nx+1)+i;
        if 1<i && i<Nx+1
            xn(inod) = x1 + (i-1)*hx + si*sj*d*hx;
        else
            xn(inod) = x1 + (i-1)*hx;
        end
        %
        yn(inod) = y1 + (j-1)*hy;
    end
end
%
% elementi
%
for k=1:Nx
    for l=1:Ny
        iele = (l-1)*Nx + k;
        % nodi corrispondenti all'elemento iele = (k,l)
        nodi(iele,1) = (l-1)*(Nx+1)+k;
        nodi(iele,2) = (l-1)*(Nx+1)+k+1;
        nodi(iele,3) = l*(Nx+1)+k+1;
        nodi(iele,4) = l*(Nx+1)+k;
    end
end
%
end

function [xn,yn,nodi] = deformedsquares2(x1,y1,x2,y2,Nx,Ny,d)
%
hx = (x2-x1)/Nx;
hy = (y2-y1)/Ny;
%
nnod = (Nx+1)*(Ny+1);
%
nele = Nx*Ny;
%
xn = zeros(nnod,1);
yn = zeros(nnod,1);
%
% nodi
%
nodi = zeros(nele,4);
%
for i=1:Nx+1
    %
    for j=1:Ny+1
        %
        inod = (j-1)*(Nx+1)+i;
        %
        if i==1 || i==Nx+1
            xn(inod) = x1 + (i-1)*hx;
        else
            xn(inod) = x1 + (i-1)*hx + d*2*(rand-0.5)*hx;
        end    
        %
        if j==1 || j==Ny+1
            yn(inod) = y1 + (j-1)*hy;
        else
            yn(inod) = y1 + (j-1)*hy + d*2*(rand-0.5)*hy;
        end
        %
    end
end
%
% elementi
%
for k=1:Nx
    for l=1:Ny
        iele = (l-1)*Nx + k;
        % nodi corrispondenti all'elemento iele = (k,l)
        nodi(iele,1) = (l-1)*(Nx+1)+k;
        nodi(iele,2) = (l-1)*(Nx+1)+k+1;
        nodi(iele,3) = l*(Nx+1)+k+1;
        nodi(iele,4) = l*(Nx+1)+k;
    end
end
%
end

function [xn,yn,nodi]=uniformtriangles(x1,y1,x2,y2,Nx,Ny,direction)

nnod = (Nx+1)*(Ny+1);
nele = 2*Nx*Ny;

xn = zeros(nnod,1);
yn = zeros(nnod,1);
nodi = zeros(nele,3);

ug=inline('j*(Nx+1)+i+1','i','j','Nx');

for i=0:Nx
    for j=0:Ny
        inod = ug(i,j,Nx);
        %
        xn(inod) = x1 + i/Nx * (x2-x1);
        yn(inod) = y1 + j/Ny * (y2-y1);
    end
end

for i=1:Nx
    for j=1:Ny
        q = (j-1)*Nx+i;
        iele1 = 2*q-1;
        iele2 = 2*q;
        
        xb = (xn(ug(i-1,j-1,Nx))+xn(ug(i,j-1,Nx))+xn(ug(i,j,Nx))+xn(ug(i-1,j,Nx)))/4;
        yb = (yn(ug(i-1,j-1,Nx))+yn(ug(i,j-1,Nx))+yn(ug(i,j,Nx))+yn(ug(i-1,j,Nx)))/4;
        
        nodi(iele1,1) = ug(i-1,j-1,Nx);
        nodi(iele1,2) = ug( i ,j-1,Nx);
        nodi(iele2,1) = ug(i-1, j ,Nx);
        nodi(iele2,3) = ug( i , j ,Nx);
        
        switch direction
            case 1
                nodi(iele1,3) = ug(i-1, j ,Nx);
                nodi(iele2,2) = ug( i ,j-1,Nx);
            case -1
                nodi(iele1,3) = ug( i , j ,Nx);
                nodi(iele2,2) = ug(i-1,j-1,Nx);
            case 0
                if (xb>(x1+x2)/2 && yb>(y1+y2)/2) || (xb<(x1+x2)/2 && yb<(y1+y2)/2)
                    nodi(iele1,3) = ug( i , j ,Nx);
                    nodi(iele2,2) = ug(i-1,j-1,Nx);
                else
                    nodi(iele1,3) = ug(i-1, j ,Nx);
                    nodi(iele2,2) = ug( i ,j-1,Nx);
                end
        end
    end
end

end

function [xn,yn,nodi]=uniformtriangles_babuaziz(x1,y1,x2,y2,Nx,Ny)

nnod = (Nx+1)*(Ny+1);
nele = 2*Nx*Ny;

xn = zeros(nnod,1);
yn = zeros(nnod,1);
nodi = zeros(nele,3);

ug=inline('j*(Nx+1)+i+1','i','j','Nx');

for i=0:Nx
    for j=0:Ny
        inod = ug(i,j,Nx);
        %
        xn(inod) = x1 + i/Nx * (x2-x1);
        yn(inod) = y1 + j/Ny * (y2-y1);
    end
end

for i=1:Nx
    for j=1:Ny
        q = (j-1)*Nx+i;
        iele1 = 2*q-1;
        iele2 = 2*q;
        
        if mod(i,2)~=mod(j,2)
            
            nodi(iele1,1) = ug(i-1,j-1,Nx);
            nodi(iele1,2) = ug( i ,j-1,Nx);
            nodi(iele1,3) = ug( i , j ,Nx);
            
            nodi(iele2,1) = ug(i-1,j-1 ,Nx);
            nodi(iele2,2) = ug( i , j ,Nx);
            nodi(iele2,3) = ug(i-1 , j ,Nx);
            
        else
            
            nodi(iele1,1) = ug(i-1,j-1,Nx);
            nodi(iele1,2) = ug( i ,j-1,Nx);
            nodi(iele1,3) = ug(i-1, j ,Nx);
            
            nodi(iele2,1) = ug(i-1, j ,Nx);
            nodi(iele2,2) = ug( i ,j-1,Nx);
            nodi(iele2,3) = ug( i , j ,Nx);
            
        end
    end
end

end

function [xn,yn,nodi] = unstructuredtriangles(x1,y1,x2,y2,h)

omega = './triangle/tmp/omega';

area = h^2*sqrt(3)/4;

% create poly file with our data

fid = fopen([omega '.poly'],'w');

fprintf(fid,'# numero vertici, 2, 0, marker imposto (1 si, 0 no)\n');
fprintf(fid,'4 2 0 1\n');
fprintf(fid,'# vertice, x, y, marker\n');
fprintf(fid,'1 %e %e +1\n',x1,y1);
fprintf(fid,'2 %e %e +1\n',x2,y1);
fprintf(fid,'3 %e %e +1\n',x2,y2);
fprintf(fid,'4 %e %e +1\n',x1,y2);
fprintf(fid,'# number of segments, marker imposto (1 si, 0 no)\n');
fprintf(fid,'4 1\n');
fprintf(fid,'# lato del bordo, estremo 1, estremo 2, marker\n');
fprintf(fid,'1 1 2 +1\n');
fprintf(fid,'2 2 3 +1\n');
fprintf(fid,'3 3 4 +1\n');
fprintf(fid,'4 4 1 +1\n');
fprintf(fid,'# numero di buchi\n');
fprintf(fid,'0\n');
fprintf(fid,'# numero di vincoli di area\n');
fprintf(fid,'1\n');
fprintf(fid,'# regione, punto nella regione, vincolo di area\n');
fprintf(fid,'1 %e %e %e\n',(x1+x2)/2,(y1+y2)/2,area);

fclose(fid);

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
        system(['.\triangle\PCWIN64\triangle.exe -pqIa ' omega]);
        %
    case 'GLNXA64'
        %
        % Linux on x86_64
        %
        system(['./triangle/GLNXA64/triangle -pqIa ' omega]);
        %
    case 'MACI64'
        %
        % Apple Mac OS X on x86_64
        %
        system(['./triangle/MACI64/triangle -pqIa ' omega]);
        %
end

% legge la struttura dati dentro matlab

% leggiamo i nodi

fid = fopen(strcat(omega,'.node'));

% leggiamo la prima riga

tmp = fscanf(fid,'%d',4); % tmp = vettore

% numero dei nodi

nnod = tmp(1);

% creiamo i vettori che conterranno le coord dei nodi

xn = zeros(nnod,1); % vettore colonna di nnod zeri
yn = zeros(nnod,1);

% creiamo un vettore per il flag fuori-dentro

nodemarker = zeros(nnod,1);

% leggiamo le coord e i flag dei nodi

for inod=1:nnod
    
    tmp = fscanf(fid,'%f',4);
    
    xn(inod) = tmp(2);         % coord x
    yn(inod) = tmp(3);         % coord y
    
    nodemarker(inod) = tmp(4); % flag
end

fclose(fid); % chiudo il file dei nodi

% adesso leggo i triangoli

fid = fopen(strcat(omega,'.ele'));

% prima riga

tmp = fscanf(fid,'%d',3);

% numero di triangoli

nele = tmp(1);

% prepariamo la matrice triangoli -- nodi

nodi = zeros(nele,3);

for iele=1:nele
    tmp = fscanf(fid,'%d',4);
    nodi(iele,1) = tmp(2);
    nodi(iele,2) = tmp(3);
    nodi(iele,3) = tmp(4);
end

fclose(fid);

end

function [xn,yn,nodi]=voronoi(x1,y1,x2,y2,NP,niter)

VERBOSE = false;

dx = x2-x1;
dy = y2-y1;

rx = rand(1,NP);
ry = rand(1,NP);

%rx1 = rand(1,100*NP);
%ry1 = rand(1,100*NP);

%rrx = [];
%rry = [];

%for i=1:length(rx)
%    if rx(i)<0.9
%        rrx = [rrx rx(i)];
%        rry = [rry ry(i)];
%    end
%end

%for i=1:length(rx1)
%    if rx1(i)>0.9
%        rrx = [rrx rx1(i)];
%        rry = [rry ry1(i)];
%    end
%end

%rx = rrx;
%ry = rry;

for remove=[]
    %for remove=[37 56 62]
    rx(remove) = [];
    ry(remove) = [];
end

NP = length(rx);

xs = dx * rx;
ys = dy * ry;

for iter=0:niter
    %
    % specchio i punti
    %
    xss = [xs xs dx+(dx-xs) dx+(dx-xs) dx+(dx-xs) xs -xs -xs -xs];
    yss = [ys -ys -ys ys dy+(dy-ys) dy+(dy-ys) dy+(dy-ys) ys -ys];
    %
    % costruisco voronoi
    %
    if VERBOSE
        fprintf('%s: calling voronoi\n\n',mfilename)
    end
    dt = delaunayTriangulation(xss',yss');
    [V,C] = voronoiDiagram(dt);
    %
    % disegno voronoi
    %
    %if NP<MAXNPDRAWING
    %    figure(1)
    %    voronoi(dt)
    %    ontop
    %    atleft
    %end
    %
    for ip=1:NP
        vv = C{ip};
        xv = V(vv,1);
        yv = V(vv,2);
        %
        xv = [xv; xv(1)];
        yv = [yv; yv(1)];
        %
        NV = length(vv);
        %
        % area of polygon
        %
        tmp = zeros(NV,1);
        %
        for iv=1:NV
            tmp(iv) = xv(iv)*yv(iv+1)-xv(iv+1)*yv(iv);
        end
        %
        area = 0.5*sum(tmp);
        %
        % calcolo baricentro
        %
        xb = 0;
        yb = 0;
        %
        for iv=1:NV
            xb = xb + (xv(iv)+xv(iv+1))*tmp(iv);
            yb = yb + (yv(iv)+yv(iv+1))*tmp(iv);
        end
        %
        xb = xb/(6*area);
        yb = yb/(6*area);
        %
        xs(ip) = xb;
        ys(ip) = yb;
    end
end

% coordinate dei vertici
xt = V(:,1);
yt = V(:,2);
Nt = length(xt);

% calcolo i nodi effettivi - quelli dei poligoni da 1 a NP
%
if VERBOSE
    fprintf('%s: computing effective vertices\n\n',mfilename)
end
NE = [];
%figure(1)
for ip=1:NP
    %    if mod(ip,1000)==1
    %        fprintf('%s: polygon %d out of %d\n',mfilename,ip,NP)
    %    end
    %    patch(xt(C{ip}),yt(C{ip}),zeros(size(V(C{ip},2))),'w')
    NE = union(NE,C{ip});
end
NE = NE';

% elimino i nodi non effettivi
%
xt(setdiff(1:Nt,NE)) = NaN;
yt(setdiff(1:Nt,NE)) = NaN;

% aggiusto i nodi di bordo
%
TOL = 1e-14;
%
xt(abs(xt)<TOL) = 0;
yt(abs(yt)<TOL) = 0;
%
xt(abs(xt-dx)<TOL) = dx;
yt(abs(yt-dy)<TOL) = dy;

% definisco i nodi complessi
%
zt = xt+sqrt(-1)*yt;

if VERBOSE
    fprintf('%s: remove duplicate vertices\n\n',mfilename)
end
for inod=NE
    if xt(inod) == 0 || xt(inod) == dx || yt(inod) == 0 || yt(inod) == dy
        equals = find(abs(zt(inod)-zt)<TOL);
        if length(equals)>1
            keep = equals(1);
            discardv = equals(2:end)';
            for ip=1:NP
                Cip = C{ip};
                for discard=discardv
                    Cip(Cip==discard) = keep;
                end
                C{ip} = unique(Cip,'stable');
            end
            xt(discardv) = NaN;
            yt(discardv) = NaN;
            zt(discardv) = NaN;
        end
    end
end

% trovo i nodi finali
%
NF = [];
for ip=1:NP
    NF = union(NF,C{ip});
end
NF = NF';

if VERBOSE
    fprintf('%s: computing final vertices\n',mfilename)
end

INVNF = zeros(1,max(NF));
for i=1:length(NF)
    INVNF(NF(i)) = i;
end
for ip=1:NP
    if mod(ip,101)==1
        if VERBOSE
            fprintf('%s: polygon %d out of %d\n',mfilename,ip,NP)
        end
    end
    C{ip} = INVNF(C{ip});
end

nnod = length(NF);
nele = NP;

NVMAX = 0;

for ip=1:NP
    NV = length(C{ip});
    if NV>NVMAX
        NVMAX = NV;
    end
end

nodi = zeros(NP,NVMAX);

for ip=1:NP
    nodi(ip,1:length(C{ip})) = C{ip};
end

xn = xt(NF)+x1;
yn = yt(NF)+y1;

end

