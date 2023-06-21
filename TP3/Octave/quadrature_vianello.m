function [xq,yq,wq]=quadrature_vianello(n,ip,mesh,COMPRESS)
%
% [XQ,YQ,WQ]=QUADRATURE_VIANELLO(N,IP,MESH,COMPRESS)
%
% Returns Vianello quadrature nodes (XQ,YQ) and weights (WQ) which
% integrate exactly polynomials of degree N on polygon IP of mesh MESH
%
% If the argument COMPRESS exists and is FALSE, *does not* performs Vianello
% compression. Otherwise, Vianello compression is always done.
%
if nargin==3
    COMPRESS = true;
end
%
WARNINGON = false;
%
% se il lato è un segmento se ho grado k dentro mi serve grado k+1 sull'edge
%
% 2*Nge-1 = k+1 da cui Nge = ceil(k/2)+1
%
% e poi mi serve grado k dentro
%
% 2*Ngi-1 = k da cui Ngi = ceil((k+1)/2)
%
% quindi
%
Nge = ceil(n/2)+1;
Ngi = ceil((n+1)/2);
%
if WARNINGON
    fprintf('%s: using Vianello with prescribed degree of precision k = %d\n',mfilename,n);
end
%
% nge Gauss points
%
[tge,twge] = mylgwt(Nge,-1,1);
%
% ngi Gauss points
%
[tgi,twgi] = mylgwt(Ngi,-1,1);
%
xge = zeros(size(tge'));
yge = zeros(size(tge'));
%
xgi = zeros(size(tgi'));
ygi = zeros(size(tgi'));
%
% vertex center of polygon
%
xvc = sum([mesh.vertex([mesh.polygon(ip).vertices]).x])/mesh.polygon(ip).NV;
yvc = sum([mesh.vertex([mesh.polygon(ip).vertices]).y])/mesh.polygon(ip).NV;
%
xq = [];
yq = [];
wq = [];
%
for iv=1:mesh.polygon(ip).NV
    %
    % get global number of edge iv
    %
    ieg = mesh.polygon(ip).edges(iv);
    %
    % length
    %
    ledge = mesh.edge(ieg).length;
    %
    % get second local extreme of the edge
    %
    ivp = iv+1;
    if ivp==mesh.polygon(ip).NV+1
        ivp = 1;
    end
    %
    % get global extrema of the edge (locally ordered)
    %
    v1 = mesh.polygon(ip).vertices(iv);
    v2 = mesh.polygon(ip).vertices(ivp);
    %
    % get coordinates of extrema
    %
    x1 = mesh.vertex(v1).x;
    y1 = mesh.vertex(v1).y;
    %
    x2 = mesh.vertex(v2).x;
    y2 = mesh.vertex(v2).y;
    %
    % normal to edge with respect to the polygon
    %
    nxe = +(y2-y1)/ledge;
    nye = -(x2-x1)/ledge;
    %
    % gauss points and weight on edge ieg
    %
    xge = (x2-x1)/2*tge+(x2+x1)/2;
    yge = (y2-y1)/2*tge+(y2+y1)/2;
    wge = ledge/2*twge;
    %
    % Gauss points and weigths on the internal line
    %
    for ige=1:Nge
        %
        xgi = (xge(ige)-xvc)/2*tgi+(xge(ige)+xvc)/2;
        ygi(:) = yge(ige);
        wgi = (xge(ige)-xvc)/2*twgi;
        %
        xq = [xq, xgi'];
        yq = [yq, ygi];
        wq = [wq, wge(ige)*nxe*wgi'];
        
%         [ygi,wgi] = mylgwt(Ngi,yvc,yge(ige));
%         xgi(:) = xge(ige);
%     
%         xq = [xq; xgi];
%         yq = [yq; ygi];
%         wq = [wq; wgi*(wge(ige)*ledge/2)*nyge(ige)];
        
    end
    
end

Nq = length(xq);

if WARNINGON
    fprintf('%s: found %d quadrature nodes\n',mfilename,Nq);
end

if COMPRESS
    %
    Nqold = Nq;
    %
    deg = n;
    X = [xq' yq'];
    omega = wq';
    pos = 0;
    %
    [pts,w,momerr] = comprexcub(deg,X,omega,pos);
    %
    % INPUT:
    % deg: polynomial exactness degree
    % X: 2-column array of point coordinates
    % omega: 1-column array of weights
    % pos: NNLS for pos=1, QR with column pivoting for pos=0
    % OUTPUT:
    % pts: 2-column array of extracted points
    % w: 1-column array of corresponding weights (positive for pos=1)
    % momerr: moment reconstruction error
    % total-degree Chebyshev-Vandermonde matrix at X
    %
    xq = pts(:,1)';
    yq = pts(:,2)';
    wq = w';
    %
    Nq = length(xq);
    %
    perc = Nq/Nqold*100;
    %
    if WARNINGON
        fprintf('%s: compression: quadrature nodes reduced from %d to %d (%f%%)\n',mfilename,Nqold,Nq,perc);
    end
    %
end

end
    


function [pts,w,momerr] = compresscub(deg,X,omega,pos)
% INPUT:
% deg: polynomial exactness degree
% X: 2-column array of point coordinates
% omega: 1-column array of weights
% pos: NNLS for pos=1, QR with column pivoting for pos=0
% OUTPUT:
% pts: 2-column array of extracted points
% w: 1-column array of corresponding weights (positive for pos=1)
% momerr: moment reconstruction error
% total-degree Chebyshev-Vandermonde matrix at X
rect=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];
V=chebvand(deg,X,rect);
% polynomial basis orthogonalization
[Q,R]=qr(V,0);
% tiny imaginary parts may appear
Q=real(Q);
% possible re-orthogonalization
% [Q,R1]=qr(Q,0);
% moments of the orthogonal basis
orthmom=Q'*omega;
% weigths computation
if pos == 1
weights=lsqnonneg(Q',orthmom);
else
weights=Q'\orthmom;
end
% indexes of nonvanishing weights and compression
ind=find(abs(weights)>0);
pts=X(ind,:);
w=weights(ind);
% moment reconstruction error
% bivariate Chebyshev basis
mom=V'*omega;
momerr=norm(V(ind,:)'*w-mom);
% discrete OP basis
% momerr=norm(Q(ind,:)’*w-orthmom);
end

function V = chebvand(deg,X,rect);
% INPUT:
% deg = polynomial degree
% X = 2-column array of point coordinates
% rect = 4-component vector such that the rectangle
% [rect(1),rect(2)] x [rect(3),rect(4)] contains X
% OUTPUT:
% V = Chebyshev-Vandermonde matrix at X
% couples with length less or equal to deg
j=(0:1:deg);
[j1,j2]=meshgrid(j);
good=find(j1(:)+j2(:)<deg+1);
couples=[j1(good) j2(good)];
% mapping X in the square [-1,1]ˆ2
a=rect(1);b=rect(2);c=rect(3);d=rect(4);
map=[(2*X(:,1)-b-a)/(b-a) (2*X(:,2)-d-c)/(d-c)];
% total-degree Chebyshev-Vandermonde matrix at X
V1=cos(couples(:,1)*acos(map(:,1)'));
V2=cos(couples(:,2)*acos(map(:,2)'));
V=(V1.*V2)';
end


function [pts,w,momerr] = comprexcub(deg,X,omega,pos)

% compression of bivariate cubature formulas by Tchakaloff points 
% or approximate Fekete points
% useful, for example, in node reduction of algebraic cubature formulas
% see the web page: http://www.math.unipd.it/~marcov/Tchakaloff.html

% by Federico Piazzon, Alvise Sommariva and Marco Vianello
% University of Padova, May 2016


% INPUT:   
% deg: polynomial exactness degree
% X: 2-column array of point coordinates
% omega: 1-column array of weights
% pos: NNLS for pos=1, QR with column pivoting for pos=0

% OUTPUT:
% pts: 2-column array of extracted points
% w: 1-column array of corresponding weights (positive for pos=1)
% momerr: moment reconstruction error


% FUNCTION BODY 
% total-degree Chebyshev-Vandermonde matrix at X
rect=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];
V=chebvandx(deg,X,rect);
% polynomial basis orthogonalization 
[Q,R]=qr(V,0);
% tiny imaginary parts could appear
Q=real(Q);
% possible re-orthogonalization
% [Q,R]=qr(Q,0);

% moments of the orthogonal basis
orthmom=Q'*omega;
% weigths computation
if pos == 1
% Tchakaloff points (positive weights)
weights=lsqnonneg(Q',orthmom);
else
% approximate Fekete points (possible negative weights) 
weights=Q'\orthmom;
end
% indexes of nonvanishing weights and compression
ind=find(abs(weights)>0);
pts=X(ind,:);
w=weights(ind);

% moment reconstruction error
% bivariate Chebyshev basis
% mom=V'*omega;
% momerr=norm(V(ind,:)'*w-mom);
% discrete OP basis
momerr=norm(Q(ind,:)'*w-orthmom);

end


function V = chebvandx(deg,X,rect)

% INPUT:
% deg = polynomial degree
% X = 2-column array of point coordinates
% rect = 4-component vector such that the rectangle
% [rect(1),rect(2)] x [rect(3),rect(4)] contains X

% OUTPUT:
% V = Chebyshev-Vandermonde matrix at X, graded lexic. order

% FUNCTION BODY 
% rectangle containing the mesh 
if isempty(rect) 
rect=[min(X(:,1)) max(X(:,1)) min(X(:,2)) max(X(:,2))];
end
    
% couples with length less or equal to deg
% graded lexicographical order
j=(0:1:deg);
[j1,j2]=meshgrid(j);
dim=(deg+1)*(deg+2)/2;
couples=zeros(dim,2);
for s=0:deg
good=find(j1(:)+j2(:)==s);
couples(1+s*(s+1)/2:(s+1)*(s+2)/2,:)=[j1(good) j2(good)];
end

% mapping the mesh in the square [-1,1]^2
a=rect(1);b=rect(2);c=rect(3);d=rect(4);
map=[(2*X(:,1)-b-a)/(b-a) (2*X(:,2)-d-c)/(d-c)];

% Chebyshev-Vandermonde matrix on the mesh
T1=chebpolys(deg,map(:,1));
T2=chebpolys(deg,map(:,2));
V=T1(:,couples(:,1)+1).*T2(:,couples(:,2)+1);

end



function T=chebpolys(deg,x)

% computes the Chebyshev-Vandermonde matrix on the real line by recurrence

% INPUT:
% deg = maximum polynomial degree
% x = 1-column array of abscissas

% OUTPUT
% T: Chebyshev-Vandermonde matrix at x, T(i,j+1)=T_j(x_i), j=0,...,deg

T=zeros(length(x),deg+1);
t0=ones(length(x),1);
T(:,1)=t0;
t1=x;
T(:,2)=t1;

for j=2:deg
t2=2*x.*t1-t0;
T(:,j+1)=t2;
t0=t1;
t1=t2;
end

end


