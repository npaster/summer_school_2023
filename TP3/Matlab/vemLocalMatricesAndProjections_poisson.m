%
%%% VEM local Matrices and Projections for the polygon ip
%
function mat=vemLocalMatricesAndProjections_poisson(k,ip,mesh)
%
% choosing Pz for k=1
%
whichPz = 'vertices';
%whichPz = 'perimeter';
%
%%% --------------------------- QUADRATURE -------------------------------
%
xq = mesh.polygon(ip).xq;
yq = mesh.polygon(ip).yq;
wq = mesh.polygon(ip).wq;
%
%%% GAUSS-LOBATTO STUFF
%
% (tgl,wgl) = Gauss-Lobatto quadrature nodes and weigths on [-1 +1]
%
[tgl,wgl] = mylglnodes(k);
% ------------------------------------------------------------------------
%
%%% POLYGONAL QUANTITIES
%
% vertex coordinates; add first vertex at bottom
%
x = [mesh.vertex([mesh.polygon(ip).vertices]).x mesh.vertex(mesh.polygon(ip).vertices(1)).x];
y = [mesh.vertex([mesh.polygon(ip).vertices]).y mesh.vertex(mesh.polygon(ip).vertices(1)).y];
%
% centroid of the polygon
%
xb = mesh.polygon(ip).xb;
yb = mesh.polygon(ip).yb;
%
% area of the polygon
%
area = mesh.polygon(ip).area;
%
% diameter of the polygon
%
h = mesh.polygon(ip).h;
%
% number of edges and vertices (must be equal!)
%
NE = mesh.polygon(ip).NE;
NV = mesh.polygon(ip).NV;
%
% edges of the polygon locally oriented
%
edge = zeros(2,NE);
%
for ie=1:NE
    edge(:,ie) = [x(ie+1)-x(ie) y(ie+1)-y(ie)]';
end
%
% length of the edges
%
ledge = zeros(NE,1);
%
for ie=1:NE
    ledge(ie) = norm(edge(:,ie));
end
%
% normal vectors to the polygon (one every edge)
%
normal = zeros(2,NE);
%
for ie=1:NE
    %
    % rotating edge
    %
    normal(1,ie) = +edge(2,ie);
    normal(2,ie) = -edge(1,ie);
    %
    % normalize
    %
    normal(:,ie) = normal(:,ie)/ledge(ie);
    %
end
%
%%% DEGREES OF FREEDOM
%
NumberVertexDof = NV;
NumberInternalNodeDof = NE*(k-1);
NumberInternalElementDof = dimP(k-2);
Ndof = NumberVertexDof + NumberInternalNodeDof + NumberInternalElementDof;
%
%%% monomial scaled basis
%
[mono, dimPk] = generateMonomialBasis(k,ip,mesh);
%
%%% first derivative of the monomial scaled basis
%
[mono_x, mono_y] = computeFirstDerivativeForMonomialBasis(k,ip,mesh);
%
%%% second derivative of the monomial scaled basis
%
[mono_xx, mono_yy] = computeSecondDerivativeForMonomialBasis(k,ip,mesh);
%%%-----------------------------------------------------------------------
%
%% --------------------------- LOCAL MATRICES ----------------------------
%
%%% ------------------------------- MATRIX G -----------------------------
%
% G = int_K (grad(mono_i)*grad(mono_j))
%
G = zeros(dimPk,dimPk);
%
for ialpha=1:dimPk
    %
    for ibeta=1:ialpha-1
        G(ialpha,ibeta) = G(ibeta,ialpha);
    end
    %
    for ibeta=ialpha:dimPk
        %
         G(ialpha,ibeta) = sum((...
            mono_x(ialpha).p(xq,yq).*mono_x(ibeta).p(xq,yq)+...
            mono_y(ialpha).p(xq,yq).*mono_y(ibeta).p(xq,yq)...
        ).*wq);
        %
    end
end
%
% We need to fix the first row of the matrix G
%
for ialpha=1:dimPk
    %
    G(1,ialpha) = 0;
    %
    if k==1 
        if strcmp(whichPz,'vertices') 
            %
            % Pz m_alpha = mean on the vertices
            %
            for iv=1:NV
                G(1,ialpha) = G(1,ialpha) + mono(ialpha).p(x(iv),y(iv));
            end
            %
            % divide by the number of vertices
            %
            G(1,ialpha) = G(1,ialpha)/NV; 
            %
        elseif strcmp(whichPz,'perimeter') 
            %
            % Pz m_alpha = media on the boundary of the polygon             
            %
            for ie=1:NE
                iep = ie+1;
                G(1,ialpha) = G(1,ialpha) + ...
                    ledge(ie)*(mono(ialpha).p(x(ie),y(ie))+mono(ialpha).p(x(iep),y(iep)))/2;
            end
            %
            % divide by the perimeter of the polygon
            %
            G(1,ialpha) = G(1,ialpha)/sum(ledge);  
            %
        else
            %
            fprintf('%s: Gc: wrong Pz for k=1\n',mfilename);
            return
            %
        end   
    else
        %
        % for k>=2 always Pz m_alpha = mean on the polygon
        %
        G(1,ialpha) = sum(mono(ialpha).p(xq,yq).*wq);
        %
        % divide by the area of the polygon
        %
        G(1,ialpha) = G(1,ialpha)/area;  
        %
    end
end
mat.G=G;
%
%%% ------------------------------- MATRIX B -----------------------------
%
B = zeros(dimPk,Ndof);
%% first row
if k==1
    for idof=1:NV
        %
        if strcmp(whichPz,'vertices')
            %
            % Pz = mean on the vertices
            %
            B(1,idof) = 1/NV;
            %
        elseif strcmp(whichPz,'perimeter') 
            %
            % Pz = media on the boundary of the polygon   
            %
            im = i-1;
            if im==0
                im = NV;
            end
            %
            B(1,idof) = 0.5*(ledge(im)+ledge(i))/sum(ledge);
            %
        else
            %
            fprintf('%s: Gc: wrong Pz for k=1',mfilename);
            return
        end
    end
else
    %
    % k>=2
    %
    % mean value on the polygon
    % it is always zero except at the first internal dof
    %
    B(1,NumberVertexDof+NumberInternalNodeDof+1) = 1;
    %
end
%
% for the second row 
%
%
for ialpha=2:dimPk
    idof=1;
    %
    % vertex dofs get contribution from two edges
    %
    for iv=1:NV
        %
        % phi_iv at vertex iv
        %
        %        iv         ivm
        %        *-----------*
        %       /     ivm
        %      /iv
        %     /
        %    *
        %   ivp
        %
        ivm = iv-1;
        ivp = iv+1;
        if ivm == 0
            ivm = NE;
        end
        %
        B(ialpha,idof) = ...
            (mono_x(ialpha).p(x(iv),y(iv))*normal(1,iv)+...
            mono_y(ialpha).p(x(iv),y(iv))*normal(2,iv))*...
            ledge(iv)/2*wgl(1);
        %
        % we integrate on edge ivm, from vertex iv to vertex ivm
        %
        % only the first one survives!
        %
        B(ialpha,idof) = B(ialpha,idof) + ...
            (mono_x(ialpha).p(x(iv),y(iv))*normal(1,ivm)+...
            mono_y(ialpha).p(x(iv),y(iv))*normal(2,ivm))*...
            ledge(ivm)/2*wgl(1); 
        idof=idof+1;
    end
    if k>=2
        %
        % edge dofs exist only of k=2 and get contribution from one edge
        %
        for ie=1:NE
            %
            iep = ie+1;
            %
            %  *--------------*
            % iep     ie      ie
            %
            for ik=1:k-1
                %
                % phi_ik with ik = internal node on edge ie
                %
                % this is the only node surviving
                %
                xs = (x(iep)-x(ie))/2*tgl(ik+1)+(x(iep)+x(ie))/2;
                ys = (y(iep)-y(ie))/2*tgl(ik+1)+(y(iep)+y(ie))/2;
                ws = wgl(ik+1)*ledge(ie)/2;
                %
                %
                % we integrate on edge ie, from vertex ie to vertex iep
                %
                %
                B(ialpha,idof) = ...
                    (mono_x(ialpha).p(xs,ys)*normal(1,ie)+...
                    mono_y(ialpha).p(xs,ys)*normal(2,ie))*...
                    ws;
                %   
                idof=idof+1;
                %
            end
            %
        end
        %
        % internal dofs exist only if k>=2 
        %
        % phi_{k*NV+k*(NE-1)+ialpha} with ialpha=1:dimP(k-2)
        %
        for ibeta=1:dimP(k-2)
            %
            multiIndex= o2m(ibeta);
            %
            if ((multiIndex(1)==mono_xx(ialpha).indexX) && (multiIndex(2)==mono_xx(ialpha).indexY))
                    B(ialpha,idof) = B(ialpha,idof) -area * mono_xx(ialpha).coeff;
            end
            if ((multiIndex(1)==mono_yy(ialpha).indexX) && (multiIndex(2)==mono_yy(ialpha).indexY))
                    B(ialpha,idof) = B(ialpha,idof) -area * mono_yy(ialpha).coeff;
            end
            %
            idof=idof+1;
        end
    end
    %
end 
mat.B=B;
%
%%% ----------- MATRIX D of the degrees of freedom -----------------------
%
D = zeros(Ndof,dimPk);
%
for ialpha=1:dimPk
    %
    idof=1;
    %
    % vertex dof
    %
    for iv=1:NV
        %
        D(idof,ialpha) = mono(ialpha).p(x(iv),y(iv));
        idof=idof+1;
        %
    end
    %
    % next dofs exist only if k>=2
    %
    if k>=2
        %
        % internal Gauss-Lobatto points dofs
        %
        for ie=1:NE
            %
            % ie edge, ie-iep are the local extrema
            %
            iep = ie+1;
            %
            % value of m_alpha at the internal Gauss-Lobatto points
            %
            for ik=1:k-1
                %
                % ik Gauss-Lobatto point on the edge
                %
                xik = (x(iep)-x(ie))/2*tgl(ik+1)+(x(iep)+x(ie))/2;
                yik = (y(iep)-y(ie))/2*tgl(ik+1)+(y(iep)+y(ie))/2;
                %
                % m_alpha evaluated at (xik,yik)
                %
                D(idof,ialpha) = mono(ialpha).p(xik,yik);
                %
                idof=idof+1;
                %
            end
                %
                %
        end
        %
        % internal dofs
        %
        for ibeta=1:dimP(k-2)
            %
            % order |beta| moment of m_alpha
            %
            D(idof,ialpha) = 1/area * sum(...
                mono(ialpha).p(xq,yq).*...
                mono(ibeta).p(xq,yq)...
                .*wq);
            %
            idof=idof+1;
            %
        end
    end
end
mat.D=D;
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% check that G=B*D
%
errGminusBD = max(max(abs(G-B*D)));
%
if errGminusBD>1e-12
    fprintf('%s: G=BD check - error is %e\n',mfilename,errGminusBD)
end
%
% Projection matrices
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% matrix PI NABLA *
%
mat.PNkstar = G\B;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% matrix PI nabla
%
mat.PNk = D*mat.PNkstar;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% matrix G tilde
%
mat.Gt = G;
mat.Gt(1,:) = 0;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX H
%
H = zeros(dimPk,dimPk);
%
for ialpha=1:dimPk
    
    for ibeta=1:ialpha-1
        H(ialpha,ibeta) = H(ibeta,ialpha);
    end
    
    for ibeta=ialpha:dimPk    
        H(ialpha,ibeta) = sum(mono(ialpha).p(xq,yq).*mono(ibeta).p(xq,yq).*wq);
    end
end
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% matrix C
%
C = zeros(dimPk,Ndof);
%
% Identity block
%
for ialpha=1:dimP(k-2)
    C(ialpha,NumberVertexDof+NumberInternalNodeDof+ialpha) = area;
end
%
HPNkstar = H * mat.PNkstar;
%
for ialpha=dimP(k-2)+1:dimPk
    for i=1:Ndof
        C(ialpha,i) = HPNkstar(ialpha,i);
    end
end
%
mat.C = C;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% matrix P0kstar
%
mat.P0kstar = H\C;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% matrix P0k
%
mat.P0k = D*mat.P0kstar;
return
