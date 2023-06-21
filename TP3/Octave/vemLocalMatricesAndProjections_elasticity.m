%
%%% VEM local Matrices and Projections for the polygon ip for elasticity
%%% problems
%
function mat=vemLocalMatricesAndProjections_elasticity(k,ip,mesh)
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
%%% ELASTICITY INFORMATION -> lamè coefficients
%
lambda=mesh.ELASTICITYVEM.lambda;
mu=mesh.ELASTICITYVEM.mu;
%
%%% local DEGREES OF FREEDOM (outline)
% 
%  vertexDof - first component
%  vertexDof - second component
%  internalEdgeDof - first component
%  internalEdgeDof - second component
%  internalElementDof - first component
%  internalElementDof - second component
% -----------------------------------------
NumberOfLocalDofs = 2*(NV+(k-1)*NE+dimP(k-2));
%
%% --------------------------- LOCAL MATRICES ----------------------------
% Basis of Polynomial
[monoVect, dimPkVect] = generateMonomialBasisVect(k,ip,mesh);
%
grad_eps = computeSymmetricGradientForMonomialBasisVect(k,ip,mesh);
%
[symmetricMatBasis, dimMatBasis] = generateSymmetricMatBasis(k-1,ip,mesh);
% 
e1=[1; 0];
e2=[0; 1];
%
% This matrix is useful only to check the test G=B*D
%
EpsilonHSym = zeros(dimMatBasis,dimPkVect);
for ialpha=1:dimMatBasis
    for ibeta=1:dimPkVect
        for i=1:length(xq)
			EpsilonHSym(ialpha,ibeta) = EpsilonHSym(ialpha,ibeta) + sum(sum(symmetricMatBasis(ialpha).p(xq(i),yq(i)).*grad_eps(ibeta).p(xq(i),yq(i))))*wq(i);
        end
    end
end
%
%%% ------------------------------- MATRIX B -----------------------------
%
% for the boundary dofs (vertices and edges) the philosopy is the same of
% the poisson problem. The only difference is that we consider vectorial
% dofs, for the internal (element) we need to compute the div(eps(u))
% instead of laplacian
[mono_x, mono_y, dimPkm1] = computeFirstDerivativeForMonomialBasis(k-1,ip,mesh);
%
B=zeros(dimMatBasis,NumberOfLocalDofs);
%
for ialpha=1:dimMatBasis
    %
    % vertices
    %
    for iv=1:NV
        ivm = iv-1;
        if ivm==0
            ivm =NE;
        end
        % first component
        B(ialpha,iv)=  B(ialpha,iv) +  (symmetricMatBasis(ialpha).p(x(iv),y(iv))*normal(:,iv))'*e1*ledge(iv)/2*wgl(1);
        B(ialpha,iv) = B(ialpha,iv) + (symmetricMatBasis(ialpha).p(x(iv),y(iv))*normal(:,ivm))'*e1*ledge(ivm)/2*wgl(1);
        % second component
        B(ialpha,NV+iv)=  B(ialpha,NV+iv) +  (symmetricMatBasis(ialpha).p(x(iv),y(iv))*normal(:,iv))'*e2*ledge(iv)/2*wgl(1);
        B(ialpha,NV+iv) = B(ialpha,NV+iv) + (symmetricMatBasis(ialpha).p(x(iv),y(iv))*normal(:,ivm))'*e2*ledge(ivm)/2*wgl(1);
    end
    %
    if k>=2
        %
        % Internal edge dofs
        %
        for ie=1:NE
            % ie +1
            iep = ie+1;
            for ik=1:k-1
                % internal Gauss-Lobatto point
                xik = (x(iep)-x(ie))/2*tgl(ik+1)+(x(iep)+x(ie))/2;
                yik = (y(iep)-y(ie))/2*tgl(ik+1)+(y(iep)+y(ie))/2;
                wik = wgl(ik+1)*ledge(ie)/2;
                
                % first component
                B(ialpha,2*NV+(k-1)*(ie-1)+ik)=(symmetricMatBasis(ialpha).p(xik,yik)*normal(:,ie))'*e1*wik;
                % second component
                B(ialpha,2*NV+(k-1)*NE+(k-1)*(ie-1)+ik)=(symmetricMatBasis(ialpha).p(xik,yik)*normal(:,ie))'*e2*wik;
            end
        end

    end
end
%
if k>=2
%
% Internal moments
%
for idof=1:dimP(k-2)
    multiIndex=o2m(idof);
    for ibeta=1:dimPkm1
        if ((multiIndex(1)==mono_x(ibeta).indexX) && (multiIndex(2)==mono_x(ibeta).indexY))
            % first component
            B(ibeta, 2*(NV+(k-1)*NE)+idof) = B(ibeta, 2*(NV+(k-1)*NE)+idof) - mono_x(ibeta).coeff *area;
            % second component
            B(ibeta+2*dimPkm1, 2*(NV+(k-1)*NE)+dimP(k-2)+idof) = B(ibeta+2*dimPkm1, 2*(NV+(k-1)*NE)+dimP(k-2)+idof) - mono_x(ibeta).coeff *area;
        end
        if ((multiIndex(1)==mono_y(ibeta).indexX) && (multiIndex(2)==mono_y(ibeta).indexY))
            % first component
            B(ibeta+2*dimPkm1, 2*(NV+(k-1)*NE)+idof) = B(ibeta+2*dimPkm1, 2*(NV+(k-1)*NE)+idof) - mono_y(ibeta).coeff *area;
            % second component
            B(ibeta+dimPkm1, 2*(NV+(k-1)*NE)+dimP(k-2)+idof) = B(ibeta+dimPkm1, 2*(NV+(k-1)*NE)+dimP(k-2)+idof) - mono_y(ibeta).coeff *area;
        end
    end
end

end
mat.B=B;
%
if k>=2
    [internalMoment] = generateMonomialBasis(k-2,ip,mesh);
end
%
%%% ------------------------------- MATRIX D -----------------------------
%
D=zeros(NumberOfLocalDofs,dimPkVect);
%
for ialpha=1:dimPkVect  
    % Nodal dofs
    for iv=1:NV
        D(iv,ialpha)=monoVect(ialpha).p(x(iv),y(iv))'*e1;
        D(NV+iv,ialpha)=monoVect(ialpha).p(x(iv),y(iv))'*e2;
    end
    if k>=2
        %
        % Internal edge dofs
        %
        for ie=1:NE
        % ie +1
        iep = ie+1;
            for ik=1:k-1
                % internal Gauss-Lobatto point
                xik = (x(iep)-x(ie))/2*tgl(ik+1)+(x(iep)+x(ie))/2;
                yik = (y(iep)-y(ie))/2*tgl(ik+1)+(y(iep)+y(ie))/2;
                %
                D(2*NV+(k-1)*(ie-1)+ik,ialpha)=monoVect(ialpha).p(xik,yik)'*e1;
                D(2*NV+(k-1)*NE+(k-1)*(ie-1)+ik,ialpha)=monoVect(ialpha).p(xik,yik)'*e2;
            end
        end
        % internal moment
        for idof=1:dimP(k-2)
            for i=1:length(wq)
                D(2*(NV+(k-1)*NE)+idof,ialpha)=D(2*(NV+(k-1)*NE)+idof,ialpha) + 1/area*(monoVect(ialpha).p(xq(i),yq(i)))'*e1*internalMoment(idof).p(xq(i),yq(i))*wq(i);
                D(2*(NV+(k-1)*NE)+dimP(k-2)+idof,ialpha)=D(2*(NV+(k-1)*NE)+dimP(k-2)+idof,ialpha) + 1/area*(monoVect(ialpha).p(xq(i),yq(i)))'*e2*internalMoment(idof).p(xq(i),yq(i))*wq(i);
            end      
        end

    end
end
mat.D=D;
%
%
% test EpsilonHsym - B*D  in jargon G=B*D
%
errGminusBD = max(max(abs(EpsilonHSym-B*D)));
%
if errGminusBD>1e-12
    fprintf('%s: EpsilonHsym=BD check - error is %e\n',mfilename,errGminusBD);
end
%
%
%%% ------------------------------- MATRIX H sym -------------------------
%
Hsym=zeros(dimMatBasis,dimMatBasis);
%
for ialpha=1:dimMatBasis
    for ibeta=1:dimMatBasis
        for j=1:length(wq)
			Hsym(ialpha,ibeta) = Hsym(ialpha,ibeta) + sum(sum((symmetricMatBasis(ialpha).p(xq(j),yq(j))).*symmetricMatBasis(ibeta).p(xq(j),yq(j))))*wq(j);
        end
    end
end
%
% Projection matrix
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% matrix PI Epsilon * 
mat.PiEpsilonStar = Hsym\B;
%
%%%%%%%%%%%%%%%%%%%%%% Other useful matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%% ------------------------ MATRIX HWithCTensor ------------------------
[voigtVectorBasis, dimVoigtVectorBasis] = generateVoigtVectorBasis(k-1,ip,mesh);
%
tensorC=computeElasticityTensorC(lambda, mu, true);
%
HWithCTensor=zeros(dimVoigtVectorBasis, dimVoigtVectorBasis);
%
for ialpha=1:dimVoigtVectorBasis
    for ibeta=1:dimVoigtVectorBasis
        for j=1:length(xq)
            HWithCTensor(ialpha,ibeta) = HWithCTensor(ialpha,ibeta) +...
             (tensorC*voigtVectorBasis(ialpha).p(xq(j),yq(j)))'*voigtVectorBasis(ibeta).p(xq(j),yq(j))...
               *wq(j);
        end
    end
end
mat.HWithCTensor=HWithCTensor;
%
%
%%% ------------------------ MATRIX HVect ------------------------
%
% this matrix is useful only for RHS for k>2
if k>=2
    %
    [mono_km2,dimPkm2] = generateMonomialBasis(k-2,ip,mesh);
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MATRIX H
    %
    H_km2 = zeros(dimPkm2,dimPkm2);
    %
    for ialpha=1:dimPkm2

        for ibeta=1:ialpha-1
            H_km2(ialpha,ibeta) = H_km2(ibeta,ialpha);
        end

        for ibeta=ialpha:dimPkm2    
            H_km2(ialpha,ibeta) = sum(mono_km2(ialpha).p(xq,yq).*mono_km2(ibeta).p(xq,yq).*wq);
        end
    end
    %
    mat.H_km2 = H_km2;
end
return
