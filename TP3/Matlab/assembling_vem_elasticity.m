% 
for ip=1:mesh.NP
    % counter of polygon
    if mod(ip,100)==0
        fprintf('%s: polygon %d out of %d\n',mesh.name,ip,mesh.NP);
    end
    %
    % Select quadrature to integrate polynomials of degree 2k
    %
    switch mesh.ELASTICITYVEM.QUADRATURE
        case 'subtriangulation'
            [xq,yq,wq] = quadrature_subtriangulation(2*k,ip,mesh);
        case 'vianello'
            [xq,yq,wq] = quadrature_vianello(2*k,ip,mesh);
        otherwise
            fprintf('%s: error in quadrature method\n',mfilename);
            return
    end
    %
    % storing quadrature nodes and weights in the mesh structure
    %
    mesh.polygon(ip).xq = xq;
    mesh.polygon(ip).yq = yq;
    mesh.polygon(ip).wq = wq;
    %
    % main call to vemprojectors
    %
    mat = vemLocalMatricesAndProjections_elasticity(k,ip,mesh);
    %
    % number of edges and vertices of the polygon
    %
    NE = mesh.polygon(ip).NE;
    NV = mesh.polygon(ip).NV;
    %
    % centroid
    %
    xb = mesh.polygon(ip).xb;
    yb = mesh.polygon(ip).yb;
    %
    % area
    %
    area = mesh.polygon(ip).area;
    %
    % number of local dofs
    %
    Ndof_local = 2*(NV+(k-1)*NE+dimP(k-2));
    %
    % local dofs
    %
    local_dofs = 1:Ndof_local;
    %
    % global dofs of polygon ip
    %
    global_dofs = get_global_dofs_elasticity(k,ip,mesh);
    %
    % identity matrix of dimension Ndof x Ndof (local)
    %
    Id = eye(Ndof_local,Ndof_local);
    %
    %%% CONSISTENCY MATRIX FOR THE DIFFUSION TERM
    %
    Kc = mat.PiEpsilonStar'*mat.HWithCTensor*mat.PiEpsilonStar;
    %
    % force symmetry
    %
     Kc = (Kc+Kc')/2;
    %
    %%% STABILITY MATRIX
    %
    %Ks = (Id-(mat.D*inv(mat.D'*mat.D))*mat.D');
    %
    Ks = (Id-(mat.D*inv(mat.D'*mat.D))*mat.D');
    tau=1/2;
    traceKc=0;
    for i=1:length(Kc)
            traceKc = traceKc + Kc(i,i);
    end
    %Ks = Ks*tau*traceKc;
        %%% D recipe stabilization
%             Ks1 = (Id-(mat.D*inv(mat.D'*mat.D))*mat.D');
%             DRecipe=zeros(size(Ks1));
%             for i=1:length(DRecipe)
%                 DRecipe(i,i) = max(Kc(i,i),1);
%             end
%             Ks=Ks1'*DRecipe*Ks1;
    %
    % force symmetry
    %
     Ks = (Ks+Ks')/2;
    %
    % save matrices on the mesh
    %
    mesh.polygon(ip).Kc = Kc;
    mesh.polygon(ip).Ks = Ks;
    %
    % assembling
    %
    A(global_dofs,global_dofs) = A(global_dofs,global_dofs) + Kc + Ks;
    %
    % right hand side \int f \pizero_k\phi_i
    %
    bE = zeros(Ndof_local,1);
    %
    % integrate right-hand-side
    %
    if k==1
        for idof=1:NV % <---------------------------------
            %
            % computing the L^2 projection 
            % of the basis function phi_i
            % in the quadrature points
            %
            % first component
            bE(idof) = sum(f1(xq,yq,lambda,mu).*wq)/NV;
            % second component
            bE(NV+idof) = sum(f2(xq,yq,lambda,mu).*wq)/NV;
            %
        end
    elseif k==2
        bE(2*NV+2*(k-1)*NE+1)=sum(f1(xq,yq,lambda,mu).*wq);
        bE(2*NV+2*(k-1)*NE+2)=sum(f2(xq,yq,lambda,mu).*wq);
    else
        H_km2=mat.H_km2;
        [mono, dimPkm2]=generateMonomialBasis(k-2,ip,mesh);
        for ialpha=1:dimPkm2
            b1(ialpha,1) =  sum(mono(ialpha).p(xq,yq).*f1(xq,yq,lambda,mu).*wq);
            b2(ialpha,1) =  sum(mono(ialpha).p(xq,yq).*f2(xq,yq,lambda,mu).*wq);
        end
        bProjectionCoeff1= H_km2\b1;
        bProjectionCoeff2= H_km2\b2;
        
        bProjectionCoeff= [bProjectionCoeff1;bProjectionCoeff2];
        
        for idof=1:length(bProjectionCoeff)
            bE(2*NV+2*(k-1)*NE+idof) = bProjectionCoeff(idof)*area;
        end
    end
    %
    %%% ASSEMBLING THE RIGHT-HAND_SIDE
    %
    b(global_dofs) = b(global_dofs) + bE;
    %
    % dump projectors on the mesh for error computation
    %
    mesh.polygon(ip).PiEpsilonStar = mat.PiEpsilonStar;
    %
end