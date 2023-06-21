%
for ip=1:mesh.NP
    % counter of polygon
    if mod(ip,100)==0
        fprintf('%s: polygon %d out of %d\n',mesh.name,ip,mesh.NP);
    end
    %
    % Select quadrature to integrate polynomials of degree 2k
    %
    switch mesh.POISSONVEM.QUADRATURE
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
    mat = vemLocalMatricesAndProjections_poisson(k,ip,mesh);
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
    %
    [mono, dimPk] = generateMonomialBasis(k,ip,mesh);
    %
    %
    %% number of local dofs -> need to define
    %
    Ndof_local = NV + (k-1)*NE + dimP(k-2);
    %
    % local dofs 
    %
    local_dofs = 1:Ndof_local;
    %
    %% global dofs of polygon ip -> need to define 
    %
    global_dofs = get_global_dofs_poisson(k,ip,mesh);
    %
    % identity matrix of dimension Ndof x Ndof (local)
    %
    Id = eye(Ndof_local,Ndof_local);
    %
    %% CONSISTENCY MATRIX FOR THE DIFFUSION TERM -> need to define
    %
    Kc = mat.PNkstar'*mat.Gt*mat.PNkstar;
    %
    % force symmetry
    %
    Kc = (Kc+Kc')/2;
    %
    %% STABILITY MATRIX -> need to define
    %
    Ks = (Id-mat.PNk)'*Id*(Id-mat.PNk);
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
    for i=1:Ndof_local
        %
        % computing the L^2 projection 
        % of the basis function phi_i
        % in the quadrature points
        %
        P0kphiq = 0;
        %
        for ialpha=1:dimPk
            P0kphiq = P0kphiq + mat.P0kstar(ialpha,i)*mono(ialpha).p(xq,yq);
        end
        %
        bE(i) = sum(P0kphiq.*f(xq,yq).*wq);
        %
    end
    %
    %%% ASSEMBLING THE RIGHT-HAND_SIDE
    %
    b(global_dofs) = b(global_dofs) + bE;
    %
    % dump projectors on the mesh for error computation
    %
    mesh.polygon(ip).P0kstar = mat.P0kstar;
    mesh.polygon(ip).PNkstar = mat.PNkstar;
    %
end