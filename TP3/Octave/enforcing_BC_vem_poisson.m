%
% array with free global dofs
%
free_global_dofs = [];
%
% boundary condition
%
% Dirichlet bc on vertices - update load term
%
for iv=1:mesh.NV
    if mesh.vertex(iv).marker == 1
        %
        % boundary vertex
        %
        uh(iv) = g(mesh.vertex(iv).x,mesh.vertex(iv).y);
        %
        % update rhs
        %
        b = b - uh(iv)*A(:,iv);
        %
    else
        %
        % internal vertex - update free global dofs
        %
        free_global_dofs = [free_global_dofs iv];
        %
    end
end
%
% tgl,wgl = Gauss-Lobatto quadrature nodes and weights on [-1 +1]
%
[tgl,~] = mylglnodes(k);
%
% array for the nodal dofs (vertices and edges) 
%
xn = zeros(mesh.NV+(k-1)*mesh.NE,1);
yn = zeros(mesh.NV+(k-1)*mesh.NE,1);
%
% vertices 
%
xn(1:mesh.NV) = [mesh.vertex.x];
yn(1:mesh.NV) = [mesh.vertex.y];
%
for ieg=1:mesh.NE
    %
    % extrema of the edge globally ordered
    %
    x1 = mesh.vertex(mesh.edge(ieg).v1).x;
    y1 = mesh.vertex(mesh.edge(ieg).v1).y;
    %
    x2 = mesh.vertex(mesh.edge(ieg).v2).x;
    y2 = mesh.vertex(mesh.edge(ieg).v2).y;
    %
    % loop on the internal nodes of the edge
    %
    for ik=1:k-1
        %
        % global index of the internal node
        %
        ik_global = mesh.NV + (ieg-1)*(k-1)+ik;
        %
        % Gauss-Lobatto points on the edge
        %
        xik = (x2-x1)/2*tgl(ik+1)+(x2+x1)/2;
        yik = (y2-y1)/2*tgl(ik+1)+(y2+y1)/2;
        %
        xn(mesh.NV+(ieg-1)*(k-1)+ik) = xik;
        yn(mesh.NV+(ieg-1)*(k-1)+ik) = yik;
        %
        if mesh.edge(ieg).marker == 1
            %
            % Dirichlet boundary condition on the edge - update load term
            %
            uh(ik_global) = g(xik,yik);
            b = b - uh(ik_global)*A(:,ik_global);
            %
        else
            %
            % internal or Neumann boundary condition - ik_global is a free node
            %
            free_global_dofs = [free_global_dofs ik_global];
            %
        end
    end    
end
%
% adding all internal global dofs to the free global dofs
%
internal_global_dofs = mesh.NV + (k-1)*mesh.NE + (1:dimP(k-2)*mesh.NP);
%
free_global_dofs = [free_global_dofs internal_global_dofs];
