%
% array with free global dofs
%
free_global_dofs = [];
%
% boundary condition
free_global_dofs_vertices1=[];
free_global_dofs_vertices2=[];
%
% Dirichlet bc on vertices - update load term
%
for iv=1:mesh.NV
    if mesh.vertex(iv).marker == 1
        %
        % boundary vertex
        %
        % first component
        %
        uh(iv) = g1(mesh.vertex(iv).x,mesh.vertex(iv).y);
        %
        % update rhs
        %
        b = b - uh(iv)*A(:,iv);
        %
        % second component
        %
        uh(mesh.NV+iv) = g2(mesh.vertex(iv).x,mesh.vertex(iv).y);
        %
        % update rhs
        %
        b = b - uh(mesh.NV+iv)*A(:,mesh.NV+iv);
        %
    else
        %
        % internal vertex - update free global dofs
        %
        free_global_dofs_vertices1 = [free_global_dofs_vertices1  iv];
        free_global_dofs_vertices2 = [free_global_dofs_vertices2  mesh.NV+iv];
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
xn = zeros(mesh.NV,1);
yn = zeros(mesh.NV,1);
%
% vertices 
%
xn(1:mesh.NV) = [mesh.vertex.x];
yn(1:mesh.NV) = [mesh.vertex.y];
% idof=2*mesh.NV+1;
free_global_dofs_edges1=[];
free_global_dofs_edges2=[];
idof=1;
for ieg=1:mesh.NE
    %     %
    %     % extrema of the edge globally ordered
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
                % first component
                %
                uh(2*mesh.NV+idof) = g1(xik,yik);
                %
                % update rhs
                %
                b = b - uh(2*mesh.NV+idof)*A(:,2*mesh.NV+idof);
                %
                % second component 
                %
                uh(2*mesh.NV+(k-1)*mesh.NE+idof) = g2(xik,yik);
                %
                % update rhs
                %
                b = b - uh(2*mesh.NV+(k-1)*mesh.NE+idof)*A(:,2*mesh.NV+(k-1)*mesh.NE+idof);
                %
            else
                %
                % internal or Neumann boundary condition - ik_global is a free node
                %
                free_global_dofs_edges1 = [free_global_dofs_edges1 2*mesh.NV+idof];
                free_global_dofs_edges2 = [free_global_dofs_edges2 2*mesh.NV+(k-1)*mesh.NE+idof];
                %
            end
            idof=idof+1;
        end 
end
free_global_dofs = [free_global_dofs_vertices1 free_global_dofs_vertices2 free_global_dofs_edges1 free_global_dofs_edges2];
% adding all internal global dofs to the free global dofs
internal_global_dofs = 2*mesh.NV + 2*(k-1)*mesh.NE + (1:2*dimP(k-2)*mesh.NP);
%
free_global_dofs = [free_global_dofs internal_global_dofs];