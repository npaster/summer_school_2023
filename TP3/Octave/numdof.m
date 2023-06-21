function numdof(mesh,k)
%
% numeration of global dofs
%
if k>1
    [tgl,~] = mylglnodes(k);
end
%
% vertex and edge nodes (boundary of polygons)
%
Ndof_boundary = mesh.NV+(k-1)*mesh.NE;
%
xn = zeros(Ndof_boundary,1);
yn = zeros(Ndof_boundary,1);
%
% vertices
%
xn(1:mesh.NV) = [mesh.vertex.x];
yn(1:mesh.NV) = [mesh.vertex.y];
%
% edge nodes
%
for ieg=1:mesh.NE
    %
    % extrema of the edge
    %
    x1 = mesh.vertex(mesh.edge(ieg).v1).x;
    y1 = mesh.vertex(mesh.edge(ieg).v1).y;
    x2 = mesh.vertex(mesh.edge(ieg).v2).x;
    y2 = mesh.vertex(mesh.edge(ieg).v2).y;
    %
    % internal edge nodes
    %
    for ik=1:k-1
        %
        % global index of node
        %
        ig = mesh.NV + (ieg-1)*(k-1)+ik;
        %
        % Gauss-Lobatto points on the edge
        %
        xik = (x2-x1)/2*tgl(ik+1)+(x2+x1)/2;
        yik = (y2-y1)/2*tgl(ik+1)+(y2+y1)/2;
        %
        xn(mesh.NV+(ieg-1)*(k-1)+ik) = xik;
        yn(mesh.NV+(ieg-1)*(k-1)+ik) = yik;
        
    end    
end
%
str= '';
for iv=1:Ndof_boundary
  str = strvcat(str,int2str(iv));
end
%
text(xn,yn,1e-3*ones(Ndof_boundary,1),str,'Color','b','BackgroundColor','w','HorizontalAlignment','center')
%
if k>1
    %
    for ip=1:mesh.NP
        str = [];
        for j=1:dimP(k-2)
            dof = Ndof_boundary + (ip-1)*dimP(k-2) + j;
            str = strvcat(str,int2str(dof));
        end
        text(mesh.polygon(ip).xb,mesh.polygon(ip).yb,1.e-3,str,'Color','b','BackgroundColor','w','HorizontalAlignment','center')
    end
end

    






