function global_dofs=get_global_dofs_ezvem(k,ip,mesh,type)
%%%
%%% This functions returns the array of global dofs of polygon ip
%%%
%
% We start by computing boundary dofs
%
% vertices dofs
%
boundary_global_dofs = mesh.polygon(ip).vertices;
%
% k-1 dofs on each edge
%
if k>=2
    %
    % loop on edges
    %
    for ie=1:mesh.polygon(ip).NE
        %
        % global number of the edge
        %
        ie_global = mesh.polygon(ip).edges(ie);
        %
        % number of edge dofs of previous edges
        %
        Ndof_previous_edges = (ie_global-1)*(k-1);
        %
        % check orientation of edge with respect to polygon
        %
        if mesh.polygon(ip).vertices(ie) == mesh.edge(ie_global).v1
            %
            % same orientation
            %
            edge_global_dofs = mesh.NV+Ndof_previous_edges+(1:k-1);
            %
        elseif mesh.polygon(ip).vertices(ie) == mesh.edge(ie_global).v2
            %
            % opposite orientation
            %
            edge_global_dofs = mesh.NV+Ndof_previous_edges+(k-1:-1:1);
            %
        else
            %
            % something was wrong!!
            %
            fprintf('%s: something was wrong!',mfilename)
            %
        end
        %
        % update boundary dofs
        %
        boundary_global_dofs = [boundary_global_dofs edge_global_dofs];
        %
    end
end
%
% dimP(k-2) internal dofs on ech polygon
%
if k>=2
    %
    % dimP(k-2) for each polygon
    %
    % number of internal dofs of previous polygons
    %
    Ndof_previous_polygons = (ip-1)*dimP(k-2);
    %
    % internal dofs of polygon ip
    %    
    internal_global_dofs = mesh.NV + mesh.NE*(k-1) + Ndof_previous_polygons+(1:dimP(k-2));
    %
else
    %
    % there are no internal dofs for k01
    %
    internal_global_dofs = [];
    %
end
%
if nargin == 3
    %
    % default: return all dofs
    %
    global_dofs = [boundary_global_dofs internal_global_dofs];
    %
elseif nargin == 4
    %
    % switch between boundary and internal dofs
    %
    switch type
        case 'boundary'
            global_dofs = boundary_global_dofs;
        case 'internal'
            global_dofs = internal_global_dofs;
        otherwise
            fprintf('%s: bad type argument: %s',mfilename,type)
            return
    end
else
    fprintf('%s: bad number of arguments',mfilename)
end
%
end