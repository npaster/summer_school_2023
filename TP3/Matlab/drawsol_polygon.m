function drawsol_polygon(mesh)
%%%
%%% draws the L2 projection of the VEM solution
%%%
k = mesh.POISSONVEM.K;
uh = mesh.uh;
%
for ip=1:mesh.NP
    %
    %%% monomial scaled basis
    %
    [mono, ~] = generateMonomialBasis(k,ip,mesh);
    %
    % samples per edge
    %
    NSe = 20;
    %
    % global dofs of the polygon
    %
    global_dofs = get_global_dofs_poisson(k,ip,mesh);  
    %
    % coefficients of m_alpha of the two projections of uh
    %
    sP0k = mesh.polygon(ip).P0kstar*uh(global_dofs);
    sPNk = mesh.polygon(ip).PNkstar*uh(global_dofs);
    %
    % array for the patch command on polygon ip
    %
    Xip = [];
    Yip = [];
    Zip = [];
    %
    for iv=1:mesh.polygon(ip).NV
        %
        ivp = iv+1;
        if ivp == mesh.polygon(ip).NV+1
            ivp = 1;
        end
        %
        x1 = mesh.vertex(mesh.polygon(ip).vertices(iv)).x;
        y1 = mesh.vertex(mesh.polygon(ip).vertices(iv)).y;
        %
        x2 = mesh.vertex(mesh.polygon(ip).vertices(ivp)).x;
        y2 = mesh.vertex(mesh.polygon(ip).vertices(ivp)).y;
        %
        xe = linspace(x1,x2,NSe);
        ye = linspace(y1,y2,NSe);
        %
        P0kuhe = 0;
        %
        for ialpha=1:dimP(k)
            P0kuhe = P0kuhe + sP0k(ialpha)*mono(ialpha).p(xe,ye);
        end
        %
        Xip = [Xip; xe(1:end-1)'];
        Yip = [Yip; ye(1:end-1)'];
        %
        % uncomment to draw the error
        %
        %Zip = [Zip; P0kuhe(1:end-1)'-ue(xv(1:end-1),yv(1:end-1))'];
        %
        Zip = [Zip; P0kuhe(1:end-1)'];
        %
    end
    %
    patch(Xip,Yip,Zip,'w')
    %
end
