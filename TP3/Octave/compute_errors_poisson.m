%%%
%%% COMPUTE VARIOUS ERRORS
%%%
%
% size of the mesh
%
errors.hmax = max([mesh.polygon.h]);
errors.hmean = sum([mesh.polygon.h])/mesh.NP;
%
% VEM degree
%
errors.k = k;
%
fprintf('%s: computing various errors\n',mfilename);
%
% exact solution at the nodal dofs (on the boundary)
%
ueh = ue(xn,yn);
%
% error on the nodal dofs
%
% l2 error on the nodal dofs
%
errors.errl2 = norm(ueh-uh(1:mesh.NV+(k-1)*mesh.NE));
errors.norml2ue = norm(ueh);
%
% maximum error on the nodal dofs
%
errors.errmax = max(abs(ueh-uh(1:mesh.NV+(k-1)*mesh.NE)));
errors.normmaxue = max(ueh);
%
% maximum error on vertices only
%
errors.errmaxV = max(abs(ueh(1:mesh.NV)-uh(1:mesh.NV)));
errors.normmaxVue = max(ueh(1:mesh.NV));
%
% L2 norm of the exact solution
%
normL2uesq = 0;
%
% H1 norm of the exact solution
%
normH1uesq = 0;
%
% L2 and H1 error with respect to \P0k\uh
%
errL2P0ksq = 0;
errH1P0ksq = 0;
%
% L2 and H1 error with respect to \PNk\uh
%
errL2PNksq = 0;
errH1PNksq = 0;
%
% maximun error on the nodal dofs of \P0k\uh
%
errmaxP0k = 0;
%
% maximun error on the nodal dofs of \PNk\uh
%
errmaxPNk = 0;
%
% loop on polygons
%
for ip=1:mesh.NP
    %
    %%% monomial scaled basis
    %
    [mono, dimPk] = generateMonomialBasis(k,ip,mesh);
    %
    %%% first derivative of the monomial scaled basis
    %
    [mono_x, mono_y] = computeFirstDerivativeForMonomialBasis(k,ip,mesh);
    %
    NV = mesh.polygon(ip).NV;
    %
    % quadrature nodes and weigths on the polygon
    %
    xq = mesh.polygon(ip).xq;
    yq = mesh.polygon(ip).yq;
    wq = mesh.polygon(ip).wq;
    %
    % global dofs of the polygon
    %
    global_dofs = get_global_dofs_poisson(k,ip,mesh);  
    %
    % boundary global dofs
    %
    boundary_global_dofs = get_global_dofs_poisson(k,ip,mesh,'boundary');  
    %
    Nbdof = length(boundary_global_dofs);
    %
    if Nbdof ~= k*NV
        fprintf('%s: wrong bgdof!\n',mfilename)
        return
    end
    %
    % coefficients of m_alpha of the two projections of uh
    %
    sP0k = mesh.polygon(ip).P0kstar*uh(global_dofs);
    sPNk = mesh.polygon(ip).PNkstar*uh(global_dofs);
    %
    % P0kuh e PNkuh nei nodi di quadratura
    %
    P0kuhq = 0;
    P0kuhxq = 0;
    P0kuhyq = 0;
    %
    PNkuhq = 0;
    PNkuhxq = 0;
    PNkuhyq = 0;
    %
    % P0kuh e PNkuh nei boundary dofs
    %
    P0kuhnip = 0;
    PNkuhnip = 0;
    %
    for ialpha=1:dimP(k)
        %
        P0kuhq = P0kuhq + sP0k(ialpha)*mono(ialpha).p(xq,yq);
        %
        P0kuhxq = P0kuhxq + sP0k(ialpha)*mono_x(ialpha).p(xq,yq);
        P0kuhyq = P0kuhyq + sP0k(ialpha)*mono_y(ialpha).p(xq,yq);
        %
        PNkuhq = PNkuhq + sPNk(ialpha)*mono(ialpha).p(xq,yq);
        %
        PNkuhxq = PNkuhxq + sPNk(ialpha)*mono_x(ialpha).p(xq,yq);
        PNkuhyq = PNkuhyq + sPNk(ialpha)*mono_y(ialpha).p(xq,yq);
        %
    end
    %
    % errori L2 squared
    %
    errH1P0ksq = errH1P0ksq + sum(((P0kuhxq-uex(xq,yq)).^2+(P0kuhyq-uey(xq,yq)).^2).*wq);
    errH1PNksq = errH1PNksq + sum(((PNkuhxq-uex(xq,yq)).^2+(PNkuhyq-uey(xq,yq)).^2).*wq);
    %
    errL2P0ksq = errL2P0ksq + sum((P0kuhq-ue(xq,yq)).^2.*wq);
    errL2PNksq = errL2PNksq + sum((PNkuhq-ue(xq,yq)).^2.*wq);
    %
    % norma L2 squared della soluzione esatta
    %
    normL2uesq = normL2uesq + sum(ue(xq,yq).^2.*wq);
    %
    % norma H1 squared della soluzione esatta
    %
    normH1uesq = normH1uesq + sum((uex(xq,yq).^2+uey(xq,yq).^2).*wq);
    %
end
%
% estraggo le radici quadre
%
errors.normL2ue = sqrt(normL2uesq);
errors.normH1ue = sqrt(normH1uesq);
%
errors.errL2P0k = sqrt(errL2P0ksq);
errors.errH1P0k = sqrt(errH1P0ksq);
%
errors.errL2PNk = sqrt(errL2PNksq);
errors.errH1PNk = sqrt(errH1PNksq);
%
errors.errL2 = errors.errL2P0k;
errors.errH1 = errors.errH1P0k;
%
mesh.errors = errors;
%
