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
fprintf('%s: computing various errors Elasticity\n',mfilename);
%
% exact solution at the nodal dofs (on the boundary)
%
%first component
ueh1 = ue1(xn,yn);
%second component
ueh2 = ue2(xn,yn);
%
% ordering 
%
v1=1:mesh.NV;
v2=[];
 if k>=2
     v2=mesh.NV+1:mesh.NV+(k-1)*mesh.NE;
 end
% exact solution with the same order of discrete solution
ueh=[ueh1(v1); ueh2(v1); ueh1(v2); ueh2(v2)];
uehV=[ueh1(v1); ueh2(v1)];
%
% error on the nodal dofs
%
% l2 error on the nodal dofs
%
errors.errl2=norm(ueh - uh(1:2*mesh.NV+2*(k-1)*mesh.NE));
errors.norml2ue=sqrt(sum(2*ueh.*ueh));
%
% L2 norm of the exact Sigma solution
%
normL2Epsilonsq = 0;
%
% L2 norm with respect to \PEspilon\Sigmah
%
normL2PiEpsilonsq=0;
%
% L2 norm with respect to \PEspilon\Sigmah
%
normL2Sigma_hsq =0;
%
% L2 norm with respect to \PEspilon\Sigmah
%
normL2Sigmasq =0;
%
%
%
tensorC=computeElasticityTensorC(lambda,mu,false);
%
%
for ip=1:mesh.NP
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
    global_dofs = get_global_dofs_elasticity(k,ip,mesh);
    %
    PiEpsilonStar=mesh.polygon(ip).PiEpsilonStar;
    %
    % ||epsilon(u)-Pi(epsilon(uh))||_0
    % 
    [symmetricMatBasis, dimMatBasis] = generateSymmetricMatBasis(k-1,ip,mesh);
    for i=1:length(wq)
        vemSolution=zeros(2,2);
        for idof=1:size(PiEpsilonStar,2)
            tmpVal=zeros(2,2);
            for irow=1:size(PiEpsilonStar,1)
                tmpVal=tmpVal+PiEpsilonStar(irow,idof)*symmetricMatBasis(irow).p(xq(i),(yq(i)));
            end
            vemSolution = vemSolution + tmpVal*uh(global_dofs(idof));
        end
        exactSolMat = epsilonU(xq(i),yq(i));
        normEpsilonU = norm(exactSolMat);
        normDiff = norm(exactSolMat-vemSolution);        
        normL2Epsilonsq = normL2Epsilonsq + normEpsilonU*normEpsilonU*wq(i);
        normL2PiEpsilonsq =  normL2PiEpsilonsq +normDiff*normDiff*wq(i);
    end
    %
    % ||sigma(u)-C Pi(epsilon(uh))||_0
    % 
    [voigtVectorBasis, dimVoigtVectorBasis] = generateVoigtVectorBasis(k-1,ip,mesh);
    for i=1:length(wq)
        vemSolution=zeros(3,1);
        for idof=1:size(PiEpsilonStar,2)
            tmpVal=zeros(3,1);
            for irow=1:size(PiEpsilonStar,1)
                tmpVal=tmpVal+PiEpsilonStar(irow,idof)*tensorC*voigtVectorBasis(irow).p(xq(i),(yq(i)));
            end
            tmpVal = tmpVal*uh(global_dofs(idof));
            vemSolution = vemSolution + tmpVal;
        end
        exactSolMat = epsilonU(xq(i),yq(i));
        exactSolVect = [exactSolMat(1,1); exactSolMat(2,2);exactSolMat(1,2)];
        exactSol = tensorC*exactSolVect;
        normSigma = norm(exactSol);
        normDiffSigma = norm(exactSol-vemSolution);
            
        normL2Sigmasq = normL2Sigmasq + normSigma*normSigma*wq(i);
        normL2Sigma_hsq =  normL2Sigma_hsq +normDiffSigma*normDiffSigma*wq(i);
    end
end

errors.normL2Epsilonsq = sqrt(normL2Epsilonsq);
errors.errL2PiEpsilonsq = sqrt(normL2PiEpsilonsq);

errors.normL2Sigmasq = sqrt(normL2Sigmasq);
errors.errL2Sigma_hsq = sqrt(normL2Sigma_hsq);

%
mesh.errors = errors;
%
