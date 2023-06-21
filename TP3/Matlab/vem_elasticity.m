function mesh=vem_elasticity(k,lambda, mu, mesh,testcase)
%
%    MESH=ELASTICITYVEM(K,MESH,TESTCASE)
%
%    K = degree of approximation
%    MESH = mesh to be considered
%    TESTCASE = actual PDE to be approximated
%      other choices are made in the code
%
%    output:
%    MESH.ELASTICITYVEM given options
%    MESH.errors contain all errors
%    MESH.uh     is the solution
%      various computed quantities are attached to MESH
%

%
% check if called for convergence analysis or standalone
%
caller = dbstack;
if  strcmp(caller(end).file,'go_elasticity_vem.m')
    ELASTICITYVEM.STANDALONE = false;
else
    close all
    ELASTICITYVEM.STANDALONE = true;
end
%
%%%%%%%%%%%%%%%%%%%%%%% QUADRATURE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ELASTICITYVEM.QUADRATURE = 'vianello';
%ELASTICITYVEM.QUADRATURE = 'subtriangulation';
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%% TESTCASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if not(exist('testcase','var'))
%     testcase = 'test1';
%     testcase = 'test2';
%     testcase = 'test3';
     testcase = 'p2';
end
%
set_testcase_elasticity(testcase)
%
ELASTICITYVEM.TESTCASE = testcase;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VEM DEGREE %%%%%%%%%%%%%%%%%%%%%%
%
if not(exist('k','var'))
    %
    k = 2;
    %
end
%
ELASTICITYVEM.K = k;
%
%%%%%%%%%%%%%%%%%%%%%%%% LAME' COEFFICIENTS %%%%%%%%%%%%%%%%%%%%%%
%
if not(exist('lambda','var'))
    %
    lambda = 1;
    %
end
%
ELASTICITYVEM.lambda=lambda;
%
if not(exist('mu','var'))
    %
    mu = 1;
    %
end
%
ELASTICITYVEM.mu=mu;
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MESH BEGINS %%%%%%%%%%%%%%%%%%%%%%

if not(exist('mesh','var'))
    
    % STRUCTURED SQUARES
    N=8;
    mesh = square2mesh('s',N);
    %mesh = square2mesh('s',Nx,Ny);
    
    % STRUCTURED TRIANGLES
    
    %mesh = square2mesh('t',N);
    %mesh = square2mesh('t',Nx,Ny);

    % UNSTRCTURED TRIANGLES
    
    %mesh = square2mesh('u',area);    
    
    % DISTORTED SQUARES
    
    %mesh = square2mesh('dd-s',Nx,Ny,0<distortion<0.5);

    % RANDOM POLYGONS
    
    %mesh = square2mesh('v',NP,0);
    
    % REGULAR POLYGONS (LLOYD)
    
    %mesh = square2mesh('v',NP);    
        
    % (ALMOST) REGULAR HEXAGONS
    
    %mesh = square2mesh('hexagons-regular/hexagon-8x10');
    %mesh = square2mesh('hexagons-regular/hexagon-18x20');
    %mesh = square2mesh('hexagons-regular/hexagon-26x30');
    %mesh = square2mesh('hexagons-regular/hexagon-34x40');
    %mesh = square2mesh('hexagons-regular/hexagon-44x50');
    %mesh = square2mesh('hexagons-regular/hexagon-52x60');
    %mesh = square2mesh('hexagons-regular/hexagon-60x70');
    %mesh = square2mesh('hexagons-regular/hexagon-70x80');
    
    % DEFORMED HEXAGONS
    
    %mesh = square2mesh('hexagons-deformed/hexagon-8x10d');
    %mesh = square2mesh('hexagons-deformed/hexagon-18x20d');
    %mesh = square2mesh('hexagons-deformed/hexagon-26x30d');
    %mesh = square2mesh('hexagons-deformed/hexagon-34x40d');
    %mesh = square2mesh('hexagons-deformed/hexagon-44x50d');
    %mesh = square2mesh('hexagons-deformed/hexagon-52x60d');
    %mesh = square2mesh('hexagons-deformed/hexagon-60x70d');
    %mesh = square2mesh('hexagons-deformed/hexagon-70x80d');
    
    % PEGASUS
    
    %mesh = square2mesh('pegasus/pegasus-1x1');
    %mesh = square2mesh('pegasus/pegasus-2x2');
    %mesh = square2mesh('pegasus/pegasus-4x4');
    %mesh = square2mesh('pegasus/pegasus-8x8');
    %mesh = square2mesh('pegasus/pegasus-10x10');
    %mesh = square2mesh('pegasus/pegasus-16x16');
    
    % VEM
    
    %mesh = square2mesh('VEM/VEM-1x1');
    %mesh = square2mesh('VEM/VEM-2x2');
    %mesh = square2mesh('VEM/VEM-4x4');
    %mesh = square2mesh('VEM/VEM-6x6');
    %mesh = square2mesh('VEM/VEM-10x10');
    %mesh = square2mesh('VEM/VEM-12x12');
    %mesh = square2mesh('VEM/VEM-14x14');
    %mesh = square2mesh('VEM/VEM-16x16');
    
    % OTHER MESHES
    
    %mesh = mesh2concave(square2mesh('v',100));
    %mesh = mesh2concave(square2mesh('v',1000,0));
    
end
% add the name of the mesh 
ELASTICITYVEM.MESHNAME = mesh.name;
mesh.ELASTICITYVEM = ELASTICITYVEM;
%
% drawing mesh
%
figure(1)
clf
drawmesh(mesh)
drawnow
view(2)
%
if ELASTICITYVEM.STANDALONE
    movegui('northwest');
end
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE STARTS HERE %%%%%%%%%%%%%%%%%%%%%
%
% total number of global dofs (including boundary dofs)
%
Ndof_global = 2*mesh.NV + 2*(k-1)*mesh.NE+2*dimP(k-2)*mesh.NP;
%
% stiffness matrix, right-hand-side, solution
%
A = sparse(Ndof_global,Ndof_global);
b = zeros(Ndof_global,1);
uh = zeros(Ndof_global,1);
%
fprintf('%s: info: elasticity problem with lambda = %f and mu = %f\n',mfilename,lambda, mu)
fprintf('%s: info: assembling stiffness matrix\n',mfilename)
%
assembling_vem_elasticity
%
fprintf('%s: info: enforcing boundary conditions\n',mfilename)
%
%
enforcing_BC_vem_elasticity
%
% extract submatrices corresponding to the free dofs
%
Kh = A(free_global_dofs,free_global_dofs);
fh = b(free_global_dofs);
%
fprintf('%s: info: solving linear system\n',mfilename)
%
uh(free_global_dofs) = Kh\fh;
%
% saving solution on the mesh
%
mesh.uh = uh;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% COMPUTING ERRORS
%
compute_errors_elasticity
%
