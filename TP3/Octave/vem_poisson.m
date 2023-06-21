function mesh=vem_poisson(k,mesh,testcase)
%
%    MESH=POISSONVEM(K,MESH,TESTCASE)
%
%    K = degree of approximation
%    MESH = mesh to be considered
%    TESTCASE = actual PDE to be approximated
%      other choices are made in the code
%
%    output:
%    MESH.POISSONVEM given options
%    MESH.errors contain all errors
%    MESH.uh     is the solution
%      various computed quantities are attached to MESH
%

%
% check if called for convergence analysis or standalone
%
caller = dbstack;
if  strcmp(caller(end).file,'go_vem_poisson.m')
    POISSONVEM.STANDALONE = false;
else
    close all
    POISSONVEM.STANDALONE = true;
end
%
%%%%%%%%%%%%%%%%%%%%%%% QUADRATURE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%
POISSONVEM.QUADRATURE = 'vianello';
%POISSONVEM.QUADRATURE = 'subtriangulation';
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%% TESTCASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
if not(exist('testcase','var'))
    % testcase = 'p1';
    % testcase = 'p2';
    % testcase = 'p3';
     testcase = 'poisson';
    % testcase = 'paraboloid';
end
%
%set_testcase(testcase); % you can use if you have symbolic toolbox in octave
[ue,uex,uey,f,g]=set_testcase_octave(testcase);
%
POISSONVEM.TESTCASE = testcase;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VEM DEGREE %%%%%%%%%%%%%%%%%%%%%%
%
if not(exist('k','var'))
    %
    k = 1;
    %
end
%
POISSONVEM.K = k;
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MESH BEGINS %%%%%%%%%%%%%%%%%%%%%%

if not(exist('mesh','var'))
    
    % STRUCTURED SQUARES
    %mesh = square2mesh('s',N);
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
    
    mesh = square2mesh('hexagons-regular/hexagon-8x10');
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
POISSONVEM.MESHNAME = mesh.name;
%
% fill up the quantities of the problem
%
mesh.POISSONVEM = POISSONVEM;
%
% drawing mesh
%
figure(1)
clf
drawmesh(mesh)
drawnow
view(2)
%
if POISSONVEM.STANDALONE
    movegui('northwest');
end

%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CODE STARTS HERE %%%%%%%%%%%%%%%%%%%%%
%
% total number of global dofs (including boundary dofs)
%
Ndof_global = mesh.NV + (k-1)*mesh.NE + dimP(k-2)*mesh.NP;
%
% stiffness matrix, right-hand-side, solution
%
A = sparse(Ndof_global,Ndof_global);
b = zeros(Ndof_global,1);
uh = zeros(Ndof_global,1);
%
fprintf('%s: info: assembling stiffness matrix\n',mfilename)
% 
assembling_vem_poisson;
%
fprintf('%s: info: enforcing boundary conditions\n',mfilename)
%
enforcing_BC_vem_poisson;
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
compute_errors_poisson
%
% draw solution on edges
%
figure(2)
clf
drawsol_edge(mesh,ue)
view(3)
%
title(sprintf('VEM solution on edges, k=%d',k))
drawnow
%
if POISSONVEM.STANDALONE
    movegui('north')
end
%
% draw L^2 projectionb of the solution on polygons
%
figure(3)
clf
drawsol_polygon(mesh)
view(3)
%
title([...
sprintf('\\Pi^0_k of VEM solution on polygons, k=%d',k) ...
newline...
sprintf('L^2 error: ||u_e-\\Pi^0_{k}u_h||_{0,\\Omega} = %0.2e',errors.errL2P0k/errors.normL2ue) ...
newline...
sprintf('H^1 error: |u_e-\\Pi^\\nabla_{k}u_h|_{1,\\Omega} = %0.2e',errors.errH1P0k/errors.normH1ue) ...
newline...
sprintf('nodal error: l^2 = %0.2e',errors.errl2/errors.norml2ue) ...
', max = ' sprintf('%0.2e',errors.errmax/errors.normmaxue) ...
]);
drawnow
%
if POISSONVEM.STANDALONE
    movegui('northeast')
end
%
