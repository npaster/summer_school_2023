clear all
close all

figure(1)
figure(10)
figure(11)

errors = [];

lambda=1;
mu=1;

k = 1;
%
 for N=[2 4 8]
     mesh = square2mesh_octave('s',N);
     mesh = vem_elasticity(k,lambda, mu, mesh,'test1'); 
     errors = [errors mesh.errors]; 
     ploterrors_vem_elasticity
 end


%for m={'hexagon-8x10d' 'hexagon-18x20d' 'hexagon-26x30d'}
%    mesh = square2mesh_octave(['hexagons-deformed/' char(m)]);
%    mesh = vem_elasticity(k,lambda, mu, mesh,'test1'); 
%    errors = [errors mesh.errors]; 
%    ploterrors_vem_elasticity
%end

% for m={'pegasus-4x4' 'pegasus-8x8' 'pegasus-16x16'}
%     mesh = square2mesh_octave(['pegasus/' char(m)]);
%     mesh = vem_elasticity(k,lambda, mu, mesh,'test1'); 
%     errors = [errors mesh.errors]; 
%     ploterrors_vem_elasticity
% end
