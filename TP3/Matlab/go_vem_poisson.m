clear all
close all

figure(1)
figure(2)
figure(3)
figure(10)
figure(11)
autoArrangeFigures(0,0,1)

errors = [];

%%% CODE STARTS HERE %%%

k = 2;

for N=[4 8 16 32 ]
     mesh = square2mesh('s',N);
     mesh = vem_poisson(k,mesh,'poisson'); 
     errors = [errors mesh.errors]; 
     ploterrors_vem_poisson
 end

% for N=[4 8 16 32]
%      mesh = square2mesh('t*',N);
%      mesh = vem_poisson(k,mesh,'poisson'); 
%      errors = [errors mesh.errors]; 
%      ploterrors_vem_poisson
%  end
 
%  k = 2;
% % %
%  for m={'hexagon-8x10d' 'hexagon-18x20d' 'hexagon-26x30d'}
%    mesh = square2mesh_octave(['hexagons-deformed/' char(m)]);
%         mesh = vem_poisson(k,mesh,'poisson'); 
%       errors = [errors mesh.errors]; 
%       ploterrors_vem_poisson
%  end

% k = 3;
% %
% for m={'pegasus-4x4' 'pegasus-8x8' 'pegasus-16x16'}
%     mesh = square2mesh(['pegasus/' char(m)]);
%       mesh = vem_poisson(k,mesh,'poisson'); 
%       errors = [errors mesh.errors]; 
%       ploterrors_vem_poisson
% end


