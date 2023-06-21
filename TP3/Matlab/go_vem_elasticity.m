clear all
close all

figure(1)
figure(10)
figure(11)
% figure(12)
autoArrangeFigures(0,0,1)

errors = [];

lambda=1;
mu=1;

k = 2;
%
for N=[2 4 8]
    mesh = square2mesh('s',N);
    mesh = vem_elasticity(k,lambda, mu, mesh,'test1'); 
    errors = [errors mesh.errors]; 
    ploterrors_vem_elasticity
end


% for m={'hexagon-8x10d' 'hexagon-18x20d' 'hexagon-26x30d'}
%     mesh = square2mesh(['hexagons-deformed/' char(m)]);
%     mesh = vem_elasticity(k,lambda, mu, mesh,'test1'); 
%     errors = [errors mesh.errors]; 
%     ploterrors_vem_elasticity
% end

% for m={'pegasus-4x4' 'pegasus-8x8' 'pegasus-16x16'}
%     mesh = square2mesh(['pegasus/' char(m)]);
%     mesh = vem_elasticity(k,lambda, mu, mesh,'test1'); 
%     errors = [errors mesh.errors]; 
%     ploterrors_vem_elasticity
% end


% % k = 3;
% %


% mesh = vem_elasticity(k,mesh,'goniometricFunction'); 
%     errors = [errors mesh.errors]; ploterrors_elasticity
% end

% k = 2;
% %

%      mesh = vem_elasticity(k,mesh,'goniometricFunction'); 
% %    mesh = ezvemElasticity(k,mesh,'patchTest3'); 
%         errors = [errors mesh.errors]; 
%      ploterrors_elasticity
 %end

