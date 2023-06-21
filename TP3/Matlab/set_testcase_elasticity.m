function set_testcase_elasticity(testcase)

syms x y lambda mu real

switch testcase
    case 'test1'
          ue1 = sin(pi*x).*sin(pi*y);
          ue2 = sin(pi*x).*sin(pi*y);
          f1=-pi*pi*(-(3*mu+lambda)*sin(pi*x).*sin(pi*y)+(mu+lambda)*cos(pi*x).*cos(pi*y));
          f2=-pi*pi*(-(3*mu+lambda)*sin(pi*x).*sin(pi*y)+(mu+lambda)*cos(pi*x).*cos(pi*y));
    case 'test2'
          ue1 = sin(pi*y*2.0)*(cos(pi*x*2.0)-1.0);
          ue2 = -sin(pi*x*2.0)*(cos(pi*y*2.0)-1.0);
          f1 = 4*pi^2*cos(2*pi*x)*sin(2*pi*y)*(lambda + 2*mu) - 4*pi^2*lambda*cos(2*pi*x)*sin(2*pi*y) - mu*(4*pi^2*cos(2*pi*x)*sin(2*pi*y) - 4*pi^2*sin(2*pi*y)*(cos(2*pi*x) - 1));
          f2 = mu*(4*pi^2*cos(2*pi*y)*sin(2*pi*x) - 4*pi^2*sin(2*pi*x)*(cos(2*pi*y) - 1)) + 4*pi^2*lambda*cos(2*pi*y)*sin(2*pi*x) - 4*pi^2*cos(2*pi*y)*sin(2*pi*x)*(lambda + 2*mu);
    case 'test3' % we have null loading term
          ue1 = x^3-3*x*y^2;
          ue2 = y^3-3*y*x^2;
          f1=0;
          f2=0;
    case 'p2' % we have null loading term
          ue1 = x.^2;
          ue2 = y.^2;
          f1=- 2*lambda - 4*mu;
          f2=- 2*lambda - 4*mu;
    case 'p1'
          ue1 = x;
          ue2 = y;
          f1=0;
          f2=0;
end

ue1x = diff(ue1,x);
ue1y = diff(ue1,y);
ue2x = diff(ue2,x);
ue2y = diff(ue2,y);
epsilonU =[ue1x (ue1y+ue2x)/2; (ue1y+ue2x)/2 ue2y];

% loading term
matlabFunction(f1,'File','f1','Vars',[x y lambda mu]);
matlabFunction(f2,'File','f2','Vars',[x y lambda mu]);

% strain 
matlabFunction(epsilonU,'File','epsilonU','Vars',[x y]);

% exact sol
matlabFunction(ue1,'File','ue1','Vars',[x y]);
matlabFunction(ue2,'File','ue2','Vars',[x y]);

% Dirichlet boundary condition
matlabFunction(ue1,'File','g1','Vars',[x y]);
matlabFunction(ue2,'File','g2','Vars',[x y]);

end
