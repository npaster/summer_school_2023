function set_testcase(testcase)

syms x y real

switch testcase
    case 'poisson'
        %ue = log((x+1)^2+(y+1)^2)+exp(x)*sin(y);        
        %ue = x^2*y + sin(2*pi*x)*sin(2*pi*y)+2;        
        %ue = x^2+2*y^2;        
        %ue = x^3*y-x*y^3 + log((1+x)^2+(1+y)^2);        
        ue = x^3 - x*y^2 + x^2*y -x*y + x^2 -x + y -1 + sin(5*x)*sin(7*y) + log(1+x+y^3);        
        %ue = x^3*y + 3*x^4 + x^3 +2*y^4 - x^2*y^2 + x^2*y -x*y + x^2 -x + y -1;
        %ue = sin(pi*10*x)*sin(pi*10*y);
    case 'p1'
        ue = 1 - 2*x + 3*y;
    case 'p2'
        ue = 1 - 2*x + 3*y + x^2 - 2*x*y - y^2;
    case 'p3'
        ue = 1 - 2*x + 3*y + x^2 - 2*x*y - y^2 + 3*x^3 - 2*x^2*y + x*y^2 + 2*y^3;
    case 'paraboloid'
        ue = (x-1/2)^2+(y-1/2)^2;
end

uex = diff(ue,x);
uey = diff(ue,y);

f = -diff(uex,x)-diff(uey,y);
f = simplify(f);
% rhs
matlabFunction(f,'File','f','Vars',[x y]);
% exact sol
matlabFunction(ue,'File','ue','Vars',[x y]);
% first component of gradient
matlabFunction(uex,'File','uex','Vars',[x y]);
% second component of gradient
matlabFunction(uey,'File','uey','Vars',[x y]);
% boundary conditions
matlabFunction(ue,'File','g','Vars',[x y]);

%ezplot(ue,[0,1],[0,1])
%drawnow

end
