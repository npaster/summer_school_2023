function [ue1,ue2,epsilonU,f1,f2,g1,g2]=set_testcase_elasticity_octave(testcase,lambda,mu)

switch testcase
    case 'test1'
        ue1 = @(x,y) sin(pi*x).*sin(pi*y);
        ue2 = @(x,y) sin(pi*x).*sin(pi*y);
        t4 =@(x,y) cos(x.*pi);
        t5 =@(x,y) cos(y.*pi);
        t6 =@(x,y) sin(x.*pi);
        t7 =@(x,y) sin(y.*pi);
        t8 =@(x,y) (t4(x,y).*t7(x,y).*pi)./2.0;
        t9 =@(x,y) (t5(x,y).*t6(x,y).*pi)./2.0;
        t10=@(x,y) t8(x,y)+t9(x,y);
        epsilonU = @(x,y) reshape([t4(x,y).*t7(x,y).*pi,t10(x,y),t10(x,y),t5(x,y).*t6(x,y).*pi],[2,2]);        
        g1=@(x,y) ue1(x,y);
        g2=@(x,y) ue2(x,y);
        f1=@(x,y) -pi*pi*(-(3*mu+lambda)*sin(pi*x).*sin(pi*y)+(mu+lambda)*cos(pi*x).*cos(pi*y));
        f2=@(x,y) -pi*pi*(-(3*mu+lambda)*sin(pi*x).*sin(pi*y)+(mu+lambda)*cos(pi*x).*cos(pi*y));    
    case 'test2'
       ue1 = @(x,y) sin(pi*y*2.0)*(cos(pi*x*2.0)-1.0);
       ue2 = @(x,y) -sin(pi*x*2.0)*(cos(pi*y*2.0)-1.0);
        t4 =@(x,y) cos(x.*pi.*2.0);
        t5 =@(x,y) cos(y.*pi.*2.0);
        t6 =@(x,y) sin(x.*pi.*2.0);
        t7 =@(x,y) sin(y.*pi.*2.0);
        t8 =@(x,y) t4(x,y)-1.0;
        t9 =@(x,y) t5(x,y)-1.0;
        t10 =@(x,y) t6(x,y).*t7(x,y).*pi.*2.0;
        t11 =@(x,y) t4(x,y).*t9(x,y).*pi;
        t12 =@(x,y) t5(x,y).*t8(x,y).*pi;
        t13 =@(x,y) -t11(x,y);
        t14 =@(x,y) t12(x,y)+t13(x,y);
        epsilonU =@(x,y) reshape([-t10(x,y),t14(x,y),t14(x,y),t10(x,y)],[2,2]);
       f1 = @(x,y) 4*pi^2*cos(2*pi*x)*sin(2*pi*y)*(lambda + 2*mu) - 4*pi^2*lambda*cos(2*pi*x)*sin(2*pi*y) - mu*(4*pi^2*cos(2*pi*x)*sin(2*pi*y) - 4*pi^2*sin(2*pi*y)*(cos(2*pi*x) - 1));
       f2 = @(x,y) mu*(4*pi^2*cos(2*pi*y)*sin(2*pi*x) - 4*pi^2*sin(2*pi*x)*(cos(2*pi*y) - 1)) + 4*pi^2*lambda*cos(2*pi*y)*sin(2*pi*x) - 4*pi^2*cos(2*pi*y)*sin(2*pi*x)*(lambda + 2*mu);
    case 'test3' % we have null loading term
          ue1 = @(x,y) x^3-3*x*y^2;
          ue2 = @(x,y) y^3-3*y*x^2;
          t2 = x.^2;
          t3 =@(x,y) y.^2;
          t4 =@(x,y) x.*y.*6.0;
          t5 =@(x,y) t2(x,y).*3.0;
          t6 =@(x,y) t3(x,y).*3.0;
          t7 =@(x,y) -t4(x,y);
          epsilonU =@(x,y) reshape([t5(x,y)-t6(x,y),t7(x,y),t7(x,y),-t5(x,y)+t6(x,y)],[2,2]);
          g1=@(x,y) ue1(x,y);
          g2=@(x,y) ue2(x,y);
          f1= @(x,y) 0;
          f2= @(x,y) 0;
    case 'p1'
          ue1 = @(x,y) x;
          ue2 = @(x,y) y;
          epsilonU =@(x,y)reshape([1.0,0.0,0.0,1.0],[2,2]);
          g1=@(x,y) ue1(x,y);
          g2=@(x,y) ue2(x,y);
          f1= @(x,y) 0;
          f2= @(x,y) 0;
   case 'p2'
          ue1 = @(x,y) x.^2;
          ue2 = @(x,y) y.^2;
          epsilonU =  @(x,y) reshape([x.*2.0,0.0,0.0,y.*2.0],[2,2]);
          g1=@(x,y) ue1(x,y);
          g2=@(x,y) ue2(x,y);
          f1=- 2*lambda - 4*mu;
          f2=- 2*lambda - 4*mu;
end

end
