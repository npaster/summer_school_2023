function [ue,uex,uey,f,g]=set_testcase_octave(testcase)

switch testcase
    case 'poisson'
        f=@(x,y) (1.0./(y.^3+x+1.0).^2)-x.*4.0-y.*2.0-(y.*6.0)./(y.^3+x+1.0)+(1.0./(y.^3+x+1.0).^2).*y.^4.*9.0+sin(x.*5.0).*sin(y.*7.0).*7.4e+1-2.0;
        ue=@(x,y) x.^2-x+y+log(x+y.^3+1.0)+x.^2.*y-x.*y-x.*y.^2+x.^3+sin(x.*5.0).*sin(y.*7.0)-1.0;
        g=@(x,y) ue(x,y);
        uex=@(x,y) x.*2.0-y+1.0./(x+y.^3+1.0)+x.*2.0.*y-y.*y+x.^2.*3.0+cos(x.*5.0).*sin(y.*7.0).*5.0-1.0;
        uey=@(x,y) -x+(y.^2.*3.0)./(x+y.^3+1.0)-x.*y.*2.0+x.^2+cos(y.*7.0).*sin(x.*5.0).*7.0+1.0;
    case 'p1'
        ue=@(x,y) x.*-2.0+y.*3.0+1.0;
        g =@(x,y) ue(x,y);
        f =@(x,y) 0;
        uex =@(x,y) -2;
        uey =@(x,y) 3;
    case 'p2'
        ue =@(x,y) -2.*x+y.*3.0+t2.*y+x.^2-y.^2+1.0;
        g =@(x,y) ue(x,y);
        f =@(x,y) 0;
        uex =@(x,y) x.*2.0-y.*2.0-2.0;
        uey =@(x,y) x.*-2.0-y.*2.0+3.0;
    case 'p3'
        ue = x.^2-y.^2-x.*2.0+y.*3.0+t3.*x-x.^2.*y.*2.0-x.*y.*2.0+x.^3.*3.0+y.^3.*2.0+1.0;
        g =@(x,y) ue(x,y);
        f =@(x,y) x.*-2.0e+1-y.*8.0;
        uey =@(x,y) -2.*x-y.*2.0+t2.*x+x.*y.*2.0+y.^2.*6.0+3.0;
        uex =@(x,y) x.*2.0-y.*2.0-x.*y.*4.0+x.^2.*9.0+y.^2-2.0;       
end
end