function drawsol_edge(mesh,ue)
%
% number of samples per edge
%
NSe = 10;
%
k = mesh.POISSONVEM.K;
uh = mesh.uh;
%
[tgl,~] = mylglnodes(k);
%
% VEM solution at boundary nodes of each edge
%
Xhn = zeros(k+1,1);
Yhn = zeros(k+1,1);
Zhn = zeros(k+1,1);
%
% exact solution on sampled nodes of each edge
%
Xes = zeros(NSe,1);
Yes = zeros(NSe,1);
Zes = zeros(NSe,1);
%
% VEM solution on sampled nodes of each edge
%
Xhs = zeros(NSe,1);
Yhs = zeros(NSe,1);
Zhs = zeros(NSe,1);
%
for ie=1:mesh.NE
    %
    v1 = mesh.edge(ie).v1;
    v2 = mesh.edge(ie).v2;
    %
    x1 = mesh.vertex(v1).x;
    y1 = mesh.vertex(v1).y;
    %
    x2 = mesh.vertex(v2).x;
    y2 = mesh.vertex(v2).y;
    %
    Xhn(1) = x1;
    Yhn(1) = y1;
    Zhn(1) = uh(v1);
    %
    Xhn(k+1) = x2;
    Yhn(k+1) = y2;
    Zhn(k+1) = uh(v2);
    %
    for ik=1:k-1
        %
        xik = (x2-x1)/2*tgl(ik+1)+(x2+x1)/2;
        yik = (y2-y1)/2*tgl(ik+1)+(y2+y1)/2;
        %
        Xhn(ik+1) = xik;
        Yhn(ik+1) = yik;
        Zhn(ik+1) = uh(mesh.NV+(ie-1)*(k-1)+ik);
        %
    end
    %
    ss = linspace(-1,1,NSe);
    %
    for is=1:NSe
        %
        xis = (x2-x1)/2*ss(is)+(x2+x1)/2;
        yis = (y2-y1)/2*ss(is)+(y2+y1)/2;
        %
        Xes(is) = xis;
        Yes(is) = yis;
        Zes(is) = ue(xis,yis);
        %
    end
    %
    Zhs = polyval(polyfit(tgl,Zhn,k),ss);
    Xhs = (x2-x1)/2*ss+(x2+x1)/2;
    Yhs = (y2-y1)/2*ss+(y2+y1)/2;
    %
    %line(Xh,Yh,Zh,'Color','k','LineWidth',1)
    line(Xhs,Yhs,Zhs,'Color','k')
    %line(Xhn,Yhn,Zhn,'Color','k','LineWidth',1)
    %line(Xe,Ye,Ze,'Color','r','LineWidth',1)
    %
end

