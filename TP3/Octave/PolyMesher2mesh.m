function mesh=PolyMesher2mesh(Node,Element,Supp,Load,P)
%
% PolyMesher2mesh(Node,Element,Supp,Load,P)
%
% takes the output of PolyMesher and makes a mesh structure
%
NP = length(Element);
%
lipmax = 0;
%
for ip=1:NP
    %
    lip = length(Element{ip});
    if lip>lipmax
        lipmax = lip;
    end
end
%
vertices = zeros(NP,lipmax);
%
for ip=1:NP 
    lip = length(Element{ip});    
    vertices(ip,1:lip) = Element{ip};
    %
end
%
xv = Node(:,1);
yv = Node(:,2);
%
mesh = conn2mesh(xv,yv,vertices);
%
end
