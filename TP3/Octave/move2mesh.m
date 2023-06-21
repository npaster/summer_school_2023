function mesh=move2mesh(x0,y0,mesh)

for iv=1:mesh.NV
    mesh.vertex(iv).x = mesh.vertex(iv).x + x0;
    mesh.vertex(iv).y = mesh.vertex(iv).y + y0;
end

mesh = rebuildmesh(mesh);

return
