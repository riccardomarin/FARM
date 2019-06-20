%%%%%
% Code for article:
% Marin, R. and Melzi, S. and Rodolà, E. and Castellani, U., FARM: Functional Automatic Registration Method for 3D Human Bodies, CGF 2019
% Github: https://github.com/riccardomarin/FARM/
% Project Page: http://profs.scienze.univr.it/~marin/farm/index.html
%%%%%

function createcache(srcmesh,Tar)

    src = RFR_mesh(['../Testset/',srcmesh],'lqsim',6890);

    tmp = MeshInfo(src.v, src.f, 1);
    factorOK = sqrt(sum(Tar.area)/sum(tmp.area));
    src.v = src.v.*factorOK; 
    Src = MeshInfo(src.v, src.f, 200);
    landmarks1 = find_DEPlandmarks(Src,0);
    save(['../Cache/cache_', srcmesh(1:end-4)],'Src','src','landmarks1','factorOK');

end