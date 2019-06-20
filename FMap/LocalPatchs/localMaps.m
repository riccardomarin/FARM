%%%%%
% Code for article:
% Marin, R. and Melzi, S. and Rodolà, E. and Castellani, U., FARM: Functional Automatic Registration Method for 3D Human Bodies, CGF 2019
% Github: https://github.com/riccardomarin/FARM/
% Project Page: http://profs.scienze.univr.it/~marin/farm/index.html
%%%%%

function localMaps(name)

    load(['../Results/Res1/result_', name,'.mat'],'landmarks1','landmarks2','Tar','Src');
    
    [smpl_idx, dato_idx]=test_head(name);
    %Calcolo pesi patch faccia
    w_head_s=zeros(size(Tar.X,1),1);
    w_head_s(smpl_idx)=1;
    dist=dijkstra_to_all(Tar,landmarks2(1));
    dist_head=dist(smpl_idx);
    max_d=max(dist_head);
    min_d=min(dist_head);
    funct=(1-(dist_head-min_d)/(max_d-min_d));
    w_head_s(smpl_idx)=(funct*1).^(2/3);

    w_head_t=zeros(size(Src.X,1),1);
    w_head_t(dato_idx)=1;
    dist=dijkstra_to_all(Src,landmarks1(1));
    dist_head=dist(dato_idx);
    max_d=max(dist_head);
    min_d=min(dist_head);
    funct=(1-(dist_head-min_d)/(max_d-min_d));
    w_head_t(dato_idx)=(funct*1);
    
    if exist('plotflag')==1          
        [smpl_idx_r, dato_idx_r]=test_hand(name,2,1);
        [smpl_idx_l, dato_idx_l]=test_hand(name,5,1);
    else
        [smpl_idx_r, dato_idx_r]=test_hand(name,2);
        [smpl_idx_l, dato_idx_l]=test_hand(name,5);
    end
    
save(['../Results/Res2/result2_',name,'.mat'],'smpl_idx','dato_idx','w_head_s','w_head_t','smpl_idx_l','smpl_idx_r','dato_idx_l','dato_idx_r')
end

