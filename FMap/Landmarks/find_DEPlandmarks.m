%%%%%
% Code for article:
% Marin, R. and Melzi, S. and Rodolà, E. and Castellani, U., FARM: Functional Automatic Registration Method for 3D Human Bodies, CGF 2019
% Github: https://github.com/riccardomarin/FARM/
% Project Page: http://profs.scienze.univr.it/~marin/farm/index.html
%%%%%

function landm = find_DEPlandmarks(Tmp, visualize)

Tmp.shape.X = Tmp.X(:,1); Tmp.shape.Y = Tmp.X(:,2); Tmp.shape.Z = Tmp.X(:,3); Tmp.shape.T = Tmp.T;

% compute and prepare DEP
Tmp = DEPinfo(Tmp,0.1,[0.04 1]);
[~, maxDEP] =  max(Tmp.DEP(:,end));
Tmp.maxDEP = dijkstra_to_all(Tmp.shape,maxDEP); 
Tmp.maxDEP = (Tmp.maxDEP - min(Tmp.maxDEP))./(max(Tmp.maxDEP) - min(Tmp.maxDEP));

% detect the head
% [~, tmphead] = max((Tmp.DEP(:,1)).*double((Tmp.DEP(:,end) < 0.8) & (Tmp.DEP(:,end) > 0.7))); % 0.85 0.4
% pHead=dijkstra_to_all(Tmp.shape,tmphead);
thr1 = 0.9;
tmphead = find( (Tmp.DEP(:,1)>thr1*max(Tmp.DEP(:,1))) & ((Tmp.DEP(:,end) < 0.85) & (Tmp.DEP(:,end) > 0.6)) ); % 0.9
while ((length(tmphead) == 0) && (thr1>0.1))
thr1 = thr1 - 0.05;
tmphead = find( (Tmp.DEP(:,1)>thr1*max(Tmp.DEP(:,1))) & ((Tmp.DEP(:,end) < 0.85) & (Tmp.DEP(:,end) > 0.6)) ); % 0.9    
end
% spectral candidates  for the head
[~,candidates_pos]=max(Tmp.laplaceBasis(:,[3:6]));
[~,candidates_neg]=min(Tmp.laplaceBasis(:,[3:6]));
% candidates = [candidates_pos'; candidates_neg'];
tmp_candidates = [candidates_pos'; candidates_neg'];
candidates = tmp_candidates( find((Tmp.DEP(tmp_candidates,end) < 0.85) & (Tmp.DEP(tmp_candidates,end) > 0.5)) );
for u = 1:length(candidates)
    d_candidates(:,u)=dijkstra_to_all(Tmp.shape,candidates(u));
end
[~, id_head] = min(min(d_candidates(tmphead,:)));

% vecchio
% [~, id_head] = min(pHead(candidates));
Tmp.head = candidates(id_head);

if visualize
    figure; colormap(jet)
    subplot(141)
    trisurf(Tmp.T,Tmp.X(:,1),Tmp.X(:,2),Tmp.X(:,3),Tmp.DEP(:,1),'SpecularStrength',0.15);
    axis equal; axis off; view([0 90]); hold on; shading interp; hold on;
    scatter3(Tmp.X(tmphead,1),Tmp.X(tmphead,2),Tmp.X(tmphead,3),40,[1 1 1],'filled');
    subplot(142)
    trisurf(Tmp.T,Tmp.X(:,1),Tmp.X(:,2),Tmp.X(:,3),Tmp.DEP(:,end),'SpecularStrength',0.15);
    axis equal; axis off; view([0 90]); hold on; shading interp; hold on;
    scatter3(Tmp.X(maxDEP,1),Tmp.X(maxDEP,2),Tmp.X(maxDEP,3),40,[0 0 0],'filled');
    subplot(143)
    trisurf(Tmp.T,Tmp.X(:,1),Tmp.X(:,2),Tmp.X(:,3),Tmp.maxDEP,'SpecularStrength',0.15);
    axis equal; axis off; view([0 90]); hold on; shading interp;
    subplot(144)
    trisurf(Tmp.T,Tmp.X(:,1),Tmp.X(:,2),Tmp.X(:,3),double((Tmp.DEP(:,end) < 0.85) & (Tmp.DEP(:,end) > 0.6)),'SpecularStrength',0.15);
    axis equal; axis off; view([0 90]); hold on; shading interp; hold on;
    scatter3(Tmp.X(Tmp.head,1),Tmp.X(Tmp.head,2),Tmp.X(Tmp.head,3),40,[0 1 1],'filled');
end
%  compute distances from the head landmarks
Tmp.dhead=dijkstra_to_all(Tmp.shape,Tmp.head);

thr = [0.90 0.3]; % 0.88 0.33 [0.85 0.37]
Tmp.dhead = (Tmp.dhead - min(Tmp.dhead))./(max(Tmp.dhead) - min(Tmp.dhead));
Tmp.ff = double( (Tmp.DEP(:,end) < thr(1)) & (Tmp.dhead > thr(2)) );
Tmp_ind = Tmp.ff.*double( (Tmp.maxDEP<0.6*max(Tmp.maxDEP)) | ( (Tmp.maxDEP<0.75*max(Tmp.maxDEP)) & (Tmp.dhead > 0.5)) ); 

if visualize
    figure; colormap(jet)
    subplot(2,4,1)
    trisurf(Tmp.T,Tmp.X(:,1),Tmp.X(:,2),Tmp.X(:,3),double(Tmp.dhead > thr(2)),'SpecularStrength',0.15);
    axis equal; axis off; view([0 90]); shading interp; hold on; 
    scatter3(Tmp.X(Tmp.head,1),Tmp.X(Tmp.head,2),Tmp.X(Tmp.head,3),40,[1 0 0],'filled');
    subplot(2,4,2)
    trisurf(Tmp.T,Tmp.X(:,1),Tmp.X(:,2),Tmp.X(:,3),double(Tmp.DEP(:,end) < thr(1)),'SpecularStrength',0.15);
    axis equal; axis off; view([0 90]); shading interp;
    subplot(2,4,3)
    trisurf(Tmp.T,Tmp.X(:,1),Tmp.X(:,2),Tmp.X(:,3),Tmp.ff,'SpecularStrength',0.15);
    axis equal; axis off; view([0 90]); shading interp;
    subplot(2,4,4)
    trisurf(Tmp.T,Tmp.X(:,1),Tmp.X(:,2),Tmp.X(:,3),Tmp_ind,'SpecularStrength',0.15);
    axis equal; axis off; view([0 90]); shading interp;
end
for i = 1:4
    if ((i == 4) && sum(Tmp_ind)==0)
        [~, Tmp.js(i)] = max(min([Tmp.dhead, Tmp.djs],[],2)); 
    else    
        d = Tmp.dhead.*Tmp_ind; d(d==0) = 100;
        [~, Tmp.js(i)] =  min(d);
    end
    Tmp.djs(:,i)=dijkstra_to_all(Tmp.shape,Tmp.js(i)); 
    Tmp.djs(:,i) = (Tmp.djs(:,i) - min(Tmp.djs(:,i)))./(max(Tmp.djs(:,i)) - min(Tmp.djs(:,i)));
    Tmp_ind((Tmp.djs(:,i) < 0.35 ) | ((Tmp.djs(:,i) < 0.55 ) & (Tmp.DEP(:,end) < 0.78 ) & (Tmp.dhead > (Tmp.dhead(Tmp.js(i)) + 0.1 ) ))) = 0; % (1/3)
    if visualize
        subplot(2,4,4+i)
        trisurf(Tmp.T,Tmp.X(:,1),Tmp.X(:,2),Tmp.X(:,3),Tmp_ind,'SpecularStrength',0.15);
        axis equal; axis off; view([0 90]); shading interp; hold on; 
        scatter3(Tmp.X(Tmp.js(i),1),Tmp.X(Tmp.js(i),2),Tmp.X(Tmp.js(i),3),40,[1 1 0],'filled');
    end
end
%% compute the extreme of each regions
[~, column] = min(Tmp.djs,[],2);
Tmp.regs = Tmp.ff.*column;
for i = 1:4
    [~, lm_tmp(i)] = max(Tmp.dhead.*double(Tmp.regs == i));
%     if ((Tmp.djs(lm_tmp(i),i) < 0.2) && (Tmp.dhead(lm_tmp(i)) < 0.5))
%         dtm = dijkstra_to_all(Tmp.shape,lm_tmp(i)); 
%         dtm = (dtm - min(dtm))./(max(dtm) - min(dtm));
%         [~, lm_tmp(i)] = max((1./dtm).*double((Tmp.djs(:,i) - 0.18) > 0).*(Tmp.dhead>thr(2)) );
%         clear dtm;
%     end
end

control=dijkstra_to_all(Tmp,Tmp.js(1));
if min(control(find(Tmp.regs == 3))) < min(control(find(Tmp.regs == 4)))
    lm = [Tmp.head lm_tmp(1) lm_tmp(3) lm_tmp(4) lm_tmp(2)];
    js = [Tmp.head Tmp.js(1) Tmp.js(3) Tmp.js(4) Tmp.js(2)];
else
    lm = [Tmp.head lm_tmp(1) lm_tmp(4) lm_tmp(3) lm_tmp(2)];
    js = [Tmp.head Tmp.js(1) Tmp.js(4) Tmp.js(3) Tmp.js(2)];
end

%%
lm2 = lm;
testa = lm(1);
mano1 = lm(2);
piede1 = lm(3);
piede2 = lm(4);
mano2 = lm(5);
lm2 = [testa, piede1, piede2, mano1, mano2];
landm=[testa, mano1, piede1, piede2, mano2];
