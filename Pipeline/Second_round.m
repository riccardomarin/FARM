%%%%%
% Code for article:
% Marin, R. and Melzi, S. and Rodolà, E. and Castellani, U., FARM: Functional Automatic Registration Method for 3D Human Bodies, CGF 2019
% Github: https://github.com/riccardomarin/FARM/
% Project Page: http://profs.scienze.univr.it/~marin/farm/index.html
%%%%%

clear all; close all; 
addpath(genpath('..\FMap'));

load('basicmodel_m_lbs.mat');
% Number of basis vectors for computing the functional map.
% Larger is usually better (more accurate results) but somewhat slower.
numEigsSrc = 30;
numEigsTar = 50;
folder='Opt1';
keep_C=0;
    
list = dir('../Results/*.obj');
for fl=1:size(list,1)
    
srcmesh = list(fl).name;
tarmesh = ['optimized_',list(fl).name(6:end)];

disp(['STARTING... ',int2str(fl),': ', srcmesh]);
if(keep_C & exist(['../Testset/Result_SR/result2_', list(fl).name(6:end-4),'.mat']))
    disp('Already Done ... Skip');
    continue;
end


if keep_C & exist(['../Cache/Secondo_Round/cache_', tarmesh(1:end-4),'.mat'])
    disp('Loading Saved Cache Tar...');
    load(['../Cache/Secondo_Round/cache_', tarmesh(1:end-4),'.mat'])
else
    tar = RFR_mesh(['../Results/',folder,'/',tarmesh]);
    landmarks2 = load('cache_smpl_base_neutro.mat','landmarks2');
    landmarks2 = landmarks2.landmarks2;
    Tar = MeshInfo(tar.v, tar.f, 200);
    save(['../Cache/Secondo_Round/cache_', tarmesh(1:end-4),'.mat'],'tar','Tar','landmarks2');
end

if exist(['../Cache/cache_', srcmesh(6:end-4),'.mat'])
    disp('Loading Saved Cache Src...');
    load(['../Results/Res1/result_', srcmesh(6:end-4),'.mat'],'C','pF_lb2','landmarks1');
    load(['../Cache/cache_', srcmesh(6:end-4),'.mat'],'Src');
else
disp('ERROR: You shouldn''t be here... skip');
continue;
end

TarLaplaceBasis = Tar.laplaceBasis; TarEigenvalues = Tar.eigenvalues;
Tar.laplaceBasis = TarLaplaceBasis(:,1:numEigsTar); Tar.eigenvalues = TarEigenvalues(1:numEigsTar);
SrcLaplaceBasis = Src.laplaceBasis; SrcEigenvalues = Src.eigenvalues;
Src.laplaceBasis = SrcLaplaceBasis(:,1:numEigsSrc); Src.eigenvalues = SrcEigenvalues(1:numEigsSrc);
Pi=sparse([1:Tar.nv], pF_lb2,1,Tar.nv,Src.nv);
C= Tar.laplaceBasis(:,1:numEigsTar)'*Tar.Ae*Pi*Src.laplaceBasis(:,1:numEigsSrc);
[C, pF_lb2] = map_boost2(Src,Tar,C,pF_lb2);
output = Tar.laplaceBasis*C*Src.laplaceBasis'*(Src.Ae*Src.X);
Joints_Target = S*output;
save(['../Results/Res2/result2_',list(fl).name(6:end-4),'.mat'],'pF_lb2','Tar','Src','landmarks1','landmarks2',...
    'C','output','Joints_Target','-append')
end

exit;
