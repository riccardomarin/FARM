%%%%%
% Code for article:
% Marin, R. and Melzi, S. and Rodolà, E. and Castellani, U., FARM: Functional Automatic Registration Method for 3D Human Bodies, CGF 2019
% Github: https://github.com/riccardomarin/FARM/
% Project Page: http://profs.scienze.univr.it/~marin/farm/index.html
%%%%%

clear all; close all; 

addpath(genpath('..\FMap'));

% Number of basis vectors for computing the functional map.
numEigsSrc = 30;
numEigsTar = 50;

% Load SMPL_basic model (Neutral here)
load('basicmodel_m_lbs.mat');

% Load SMPL precomputed landmarks
tarmesh = './smpl_base_neutro.obj';
if exist(['cache_', tarmesh(3:end-4),'.mat'])
    load(['cache_', tarmesh(3:end-4),'.mat'])
else
    tar.v = v_template;
    tar.f = double(f)+1;
    Tar = MeshInfo(tar.v, tar.f, 200);
    landmarks2 = find_DEPlandmarks(Tar,0);%orderLandmarks(tar, Tar)';
    save(['./cache_', tarmesh(3:end-4),'.mat'],'tar','Tar','landmarks2');
end

list = dir('../Testset/*.*');
list = list(3:end);

for fl=1:size(list,1)
    disp(['STARTING... ',int2str(fl),': ', list(fl).name]);
    if(exist(['../Results/Res1/result_', list(fl).name(1:end-4),'.mat'],'file'))
        disp('Already Done ...');
    else

        if exist(['../Cache/cache_', list(fl).name(1:end-4),'.mat'],'file')
            disp('Loading Saved Cache...'); 
        else
            createcache(list(fl).name,Tar);
        end
    createresult(list(fl).name);
    end
end

disp('-------- FP2P FINISHED ---------');

exit;