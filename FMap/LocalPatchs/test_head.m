%%%%%
% Code for article:
% Marin, R. and Melzi, S. and Rodolà, E. and Castellani, U., FARM: Functional Automatic Registration Method for 3D Human Bodies, CGF 2019
% Github: https://github.com/riccardomarin/FARM/
% Project Page: http://profs.scienze.univr.it/~marin/farm/index.html
%%%%%

function [smpl_idx, Dato_idx]= test_head(name, plotflag)

rng('default') 
rng(0)

M_full = load_obj('smpl_base_neutro.obj');

N_full = load_obj(['../Results/Dato_',name,'.obj']);
load(['../Results/Res1/result_',name,'.mat']);

dist_thresh = 0.37;

% extract heads
M = extract_patch(M_full, calc_dist_map(M_full, landmarks2(1)) < dist_thresh);
N = extract_patch(N_full, calc_dist_map(N_full, landmarks1(1)) < dist_thresh);
M_headtip = find(landmarks2(1)==M.fullshape_idx);
N_headtip = find(landmarks1(1)==N.fullshape_idx);

% center at zero
M.VERT = M.VERT - repmat(mean(M.VERT), M.n, 1);
N.VERT = N.VERT - repmat(mean(N.VERT), N.n, 1);

% establish local reference frame passing throug the head tip
[M.x,M.y,M.z] = calc_head_lrf(M, M_headtip);
[N.x,N.y,N.z] = calc_head_lrf(N, N_headtip);

% align the head verticals
M.VERT = M.VERT*[M.x M.y M.z];
[M.x,M.y,M.z] = calc_head_lrf(M, M_headtip);
N.VERT = N.VERT*[N.x N.y N.z];
[N.x,N.y,N.z] = calc_head_lrf(N, N_headtip);

n_angles = 16;
angles = linspace(0, 2*pi*(1-1/n_angles), n_angles);

N_orig = N;
solutions = cell(1,n_angles);

parfor iter=1:n_angles
    
    fprintf('outer iter %d/%d\n', iter, n_angles)
    
    angle = angles(iter);
    R = [cos(angle) -sin(angle) 0 ; sin(angle) cos(angle) 0 ; 0 0 1];
    
    N_deformable = N_orig;
    N_deformable.VERT = N_orig.VERT*R;

    N_deformable = align_icp_cpd(M, N_deformable);
    matches = knnsearch(M.VERT, N_deformable.VERT);
    
    solutions{iter}.N = N_deformable;
    solutions{iter}.matches = matches;
  
    
end

best_sol = 0;
best_distortion = Inf;

for i=1:length(solutions)
    
    matches = solutions{i}.matches;
    
    N_pts = fps_euclidean(N_orig.VERT, 100, 1);
    M_pts = matches(N_pts);
    
    DM = calc_dist_matrix(M, M_pts);
    DN = calc_dist_matrix(N_orig, N_pts);
    distortion = sum(sum(abs(DM - DN)));
    
    if distortion<best_distortion
        best_distortion = distortion;
        best_sol = i;
    end
    
    fprintf('(%d) Distortion: %.4f\n', i, distortion)
    
end

fprintf('Best solution: %.4f\n', best_distortion)
matches = solutions{best_sol}.matches;

% plot dense matches on the partial meshes
if exist('plotflag')==1
colors = create_colormap(M, M);
figure
subplot(121), colormap(colors), plot_scalar_map(M, 1:M.n); shading flat; freeze_colors
subplot(122), colormap(colors(matches,:)), plot_scalar_map(N_orig, 1:N_orig.n); shading flat
end

N_pts = fps_euclidean(N_orig.VERT, 100, 1);
M_pts = matches(N_pts);

% plot point matches on the partial meshes
M2 = M;
M2.VERT = M.VERT(M_pts,:);
% colors = create_colormap(M2,M2);
% figure, colormap([1 1 1])
% subplot(121), plot_scalar_map(M, ones(M.n,1)); light; lighting phong
% hold on, plot_cloud_color(M.VERT(M_pts,:), colors, 20)
% subplot(122), plot_scalar_map(N_orig, ones(N_orig.n,1)); light; lighting phong
% hold on, plot_cloud_color(N_orig.VERT(N_pts,:), colors, 20)

% plot point matches on the full meshes
Dato_idx = N_orig.fullshape_idx;
smpl_idx = M.fullshape_idx(matches);
% figure, colormap([1 1 1])
% subplot(121), plot_scalar_map(M_full, ones(M_full.n,1)); light; lighting phong
% hold on, plot_cloud_color(M_full.VERT(M.fullshape_idx(M_pts),:), colors, 20)
% subplot(122), plot_scalar_map(N_full, ones(N_full.n,1)); light; lighting phong
% hold on, plot_cloud_color(N_full.VERT(N_orig.fullshape_idx(N_pts),:), colors, 20)
end