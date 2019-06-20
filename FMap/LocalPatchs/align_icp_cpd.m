function N_out = align_icp_cpd(M, N)

opt = struct;
opt.verbose = false;
opt.viz=0;
opt.tol=1e-8;

opt.method='rigid'; % use rigid registration
opt.outliers=0.6;   % use 0.6 noise weight
opt.normalize=1;    % normalize to unit variance and zero mean before registering (default)
opt.scale=1;        % estimate global scaling too (default)
opt.rot=1;          % estimate strictly rotational matrix (default)
opt.corresp=1;      % compute correspondence vector at the end of registration (not being estimated by default)
opt.max_it=200;     % max number of iterations

Transform = cpd_register(M.VERT, N.VERT, opt);

% figure, colormap([1 0 0; 1 1 0])
% subplot(131)
% plot_scalar_map(M, ones(M.n,1)), light; hold on; alpha(0.5)
% plot_scalar_map(N, 2*ones(N.n,1)), hold on
% view([-215 25]); view([144 10])
% title('Before')

N_out = N;
N_out.VERT = Transform.Y;

% subplot(132)
% plot_scalar_map(M, ones(M.n,1)), light; hold on; alpha(0.5)
% plot_scalar_map(N_out, 2*ones(N_out.n,1)), hold on
% view([-215 25]); view([144 10])
% title('After rigid')

opt.method='nonrigid'; % use rigid registration
opt.beta=4;         % the width of Gaussian kernel (smoothness)
opt.lambda=10;      % regularization weight
opt.outliers=0.6;   % use 0.6 noise weight

Transform = cpd_register(M.VERT, N_out.VERT, opt);
N_out.VERT = Transform.Y;

% subplot(133)
% plot_scalar_map(M, ones(M.n,1)), light; hold on; alpha(0.5)
% plot_scalar_map(N_out, 2*ones(N_out.n,1)), hold on
% view([-215 25]); view([144 10])
% title('After nonrigid')

end
