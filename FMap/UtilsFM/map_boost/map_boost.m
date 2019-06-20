%%%%%
% Code for article:
% Marin, R. and Melzi, S. and Rodolà, E. and Castellani, U., FARM: Functional Automatic Registration Method for 3D Human Bodies, CGF 2019
% Github: https://github.com/riccardomarin/FARM/
% Project Page: http://profs.scienze.univr.it/~marin/farm/index.html
%%%%%

function [C_out,matches_out] = map_boost(Src,Tar, C,pF_lb2,plotflag)

    M.evecs = Tar.laplaceBasis;
    M.S = Tar.Ae;
    M.n=Tar.nv;
    M.VERT=Tar.X;
    M.TRIV=Tar.T;

    N.evecs = Src.laplaceBasis;
    N.S = Src.Ae;
    N.n = Src.nv;
    N.VERT=Src.X;
    N.TRIV=Src.T;

    C_in = C;
    matches_in = pF_lb2;

    params = struct;
    params.fps_rand = true;
    params.variance = 0.1; % the larger, the more concentrated
    params.mu2      = 1.0;
    params.in_iter  = 1e4;

    disp('First p2p boost...');
    params.fps      = 200;
    params.variance = 0.1;
    [C_out, matches_out] = filter_21(M, N, C_in, matches_in, params);

    disp('Second p2p boost...');
    params.fps      = 300;
    params.variance = 0.5;
    [C_out, matches_out] = filter_21(M, N, C_out, matches_out, params);

    disp('Third p2p boost...');
    params.fps      = 400;
    params.variance = 1.0;
    [C_out, matches_out] = filter_21(M, N, C_out, matches_out, params);



    if exist('plotflag')==1
        figure, colormap(bluewhitered)
        subplot(121), imagesc(C_in), axis image
        subplot(122), imagesc(C_out), axis image

        colors = create_colormap(N,N);
        figure
        subplot(131), colormap(colors), plot_scalar_map(N, 1:N.n), shading flat; freeze_colors; axis off
        subplot(132), colormap(colors(matches_in,:)), plot_scalar_map(M, 1:M.n), shading flat; view([0 90]); freeze_colors; axis off
        subplot(133), colormap(colors(matches_out,:)), plot_scalar_map(M, 1:M.n), shading flat; view([0 90]); axis off
    end
end

