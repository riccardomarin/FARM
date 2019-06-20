function [C_refined, matches_refined] = refine_matches(...
    part, model, part_ind, model_ind, C_init, options)

F = part_ind;
G = model_ind;

k = size(C_init,2);

W = zeros(k);
for i=1:k
    for j=1:k
        slope = 1;
        direction = [1 slope];
        direction = direction./norm(direction);
        W(i,j) = exp(-0.03*sqrt(i.^2 + j.^2))*norm(cross([direction 0], [i,j, 0]-[1 1 0]));
    end
end
W = [W; zeros(size(C_init,1)-k,k)];

d = ones(1, size(C_init,2));
D = ones(size(C_init));

if isfield(options, 'mu1')
    mu1 = options.mu1;
else
    mu1 = 1e-2; % diagonal mask
end

if isfield(options, 'mu2')
    mu2 = options.mu2;
else
    mu2 = 1e1; % orthogonality
end

for iter=1:options.refine_iters
    
    fprintf('Outer iteration: %d/%d\n', iter, options.refine_iters);
    
    A = part.evecs'*part.S*F;
    B = model.evecs'*model.S*G;
    
    manifold = euclideanfactory(size(C_init,1), size(C_init,2));
    problem = {};
    
    problem.M = manifold;
    
    problem.cost = @(C) (...
        sum(sum((C*A-B).^2).^0.5) + ...
        mu1 * norm(C.*W,'fro')^2 + ...
        mu2 * (norm(C'*C,'fro')^2 - sum(diag(C'*C).^2) + sum((diag(C'*C) - d').^2) ));
    
    problem.egrad = @(C) (...
        norm_21_gradient(C,A,B) + ...
        mu1 * 2 * C.*W.*W + ...
        mu2 * 4*(C*C'*C - C.*D ));
    
    options.maxiter = 5e3;
    options.verbosity = 1;
    C_refined = conjugategradient(problem, C_init, options);
    
%     figure,colormap(bluewhitered)
%     subplot(121),imagesc(C_init),colorbar,axis image
%     subplot(122),imagesc(C_refined),colorbar,axis image
    
    [matches_refined, ~] = flann_search(...
        C_refined*part.evecs', ...
        model.evecs', ...
        1, struct());
    
    C_init = C_refined;
    
    fps = fps_euclidean(model.VERT, 1e3, randi(model.n));
    F = sparse(matches_refined(fps), 1:length(fps), 1, part.n, length(fps));
    G = sparse(fps, 1:length(fps), 1, model.n, length(fps));
    
end

end
