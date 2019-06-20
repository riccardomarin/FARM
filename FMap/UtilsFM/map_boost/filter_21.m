function [C_out, matches_out] = filter_21(M, N, C_in, matches_in, params)

if params.fps_rand
    seed = randi(M.n);
else
    seed = 1;
end

fps = fps_euclidean(M.VERT, params.fps, seed);
FG = compute_indicator_functions({M,N}, [fps matches_in(fps)]', params.variance);
F = FG{1};
G = FG{2};

A = N.evecs'*N.S*G;
B = M.evecs'*M.S*F;

d = ones(1, size(C_in,2));
D = ones(size(C_in));

manifold = euclideanfactory(size(C_in,1), size(C_in,2));
problem = {};

problem.M = manifold;

problem.cost = @(C) (...
    sum(sum((C*A-B).^2).^0.5) + ...
    params.mu2 * (norm(C'*C,'fro')^2 - sum(diag(C'*C).^2) + sum((diag(C'*C) - d').^2) ));

problem.egrad = @(C) (...
    norm_21_gradient(C,A,B) + ...
    params.mu2 * 4*(C*C'*C - C.*D ));

% figure, checkgradient(problem)

options = struct;
options.maxiter = params.in_iter;
% options.tolgradnorm = 1e-06;
% options.minstepsize = 1e-06;
options.verbosity = 0;

[C_out, ~, ~, ~] = conjugategradient(problem, C_in, options);
[matches_out, ~] = knnsearch(N.S*N.evecs*C_out', M.S*M.evecs);

end
