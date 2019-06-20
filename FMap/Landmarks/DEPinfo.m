function [S] = DEPinfo(S,c,deltas)
S.D = zeros(S.nv,S.nv);
for j = 1:S.nv
    S.D(:,j) = sqrt(sum( bsxfun(@times, ((1./S.eigenvalues(2:end).^2)'), (bsxfun(@minus,S.laplaceBasis(j,2:end),S.laplaceBasis(:,2:end)).^2)) ,2));    
end
S.D  = ( S.D - min(S.D(:)) )./(max(S.D(:)) - min(S.D(:)) );

n_scale = length(deltas);
S.points = 1:S.nv;
S.dists = S.D;

% compute area elements
S.Omega = abs(diag(sum(full(S.Ae),2)));
S.f0 = ones(S.nv,1);


% compute DEP descriptor
% tic
% disp('compute DEP descriptor')
tmpDEP = zeros(S.nv,n_scale);
for q = 1 : n_scale 
    S.D_norm = S.dists;
    S.D_norm(S.dists > deltas(q)) = 1; % 0 non inv, 1 inv
    S.D_norm = ones(size(S.D_norm)) - S.D_norm; 
 
%     tmpDEP(:,q) = computeDEPdesc(S.D_norm*S.Omega,c,S.f0);
    evolutionOP = S.D_norm*S.Omega;
    EyE_matrix = eye( size( evolutionOP,1 ));
    r  = (c/eigs(evolutionOP,1,'lm')); 
    Delta1 = EyE_matrix - ( (r) *  evolutionOP );
    tmpDEP(:,q) = (Delta1\S.f0);

    tmpDEP(:,q) = (tmpDEP(:,q) - min(tmpDEP(:,q)))./max(tmpDEP(:,q) - min(tmpDEP(:,q)));
end
% toc

S.DEP = tmpDEP; 
