%
% NOTE: tris_to_keep can either be a logical array or an index array
%
function P = extract_patch(M, tris_to_keep)

if size(tris_to_keep,1)==M.n
    tris_to_keep = tris_to_keep(M.TRIV(:,1)) | tris_to_keep(M.TRIV(:,2)) | tris_to_keep(M.TRIV(:,3));
end

P = {};

P.TRIV = M.TRIV(tris_to_keep,:);
P.m = size(P.TRIV,1);

verts_to_keep = unique(P.TRIV);
P.VERT = M.VERT(verts_to_keep,:);
P.n = size(P.VERT,1);

P.fullshape_idx = verts_to_keep;

old_to_new = zeros(M.n,1);
old_to_new(verts_to_keep) = 1:length(verts_to_keep);

P.TRIV = reshape(old_to_new(P.TRIV(:)), P.m, 3);

end
