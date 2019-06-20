function M = load_obj(fname)

[coord,tri] = read_obj(fname);

M = {};
M.VERT = coord';
M.TRIV = tri';
M.n = size(M.VERT,1);
M.m = size(M.TRIV,1);

end
