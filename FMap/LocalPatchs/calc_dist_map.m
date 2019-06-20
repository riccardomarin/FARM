function d = calc_dist_map(M, p)

march = fastmarchmex('init', int32(M.TRIV-1), double(M.VERT(:,1)), double(M.VERT(:,2)), double(M.VERT(:,3)));

source = inf(M.n,1);
source(p) = 0;
d = fastmarchmex('march', march, double(source));

fastmarchmex('deinit', march);

end
