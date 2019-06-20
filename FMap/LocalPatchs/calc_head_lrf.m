function [x,y,z] = calc_head_lrf(M, headtip_idx)

z = M.VERT(headtip_idx,:)';
z = z./norm(z);
y = rand(3,1);
y = y./norm(y);
y = y - (y'*z)*z;
y = y./norm(y);
x = cross(y,z);

end
