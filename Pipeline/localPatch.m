%%%%%
% Code for article:
% Marin, R. and Melzi, S. and Rodolà, E. and Castellani, U., FARM: Functional Automatic Registration Method for 3D Human Bodies, CGF 2019
% Github: https://github.com/riccardomarin/FARM/
% Project Page: http://profs.scienze.univr.it/~marin/farm/index.html
%%%%%


addpath(genpath('..\FMap'));
list = dir('../Testset/*.*');
list = list(3:end);

for fm=1:size(list,1)
    name=list(fm).name(1:end-4);
    localMaps(name);
end

exit;