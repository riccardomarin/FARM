%%%%%
% Code for article:
% Marin, R. and Melzi, S. and Rodolà, E. and Castellani, U., FARM: Functional Automatic Registration Method for 3D Human Bodies, CGF 2019
% Github: https://github.com/riccardomarin/FARM/
% Project Page: http://profs.scienze.univr.it/~marin/farm/index.html
%%%%%

clear all;
addpath(genpath('../FMap'));
list = dir('../Results/*.obj');
addpath(genpath('ARAP'));

for f=1:size(list,1)

	 name=list(f).name(6:end);
	 if(exist(['../Results/ARAP/arap_',name]))
		 continue
	 end
	 
	o_smpl = readObj(['../Results/Opt2/optimized2_', name]);
	o_target = readObj(['../Results/dato_', name]);
	o_target.n=per_vertex_normals(o_target.v,o_target.f);
	delta_t=0.01;

	a_arap=0.5;
	a_data=1;

	new=o_smpl.v;
	disp(name);
	for i=1:400
		  
		if(a_arap>0)
			[G,E] = arap_gradient(o_smpl.v,o_smpl.f,new);
			if(E>1)
				o_smpl.v=new;
				disp(E)
			end
		else
			G=zeros(6890,3);
		end

		o_new.v=new;
		o_new.f=o_smpl.f;
		o_new.n=per_vertex_normals(o_new.v,o_new.f);
		
		if (mod(i, 50) == 1)
			if (i > 350)
				targetId = knnsearch(o_target.v, new);
			else
			   targetId=myknn2(o_new,o_target);
			end
		end
		
		G2=new-o_target.v(targetId,:);
		new=new-delta_t*(a_arap*G+a_data*G2);
		if(mod(i,5)==0)

		end

		if a_arap>0.9
			a_arap=a_arap-0.005;
		end
	end

	load(['../Cache/cache_',name(1:end-4),'.mat'])
	fileid=fopen(['../Results/ARAP/arap_',name],'w');
	fprintf(fileid,'v %6.6f %6.6f %6.6f\n',new');
	fprintf(fileid,'f %d %d %d\n',o_smpl.f');
	fclose(fileid);

end