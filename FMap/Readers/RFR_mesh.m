% method:   (lqsim, N_Vertex)
%           (reduce, N_Faces)

function o = RFR_mesh(name,method, n)

    o=UniversalReader(name);

    if nargin > 1
        o_old=o;
        old_n= per_vertex_normals(o.v,o.f);
        [o.v, o.f] = meshfix(o.v,o.f);
        if strcmp(method, 'lqsim')
            % Remesh LQSIM
            params = struct;
            params.vertices=n;
            params.verbose=0;
            surface.VERT=o.v;
            surface.TRIV=o.f;
            [surface] = remesh(surface,params);
            o.v=surface.VERT;
            o.f=surface.TRIV;
        end

        if strcmp(method, 'reduce')
            % Remesh matlab
            [o.f,o.v] = reducepatch(o.f,o.v,n);
        end
        
        [o.v, o.f] = meshfix(o.v,o.f);
        nn= per_vertex_normals(o.v,o.f);
        map=knnsearch(o_old.v,o.v);
        a=sum(sign(dot(nn(1:1000,:),old_n(map(1:1000),:),2)));
        if(a<0)
            disp('Warning: Meshfix flipped faces orientation. Fixed.');
            o.f=[o.f(:,3), o.f(:,2), o.f(:,1)];
        end
    end
end