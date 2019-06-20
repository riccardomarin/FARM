function D = dijkstra_to_all(shape, sources)

if ~isfield(shape, 'surface')
    if size(shape.X,2)==3
    shape.surface.X = shape.X(:,1);
    shape.surface.Y = shape.X(:,2);
    shape.surface.Z = shape.X(:,3);
    shape.surface.TRIV = shape.T;
    else
    shape.surface.X = shape.X;
    shape.surface.Y = shape.Y;
    shape.surface.Z = shape.Z;
    shape.surface.TRIV = shape.T;    
    end
end
    S = shape.surface;
    
    if(size(sources,2) > size(sources,1))
        sources = sources';
    end
    
    if(size(sources,2) > 1)
        error('sources must be stored in a vector');
    end
    
    D = comp_geodesics_to_all(double(S.X), double(S.Y), double(S.Z), ...
                              double(S.TRIV'), sources, 1);
end