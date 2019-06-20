function o = UniversalReader(fname)

[filepath,name,ext] = fileparts(fname);
if ext=='.mat'
    load(fname);
    try
    o.v=[surface.X, surface.Y, surface.Z];
    o.f=surface.TRIV;
    catch
       o.v=[shape.X, shape.Y, shape.Z];
       o.f=shape.TRIV;
    end

end

if ext=='.obj'
    [o.v, o.f]=readOBJ(fname);
end  
if ext == '.off'
    [X T]=readOFF(fname);
    o.v=X;
    o.f=T;
end

if ext == '.ply'
    [X T]=readPLY(fname);
    o.v=X;
    o.f=T;
end