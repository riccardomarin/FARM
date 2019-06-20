function list=myknn2(shape1,shape2)
    list=zeros(1,size(shape1.v,1));
    v1=shape1.v;
    v2=[shape2.v , [1:size(shape2.v,1)]'];
    for i=1:size(shape1.v,1)
        angles = acos(shape1.n(i,:)*shape2.n(:,:)');
        mask = [(angles < 1/2*pi)'];
        if(nnz(mask))
        valid=v2(mask(:,1)==1,:);
        res=[pdist2(v1(i,:),valid(:,[1:3]))' valid(:,4)];
        [minVal minInd] = min(res(:,1));
        if(minVal > 0.1)
            mask = [(angles < 1/4*pi)'];
            valid=v2(mask(:,1)==1,:);
            res=[pdist2(v1(i,:),valid(:,[1:3]))' valid(:,4)];
            [minVal minInd] = min(res(:,1));
            if(minVal>0.005)
                mask = [(angles < 1/2*pi)'];
                valid=v2(mask(:,1)==1,:);
                res=[pdist2(v1(i,:),valid(:,[1:3]))' valid(:,4)];
                [minVal minInd] = min(res(:,1));
            end
        end
        list(i)=res(minInd,2);
        else
            l=knnsearch( shape2.v,shape1.v);
            list(i)=l(i);
        end
    end
end