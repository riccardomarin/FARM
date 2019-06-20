function [Fs] = compute_indicator_functions(Ms,matches,alpha)
    ncorr = size(matches,2);
    parfor m = 1:numel(Ms)
        M = Ms{m};
         n = size(M.VERT,1); 
        F = [];
        Mdist = fastmarchmex('init', int32(M.TRIV-1), double(M.VERT(:,1)), double(M.VERT(:,2)), double(M.VERT(:,3)));
        for i=1:ncorr
            if(matches(m,i)<=0)
                F(:,i) = zeros(n,1);
                continue;
            end
            source = inf(n, 1);
            source(matches(m,i)) = 0;
            dM = exp(-0.5*alpha*fastmarchmex('march', Mdist, double(source)));
            F(:,i)=dM;   
        end
         Fs{m}=F;
%          Fs{m}=[Fs{m} M.HKS];
%          Fs{m}=[Fs{m} M.WKS];
        fastmarchmex('deinit', Mdist);
    end
end

%         for i=1:ncorr
%                 F(:,i) = zeros(n,1);
%                 if(matches(m,i)>0)
%                  F(matches(m,i),i)=1;
%                 end
%                 
%         end
%         Fs{m}=F;