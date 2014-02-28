function [e,S,m] = dp_tseg_c(T,k,offset)
% implementation of dynamic programming of time series segmentation algorithm in :
%   R. Bellman. On the approximation of curves by line segments
%   using dynamic programming. Communications of the ACM,
%   4(6), 1961.
% Input:
% T : d by t time sequence with each sample dimension d and t samples;
% k : the specified number of segmentations the user want;
% Output:
% e : error term as representing k segments with k means;
% S : the segmentation result with length k-1;
% m : the k means;
% with cache
% Author: Chenyang Zhang;
% Data: Dec 18, 2013;
% Modified on 2/28/2014;
% Bug Fixed
global C
global V;

S = [];
if k>size(T,2)
    e = Inf;
    m = zeros(size(T,1),1);
    return;
end
if(k==1)
    m = mean(T,2);
    e = 0;
    for i = 1:size(T,2)
        e = e+norm(T(:,i)-m,2);
    end
    return;
end
m = [];
min_err = inf;
for j = 1:size(T,2)-1
    if(V(offset+1,offset+j,k-1)==1)
        e1 = C(offset+1,offset+j,k).e1;
        S1 = C(offset+1,offset+j,k).S1;
        m1 = C(offset+1,offset+j,k).m1;
    else
        [e1,S1,m1] = dp_tseg_c(T(:,1:j),k-1,offset);
        C(offset+1,offset+j,k).e1 = e1;
        C(offset+1,offset+j,k).S1 = S1;
        C(offset+1,offset+j,k).m1 = m1;
        V(offset+1,offset+j,k)=1;
    end
    
    if(V(offset+j+1,offset+size(T,2),1)==1)
        e2 = C(offset+j+1,offset+size(T,2),1).e2;
        S2 = C(offset+j+1,offset+size(T,2),1).S2;
        m2 = C(offset+j+1,offset+size(T,2),1).m2;
    else
        [e2,S2,m2] = dp_tseg_c(T(:,j+1:end),1,j+offset);
        C(offset+j+1,offset+size(T,2),1).e2 = e2;
        C(offset+j+1,offset+size(T,2),1).S2 = S2;
        C(offset+j+1,offset+size(T,2),1).m2 = m2;
        V(offset+j+1,offset+size(T,2),1)=1;
    end
    
    if e1+e2<min_err
        min_err = e1+e2;
        S = [S1,j,S2];
        m = [m1,m2];
    end
end
e = min_err;
end
