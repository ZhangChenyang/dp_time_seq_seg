function [e,S,m] = dp_tseg_c(T,k)
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
global C
global V;

S = [];
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
    if(V(j,k))==1
        e1 = C(j,k,1).e;
        e2 = C(j,k,2).e;
        S1 = C(j,k,1).S;
        S2 = C(j,k,2).S;
        m1 = C(j,k,1).m;
        m2 = C(j,k,2).m;       
    else
        [e1,S1,m1] = dp_tseg_c(T(:,1:j),k-1);
        [e2,S2,m2] = dp_tseg_c(T(:,j+1:end),1);
        V(j,k)=1;
        C(j,k,1).e = e1;
        C(j,k,2).e = e2;
        C(j,k,1).S = S1;
        C(j,k,2).S = S2;
        C(j,k,1).m = m1;
        C(j,k,2).m = m2;  
    end
    
    if e1+e2<min_err
        min_err = e1+e2;
        S = [S1,j,S2];
        m = [m1,m2];
    end
end
e = min_err;
end
