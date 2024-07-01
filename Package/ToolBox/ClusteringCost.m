%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML101
% Project Title: Evolutionary Automatic Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%
function [z, out] = ClusteringCost(s, X)
    Method = 'DB';	% DB or CS
    m = s(:,1:end-1);
    
    Threshold = 0.4;
    
    a = s(:,end);
    if sum(a>=Threshold)<2
        [~, SortOrder] = sort(a,'descend');
        a(SortOrder(1:2)) = 1;
    end
    
    m = m(a>=Threshold,:);
    
    switch Method
        case 'DB'
            [z, out] = DBIndex(m, X);
            out.m = m;
            
        case 'CS'
            [z, out] = CSIndex(m, X);
            
    end
    
    out.a = double(a>=0.5);
    
end