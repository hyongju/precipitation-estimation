function [pout,index] = nearestPntDist(pnt, pos)
    pntpos = repmat(pnt,size(pos,1),1);
    postmp = pntpos - pos;
    [pout,index] = min(sqrt(postmp(:,1).^2 + postmp(:,2).^2));
%     pout = pos(index,:);
end