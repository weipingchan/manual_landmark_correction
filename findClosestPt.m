function closestPt=findClosestPt(ptList,tarPt)
    %compute Euclidean distances:
    distances = sqrt(sum(bsxfun(@minus, ptList,tarPt).^2,2));
    %find the smallest distance and use that as an index into B:
    closestPt0=ptList(distances==min(distances),:);

    if size(closestPt0,1)>1
        closestPt=closestPt0(1,:);
    else
        closestPt=closestPt0;
    end
end