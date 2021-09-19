function [WingAxesSlopes,tipPts]=findSlopesOfWingAxesManualRef(mask,symAxis,realCen,segPts,boundingBox,newTipList)

lowerSegMaskPts=[[0,size(mask,1)];[size(mask,2),size(mask,1)];segPts(4,:);realCen;segPts(1,:)];
lowerSegMask = roipoly(mask,round(lowerSegMaskPts(:,1)),round(lowerSegMaskPts(:,2)));
lowerMask=immultiply(mask,lowerSegMask);

% skeletonMask = bwmorph(mask,'skel',Inf);
% skeletonEndPtsMap = bwmorph(skeletonMask,'endpoint',Inf);
% [skeletonEndPtX,skeletonEndPtY] = find(skeletonEndPtsMap==1);
% skeletonEndPts=[skeletonEndPtY,skeletonEndPtX];
% disp('The skeleton Image has been created.');
%figure,imshow(mask); hold on;
%plot(skeletonEndPts(:,1),skeletonEndPts(:,2),'r*');

%To prevent the result from round shape wing may be biased by the lateral
%edge, the front end of bounding box are lifted a little bit.
% liftedLength=100;
% disp(['Variable [liftedLength]: ',num2str(liftedLength)]);
%Left end
%tipL=findClosestPt(skeletonEndPts,boundingBox(1:2)-[0 liftedLength]);
tipL=newTipList(1,:);
%Right end
%tipR=findClosestPt(skeletonEndPts,[boundingBox(1)+boundingBox(3),boundingBox(2)-liftedLength]);
tipR=newTipList(2,:);
%figure,imshow(mask);hold on;plot(tipL(1),tipL(2),'r*');plot(tipR(1),tipR(2),'r*');plot(segPts(2,1),segPts(2,2),'r*');
disp('Two tips of fore-wing front ends are found.');

VectorL=tipL-segPts(2,:);
foreWingLongAxisSlopeL=VectorL(2)/VectorL(1);
VectorR=tipR-segPts(3,:);
foreWingLongAxisSlopeR=VectorR(2)/VectorR(1);

nSection=10;
partMask=lowerMask;
%Left
%tarCorner=segPts(end,:);
%LeftRight='L';
disp('Start to find the inner tangent line of Left Hind-wing.');
try
    tangentLineSlopeL=find_inner_tangent_line_slope2(partMask,symAxis,segPts(end,:),'L', nSection,boundingBox);
    if isempty(tangentLineSlopeL) tangentLineSlopeL=symAxis(2)/symAxis(1);, end; 
catch
    tangentLineSlopeL=symAxis(2)/symAxis(1);
end
disp('The slope of inner tangent line of Left hind-wing is determined.');
%Right
disp('Start to find the inner tangent line of Right hind-wing.');
try
    tangentLineSlopeR=find_inner_tangent_line_slope2(partMask,symAxis,segPts(end-1,:),'R', nSection,boundingBox);
    if isempty(tangentLineSlopeR) tangentLineSlopeR=symAxis(2)/symAxis(1);, end;
catch
    tangentLineSlopeR=symAxis(2)/symAxis(1);
end
disp('The slope of inner tangent line of Right hind wing is determined.');

tipPts=[tipL; tipR];
WingAxesSlopes=[foreWingLongAxisSlopeL, tangentLineSlopeL, foreWingLongAxisSlopeR, tangentLineSlopeR];
end