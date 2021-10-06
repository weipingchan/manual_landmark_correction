function [allComImg, allinfo]=integrate4WingsInto1Panel(shapeImg,tformWingImg,refPts)
%Forewing mix
disp('Start to combine forewings.');
Lpar=shapeImg{1,2};
Rpar=shapeImg{2,2};
bodyW=pdist([refPts(end-2,:);refPts(end,:)],'euclidean')+200; %This value determines the distance between left and right wings
Lref=Lpar(2,:)+[bodyW 0];
Rref=Rpar(2,:);
disp('The reference image is left forewing.');

LforImg=tformWingImg{1};
RforImg=tformWingImg{3};

horMvDis=Lref-Rref;
if horMvDis(1)<0 || horMvDis(2)<0 %If the moving direction is negative, keeping the full size will change the coordination of the target object
    mvRImg = imtranslate(RforImg,horMvDis);
    disp('Full output view is NOT applied since variable [horMvDis] contains negative value.');
else
    mvRImg = imtranslate(RforImg,horMvDis,'FillValues',0,'OutputView','full');
end
disp('Right forewing has been re-registered.');
FcomImg = imfuse(LforImg,mvRImg,'blend');
disp('The composition of forewings is finished.');

%Hindwing mix
disp('Start to combine hindwings.');
LhindImg=tformWingImg{2};
RhindImg=tformWingImg{4};

LverMvDis=Lpar(2,:)-(Lpar(3,:));
if LverMvDis(1)<0 || LverMvDis(2)<0
    LmvRImg = imtranslate(LhindImg,LverMvDis);
    disp('Full output view is NOT applied since variable [LverMvDis] contains negative value.');
else
    LmvRImg = imtranslate(LhindImg,LverMvDis,'OutputView','full');
end
disp('Left hindwing has been re-registered.');

RverMvDis=Rpar(2,:)-Rpar(3,:);
if RverMvDis(1)<0 || RverMvDis(2)<0
    RmvRImg0 = imtranslate(RhindImg,RverMvDis);
    disp('Full output view is NOT applied since variable [RverMvDis] contains negative value.');
else
    RmvRImg0 = imtranslate(RhindImg,RverMvDis,'OutputView','full');
end
if horMvDis(1)<0 || horMvDis(2)<0
    RmvRImg = imtranslate(RmvRImg0,horMvDis);
    disp('Full output view is NOT applied since variable [horMvDis] contains negative value.');
else
    RmvRImg = imtranslate(RmvRImg0,horMvDis,'OutputView','full');
end
disp('Right hindwing has been re-registered.');

HcomImg = imfuse(LmvRImg,RmvRImg,'blend');
disp('The composition of hindwings is finished.');
%Combine both images to derive the final one
allComImg = imfuse(FcomImg,HcomImg);
disp('The composition of all wings is finished.');

%derive new Info
allinfo=[Lpar(1,:) ; Lpar(4,:) ; Rpar(1,:)+horMvDis; Rpar(4,:)+horMvDis];
disp('The composition of all parameters has been generated.');
end