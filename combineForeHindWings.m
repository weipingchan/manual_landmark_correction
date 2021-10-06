function shapeImg=combineForeHindWings(tformWingImg,tformPtData)
shapeImg=cell(2,2);
for sideInd=1:2
obID=(sideInd-1)*2+1;
forImg=tformWingImg{obID};
hindImg=tformWingImg{obID+1};
forData=tformPtData{obID};
hindData=tformPtData{obID+1};

if sideInd==1 %Left wings
    forRefA=forData(2,:); %The corner between forewing and body
    forRefB=forData(end-2,:); %The intersection between body-wing line and the centroid-fore-hindwing line
    forRefDist=pdist([forRefA;forRefB],'euclidean');
    forRef=forRefA+[0 forRefDist];
    hindRef=hindData(end-2,:);
    hindConjo=hindData(6,:); %The corner between hindwing and body
elseif sideInd==2 %Right wings
    forRefA=forData(3,:); %The corner between forewing and body
    forRefB=forData(end,:); %The intersection between body-wing line and the centroid-fore-hindwing line
    forRefDist=pdist([forRefA;forRefB],'euclidean');
    forRef=forRefA+[0 forRefDist];
    hindRef=hindData(end,:);
    hindConjo=hindData(5,:);  %The corner between hindwing and body
end

moveDis=forRef-hindRef;
if moveDis(1)<0 || moveDis(2)<0 %If the moving direction is negative, keeping the full size will change the coordination of the target object
    mvhindImg = imtranslate(hindImg,moveDis);
    disp('Full output view is NOT applied since variable [moveDis] contains negative value.');
else
    mvhindImg = imtranslate(hindImg,moveDis,'OutputView','full');
end
comImg = imfuse(forImg,mvhindImg);
hindConj=hindConjo+moveDis;
comInfo=[forRefA ; forRef ; hindRef ; hindConj];
shapeImg{sideInd,1}=comImg;
shapeImg{sideInd,2}=comInfo;
end
end