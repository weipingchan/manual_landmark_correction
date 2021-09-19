function [bolbMorph,antFig]=oneAntennaeMeasure4(mask, forkPt00, oneAntenna,antenna1,It0,Dd0,realCen,scalelen)
     try %prevent potential error and interruption
        It = bwmorph(antenna1,'thin','inf');
        [Bbx,Bby] =  find(bwmorph(It,'branchpoints'));
        %select branch point and convert to index
        dist = sqrt(sum(bsxfun(@minus, [Bbx,Bby], flip(forkPt00)).^2,2));
        [~,minLoc]=min(abs(dist-mean(dist(dist<mean(dist)+std(dist)*1/2,:)))); %find the one closest to the mean of those points near forkPt
        Bbi=sub2ind(size(antenna1),Bbx(minLoc),Bby(minLoc)); %The forkPoint we want

        [i,j] = find(immultiply(bwmorph(It,'endpoints'),oneAntenna)); %Calculation with refined end points
        Dd = bwdistgeodesic(It,Bbi,'quasi');
        distanceAnt=[];
        for n = 1:numel(i)
            distanceAnt=[distanceAnt; [i(n),j(n),Dd(i(n),j(n))]];
        end
        %%
        if size(distanceAnt,1)>1
            %Cluster points long distances from the fork
            Dist = sqrt(sum(bsxfun(@minus, distanceAnt(:,1:2), flip(forkPt00)).^2,2));
            antennaTipDist0= sortrows(distanceAnt(Dist>mean(Dist)+std(Dist)*1/2,:),2); %first one is the left one
            antennaTipDistA0=sortrows(antennaTipDist0(antennaTipDist0(:,1)<mean(antennaTipDist0(:,1)),:),3);
            antennaTipDistB0=sortrows(antennaTipDist0(antennaTipDist0(:,1)>=mean(antennaTipDist0(:,1)),:),3);

            if ~isempty(antennaTipDistA0)
                antennaTipDistA=antennaTipDistA0(end,:);
            else
                antennaTipDistA=antennaTipDistA0;
            end
            if ~isempty(antennaTipDistB0)
                antennaTipDistB=antennaTipDistB0(end,:);
            else
                antennaTipDistB=antennaTipDistB0;
            end
            antennaTipDist=[antennaTipDistA;antennaTipDistB];
        else
            antennaTipDist=distanceAnt;
        end
    catch
        antennaTipDist=[];
    end


    if ~isempty(antennaTipDist)
        try
            antenaeEdgeMask=imdilate(oneAntenna,strel('disk',1))-imerode(oneAntenna,strel('disk',1));
            [ei,ej] = find(immultiply(It0,antenaeEdgeMask));
            antennaeBase0=sortrows([ei,ej],-1); %sort in descending order
            antennaeBase=antennaeBase0(1,:);
            fork2baseDist=Dd0(antennaeBase(1),antennaeBase(2));
        catch
            fork2baseDist=0;
            antennaeBase=forkPt00;
        end
        antL= antennaTipDist(1,3)-fork2baseDist;
    else
        antL= [];
        antennaeBase=[];
    end
    
    try
        [bolbMorph1,~,tipAntennae,tipAnt]=oneAntennaeMeasure03(mask, forkPt00, oneAntenna,antL, antennaeBase);
    catch
        bolbMorph1=-9999 + zeros(1, 4);
        tipAnt=[];
    end
    %judge if the antenna is on the right side or left side
    ss01 = regionprops(oneAntenna,'centroid');
    antCens=ss01.Centroid;
    if antCens(1)<realCen(1)
        rfAnt=1;
    else
        rfAnt=2;
    end
    bolbMorph=-9999 + zeros(2,4);
    if ~isempty(antL) && ~isempty(tipAnt)
        basetipdist=pdist([tipAnt;antennaeBase],'euclidean');
        curveDegree=antL/basetipdist;
    elseif isempty(antL) && ~isempty(tipAnt) %Added March 14, 2020
        if bolbMorph1(1,1)<0
            antL=-9999;
            curveDegree=-9999;
        else
            antL=bolbMorph1(1,1);
            curveDegree=bolbMorph1(1,4);
        end
    else  %Added March 14, 2020
        antL=-9999; %Added March 14, 2020
        curveDegree=-9999; %Added March 14, 2020
    end
    if antL>0 antLout=antL/scalelen*10;, else antLout=-9999;,  end;
    if bolbMorph1(2)>0 antWout=bolbMorph1(2:3)/scalelen*10;, else antWout=bolbMorph1(2:3);,  end;
    bolbMorph(rfAnt,:)=[antLout, antWout, curveDegree];
    antFig=oneAntenna;
            
            %%
%         %Plot antennae
%         box=regionprops(oneAntenna,'BoundingBox'); %for visualization
%         antBasePts=antennaeBase; %for visualization
%         textPos=box.BoundingBox(1:2); %for visualization
%         
%         antFigv=figure;
%         imshow(labeloverlay(double(mask*0.1+oneAntenna*0.5),It,'Colormap','autumn','Transparency',0));hold on;
% %     plot(antBasePts(:,2), antBasePts(:,1),'yO');
%         plot(antBasePts(:,2), antBasePts(:,1),'y*');
%             if rfAnt==1
%                 text(double(textPos(1))-30,double(textPos(2))-10,[num2str(round(bolbMorph(rfAnt,1),2)),'mm'],'color','g');
%             else
%                 text(double(textPos(1))+30,double(textPos(2))-10,[num2str(round(bolbMorph(rfAnt,1),2)),'mm'],'color','g');
%             end
%         hold off;   
%         
%         %%
%         %Plot tips
%         box=regionprops(oneAntenna,'BoundingBox'); %for visualization
%         antBasePts=antennaeBase; %for visualization
%         textPos=box.BoundingBox(1:2); %for visualization
%         
%         antFigv=figure;
%         imshow(labeloverlay(double(mask*0.1+oneAntenna*0.2+tipAntennae),It,'Colormap','autumn','Transparency',0));hold on;
% %     plot(antBasePts(:,2), antBasePts(:,1),'yO');
%         plot(antBasePts(:,2), antBasePts(:,1),'yx');
%             if rfAnt==1
%                 text(double(textPos(1))-30,double(textPos(2))-10,[num2str(round(bolbMorph(rfAnt,3),2)),'mm'],'color','w');
%             else
%                 text(double(textPos(1))+30,double(textPos(2))-10,[num2str(round(bolbMorph(rfAnt,3),2)),'mm'],'color','w');
%             end
%         hold off;   
end