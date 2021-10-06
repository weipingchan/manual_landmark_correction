function [bolbMorph,antFig]=twoAntennaMeasure5(mask, body, antenna1,distanceAnt0,forkPt00, bodyraw, upperSegMask, scalelen)
            %Head mask until two antennae spread
            radii=10;
            while 1
                Cmask = createCirclesMask(mask, forkPt00, radii);
                revCmask= imcomplement(immultiply(Cmask, bodyraw));
                cAntenna00=bwareaopen(immultiply(antenna1,revCmask),100);
                cAntenna0=immultiply(logical(bwareaopen(immultiply(cAntenna00, upperSegMask),200)),  imcomplement(bodyraw));
                [labeledAntenna0, numberOfBlobs0] = bwlabel(cAntenna0);                
                if numberOfBlobs0>2 %if there are more than two objects, take the two closest to the top
                    poslist=[];
                    for k0=1:numberOfBlobs0
                        obj=ismember(labeledAntenna0, k0);
                        [x0,~]=find(obj);
                        poslist=[poslist; min(x0)];
                    end
                    poslist0=poslist;
                    idxlist=[];
                     for ii=1:2
                      [~,idx] = min(poslist0);
                      % remove the smallest value for the next iteration:
                      poslist0(idx) = 9999;
                      idxlist=[idxlist,idx];
                     end                   
                    cAntenna=ismember(labeledAntenna0, idxlist(1))+ismember(labeledAntenna0, idxlist(2));
                else
                    cAntenna=cAntenna0;
                end
                
                [labeledAntenna, numberOfBlobs] = bwlabel(cAntenna);
            %     cc = bwconncomp(cAntenna);
            %     numberO  = cc.NumObjects;
                if radii>= min(size(mask))/2
                    flag=1;
                    break
                elseif numberOfBlobs>=2
                    flag=0;
                    break
                else
                    radii=radii+5;
                end
            end

            if flag==0 %successfully finish the previous step
                It = bwmorph(antenna1,'thin','inf');
                [Bbx,Bby] =  find(bwmorph(It,'branchpoints'));
                %select branch point and convert to index
                dist = sqrt(sum(bsxfun(@minus, [Bbx,Bby], flip(forkPt00)).^2,2));
                [~,minLoc]=min(abs(dist-mean(dist(dist<mean(dist)+std(dist)*1/2,:)))); %find the one closest to the mean of those points near forkPt
                Bbi=sub2ind(size(antenna1),Bbx(minLoc),Bby(minLoc)); %The forkPoint we want

                [i,j] = find(immultiply(bwmorph(It,'endpoints'),cAntenna));  %Calculation with refined end points
                Dd = bwdistgeodesic(It,Bbi,'quasi');
                distanceAnt=[];
                for n = 1:numel(i)
                    distanceAnt=[distanceAnt; [i(n),j(n),Dd(i(n),j(n))]];
                end         
                %%
                 try %prevent potential error and interruption

                    %%
                    if size(distanceAnt,1)>2
                        %Cluster points long distances from the fork
                        Dist0=sqrt(sum(bsxfun(@minus, distanceAnt0(:,1:2), flip(forkPt00)).^2,2));
                        Dist = sqrt(sum(bsxfun(@minus, distanceAnt(:,1:2), flip(forkPt00)).^2,2));
                        antennaTipDist0= sortrows(distanceAnt(Dist>mean(Dist0)+std(Dist0)*1/2,:),2); %first one is the left one
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
                
                %%
                antFig2=cell(0,2);
                antBasePts=zeros(2, 2); %for visualization
                bolbMorph=-9999 + zeros(2, 4);
                textPos=zeros(2, 2); %for visualization
                tipHighlight=cell(0,2); %for visualization
                for k = 1 : numberOfBlobs
                    oneAntenna = bwareafilt(bwareaopen(logical(immultiply(ismember(labeledAntenna, k),imcomplement(body))),100),1);

                    if nnz(oneAntenna)~=0
                    %%
                        bondAnt = bwboundaries(oneAntenna);
                        if size(antennaTipDist,1)==2
                            d2pt1=mean(sqrt(sum(bsxfun(@minus, bondAnt{1}, antennaTipDist(1,1:2)).^2,2))); %Left antenna
                            d2pt2=mean(sqrt(sum(bsxfun(@minus, bondAnt{1}, antennaTipDist(2,1:2)).^2,2))); %Right antenna
                            if d2pt1<d2pt2
                                rfAnt=1;
                            else
                                rfAnt=2;
                            end
                        else
                            %Use one antenna module to judge if the antenna is the right side or left side
                            ss01 = bondAnt{1};
                            antCens=mean(ss01);
                            if antCens(1)<forkPt00(1) %Use forkPt00 to replace realCen
                                rfAnt=1;
                            else
                                rfAnt=2;
                            end
                        end
                          
                        box=regionprops(oneAntenna,'BoundingBox'); %for visualization
        %                 tipAnt= antennaTipDist(rfAnt,1:2);
                        try
                            antenaeEdgeMask=imdilate(oneAntenna,strel('disk',1))-imerode(oneAntenna,strel('disk',1));
                            [ei,ej] = find(immultiply(It,antenaeEdgeMask));
                            antennaeBase0=sortrows([ei,ej],-1); %sort in descending order
                            antennaeBase=antennaeBase0(1,:);
                            fork2baseDist=Dd(antennaeBase(1),antennaeBase(2));
                        catch
                            fork2baseDist=0;
                            antennaeBase=forkPt00;
                        end
                        antL= antennaTipDist(rfAnt,3)-fork2baseDist;
                        try
                            [bolbMorph1,~,tipAntennae,tipAnt]=oneAntennaeMeasure03(mask, forkPt00, oneAntenna,antL,antennaeBase);
                        catch
                            bolbMorph1=-9999 + zeros(1, 4);
                            tipAntennae=zeros(size(mask));
                            tipAnt=[];
                        end
                        
                        if ~isempty(antL) && ~isempty(tipAnt)
                            basetipdist=pdist([tipAnt;antennaeBase],'euclidean');
                            curveDegree=antL/basetipdist;
                        elseif isempty(antL) && ~isempty(tipAnt) 
                            if bolbMorph1(1,1)<0 
                                antL=-9999; 
                                curveDegree=-9999; 
                            else 
                                antL=bolbMorph1(1,1); 
                                curveDegree=bolbMorph1(1,4); 
                            end 
                        else 
                            antL=-9999; 
                            curveDegree=-9999;
                        end
                        if antL>0 antLout=antL/scalelen*10;, else antLout=-9999;,  end;
                        if bolbMorph1(2)>0 antWout=bolbMorph1(2:3)/scalelen*10;, else antWout=bolbMorph1(2:3);,  end;
                        bolbMorph(rfAnt,:)=[antLout, antWout, curveDegree];
                        antFig2{k}=oneAntenna;
                        antBasePts(rfAnt,:)=antennaeBase;
                        textPos(rfAnt,:)=box.BoundingBox(1:2);
                        tipHighlight{rfAnt}=tipAntennae;
                    else
                        antFig2{k}=oneAntenna;
                        antBasePts=forkPt00;
                    end
                    clear('antL');
                end
                antFig=double(mask*0.1+ (antFig2{1}+antFig2{2})*0.8);
            else %There are actually no antennae
                disp('Oops! NO antennae were found.');
                bolbMorph=-9999 + zeros(2, 4);
                antFig=double(mask*0.1);
            end
            %%
%                 %Plot antenna
%                 antFigv=figure;
%                 imshow(labeloverlay(double(mask*0.1+(antFig2{1}+antFig2{2})*0.5),It,'Colormap','autumn','Transparency',0));hold on;
% %                 plot(antBasePts(:,2), antBasePts(:,1),'yO');
%                 plot(antBasePts(:,2), antBasePts(:,1),'y*');
%                 for n = 1:size(antennaTipDist,1)
%                     if n==1
%                         text(double(textPos(n,1))-30,double(textPos(n,2))-10,[num2str(round(bolbMorph(n,1),2)),'mm'],'color','g');
%                     else
%                         text(double(textPos(n,1))+30,double(textPos(n,2))-10,[num2str(round(bolbMorph(n,1),2)),'mm'],'color','g');
%                     end
%                 end
%                 hold off;        
%                 
%                 %%
%                 %Plot tips
%                 antFigv=figure;
%                 imshow(labeloverlay(double(mask*0.1+(antFig2{1}+antFig2{2})*0.2+(tipHighlight{1}+tipHighlight{2})),It,'Colormap','autumn','Transparency',0));hold on;
% %                 plot(antBasePts(:,2), antBasePts(:,1),'yO');
%                 plot(antBasePts(:,2), antBasePts(:,1),'yx');
%                 for n = 1:size(antennaTipDist,1)
%                     if n==1
%                         text(double(textPos(n,1))-30,double(textPos(n,2))-10,[num2str(round(bolbMorph(n,3),2)),'mm'],'color','w');
%                     else
%                         text(double(textPos(n,1))+30,double(textPos(n,2))-10,[num2str(round(bolbMorph(n,3),2)),'mm'],'color','w');
%                     end
%                 end
%                 hold off;       
end