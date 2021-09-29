function [newTipList,newRefList]=manuallyDefineKeyRefPts2(ref,tipPts,tipList,refPts,ptNameList)
    %ptNameList= {'L-F&H','L-F&B','R-F&B','R-F&H','R-H&B','L-H&B'};
    %tipList={'tip-LF','tip-RF'};
    %     {'Left corner dividing Fore & Hing Wings','Left corner dividing Fore Wing & Body','Right corner dividing  Fore Wing & Body',...
    %         'Right corner dividing Fore & Hing Wings','Right corner dividing Hind Wing & Body','Left corner dividing Hind Wing & Body'};
    [B0,~]=bwboundaries(ref);
    if length(B0)==1
        edgePt=B0{1};
    else
        edgePt=[];
        for obj=1:length(B0)
            edgePt=[edgePt; B0{obj}]; %Make sure the edges of multiple parts are all in the list
        end
    end
    
    refimg=figure;
    imshow(ref);
    hold on;

    bflag=0;
    newRefList=[];
    newTipList=[];
    hid=1;
    while 1
        if hid<=size(refPts,1)
            hit = interactPt(refPts,hid,ptNameList,'r');
        else
            hit = interactPt(tipPts,hid-size(refPts,1),tipList,'g');
        end
        
        while 1
            while 1
            k = waitforbuttonpress;
            %13 enter
            % 32 space
            %104 f
            %114 r
                if k==1 %exclude mouse click
                    returnValue = double(get(gcf,'CurrentCharacter'));
                    break                
                end
            end

            if returnValue==13 %Press Enter to confirm using this box for adjusting
                break
%             elseif returnValue==102 %Press 'f' to force exit the entire manual process
%                 bflag=1;
%                 break
            elseif returnValue==114 %Press 'r' to redo the manual process
                bflag=2;
                break
            end

        end
        
        if  bflag==0
%             pause; %Press enter if finished with the adjusting
            position0 = getPosition(hit);
            delete(hit);
            position=findClosestPt([edgePt(:,2),edgePt(:,1)], position0); %move the point to the edge
            if hid<=size(refPts,1)
                newRefList=[newRefList; position];
                plot(position(:,1),position(:,2), 'r+','MarkerSize',8);
                plot(position(:,1),position(:,2), 'ro','MarkerSize',3);
            else
                newTipList=[newTipList; position];
                plot(position(:,1),position(:,2), 'gx','MarkerSize',12);
                plot(position(:,1),position(:,2), 'go','MarkerSize',8);
            end
            
            if hid==size(refPts,1)+size(tipPts,1)
                break
            end
            hid=hid+1;
%         elseif  bflag==1
%             break
        elseif bflag==2
            close(refimg);
            refimg=figure;
            imshow(ref);
            hold on;
            newRefList=[];
            newTipList=[];
            bflag=0;
            hid=1;
        end
    end
    hold off;
    close(refimg);
    clear('refimg');
end