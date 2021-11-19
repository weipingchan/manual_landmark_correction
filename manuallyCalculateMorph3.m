function manuallyCalculateMorph3(morph_mat_directory,Code_directory,spp_json_directory,Result_directory,antTrimPx, jumbo, SphingidaeOrNot)
% Code_directory=''; %Directory in which these MATLAB morphometric segmentation scripts live
% morph_mat_directory=''; %Directory of morph-seg.mat files pulled for problematic images; It can also be the dir of allBandsMask.mat
% Result_directory=''; %Directory of JSON files created during manual wing segmentation process
addpath(genpath(Code_directory)) %Add the library to the path
SphingidaeOrNot=0; %if the pre-posed wing image is desired desperately for Sphingidae-like specimens, manually change this to 1. Please note that the resulting Sphingidae shape analyses may be very error-prone.
% Turn off this warning "Warning: Image is too big to fit on screen; displaying at 33%"
% To set the warning state, you must first know the message identifier for the one warning you want to enable. 
warning('off', 'Images:initSize:adjustingMag');

%Read the file list in the Img_directory
ana_mode='morphDat';
img_ds = struct2dataset(dir(fullfile(morph_mat_directory,'*morph*.mat')));
img_listing=img_ds(:,1);
if isempty(img_listing) %If the dataset is not morph data, we use all bands data
    img_ds = struct2dataset(dir(fullfile(morph_mat_directory,'*_AllBandsMask.mat')));
    img_listing=img_ds(:,1);
    ana_mode='allBands';
end

%cd(Result_directory); %Move to the directory where the results will be stored.
disp('Start to create / find corresponding folders.');
%Create result directory
if ~exist(fullfile(Result_directory,'Shape_analysis'), 'dir')
    mkdir(fullfile(Result_directory,'Shape_analysis'));
end

subFolderList={'primary_key_pts','mask','main_characters','wing_shape_visualization','shape_context','shape_mask','wing_shape_matrices-seg','wing_segmentation'};

for fold=1:length(subFolderList)
    if ~exist(fullfile(Result_directory,'Shape_analysis',subFolderList{fold}), 'dir')
        mkdir(fullfile(Result_directory,'Shape_analysis',subFolderList{fold}));
    end
    disp(['corresponding folder ', subFolderList{fold}, ' is created / found.']);
end


vdlist={'dorsal','ventral'};
disp('Start to read images into memory');

for matinID=1:length(img_listing)
    if length(img_listing)>1
        matinname=img_listing.name{matinID};
    elseif length(img_listing)==1
        matinname=img_listing.name;
    end
%     if contains(matinname, 'dorsal')
%         side=1;
%         template0=strsplit(matinname,['_',vdlist{side}]);
%         template=template0{1};
% 
%     elseif contains(matinname, 'ventral')
%         side=2;
%         template0=strsplit(matinname,['_',vdlist{side}]);
%         template=template0{1};
%     end
    [barcode, side, flag]=file_name_decoder(matinname);
%     flag='_m';

    disp(['Start to manually analyze specimen: ', barcode,'_',vdlist{side},flag]);
    
    bflag=0;
    try
            matin=fullfile(morph_mat_directory,matinname);
            %read original matrix
            sppmat0=load(matin);
            fieldName=cell2mat(fieldnames(sppmat0));
            sppmat=sppmat0.(fieldName);
            clear sppmat0;

            if strcmp(ana_mode,'morphDat') %Use pre-run points
                ref=sppmat{1};
                tipPts=sppmat{5};
                refPts=sppmat{6};
                mask=sppmat{1};
                scale=sppmat{12};
            elseif strcmp(ana_mode,'allBands') %Use temporary points
                ref=imclearborder(sppmat{end-1});
%                 if nnz(ref)==0
%                     ref=sppmat{end-1};
%                 end
                stat=regionprops(ref,'Centroid','BoundingBox');
                tipPts=[stat.BoundingBox(1),stat.BoundingBox(2) ; stat.BoundingBox(1)+stat.BoundingBox(3),stat.BoundingBox(2)];
                refPts=zeros(6,2)+stat.Centroid;
                mask=imclearborder(sppmat{end-1});
                scale=sppmat{end};
            end
            ptNameList= {'L-F&H','L-F&B','R-F&B','R-F&H','R-H&B','L-H&B'};
            tipList={'tip-LF','tip-RF'};
        %     {'Left corner dividing fore & hindwings','Left corner dividing forewing & body','Right corner dividing  forewing & body',...
        %         'Right corner dividing fore & hindwings','Right corner dividing hindwing & body','Left corner dividing hindwing & body'};

        %Manually define those key points
        [newTipList,newRefList, bflag]=manuallyDefineKeyRefPts2(ref,tipPts,tipList, refPts(1:6,:), ptNameList);
        disp(['Key reference points in Img No. ', num2str(matinID),' out of  ',num2str(length(img_listing)),' have been manually defined.']);

        %Adjusted analyzing process
        [realCen, symAxis, symOrtho, boundingBox,  tipPts, refPts, wingParts, WingAxesSlopes, tformWingImg, shapeImg, allComImg, allinfo]=wingShapePreprocessingManualRef2(barcode,flag,mask,spp_json_directory,Result_directory,subFolderList,side,SphingidaeOrNot,newTipList,newRefList);
        
        bodyTrimPx=5;
%         jumbo=0;
%         antTrimPx=5; %Adjust this to extract bold antennae if necessary (20-30 for big antennae; 5-10 for thin antennae)
        if jumbo==0
            [bodyMask, bodyCharacters, antennaMask,antennaCharacters]=body_antenna_module2(mask, wingParts, refPts, tipPts, bodyTrimPx, antTrimPx, scale);
        else
            [bodyMask, bodyCharacters, antennaMask,antennaCharacters]=body_antenna_module_jumbo(mask, wingParts, refPts, tipPts, bodyTrimPx, antTrimPx, scale);
        end
        %segmented_img0=bodyMask*5+antennaMask*6+wingParts{1}*1+wingParts{2}*3+wingParts{3}*2+wingParts{4}*4;
        segmented_img=zeros(size(bodyMask));
        segmented_img(mask==1)=7;
        segmented_img(antennaMask==1)=6;
        segmented_img(bodyMask==1)=5;
        if side==1 %On dorsal side, forewings are dominant
            segmented_img(wingParts{2}==1)=3;
            segmented_img(wingParts{4}==1)=4;
            segmented_img(wingParts{1}==1)=1;
            segmented_img(wingParts{3}==1)=2;
        elseif side==2 %On ventral side, hindwings are dominant
            segmented_img(wingParts{1}==1)=2;
            segmented_img(wingParts{3}==1)=1;
            segmented_img(wingParts{2}==1)=4;
            segmented_img(wingParts{4}==1)=3;
        end 
        %0 is background
        %1 is left forewing from dorsal side
        %2 is right forewing from dorsal side
        %3 is left hindwing from dorsal side
        %4 is right hindwing from dorsal side
        %5 is body
        %6 is antennae
        %7 is uncertain

        %Create visualization Reference: https://color.adobe.com/zh/explore
        color_idx=[[0,[0,0,0]],
            [1,[242, 200, 5]],
            [2,[242, 135, 5]],
            [3,[242, 75, 153]],
            [4,[166, 3, 33]],
            [5,[242, 240, 240]],
            [6,[171, 5, 242]],
            [7,[40, 40, 40]]
            ];

        segmented_img_vis=label_segmentation_visualization(segmented_img,color_idx);
        %figure,imshow(segmented_img_vis)

      %Generate output matrix
        morphinfo=cell(0,15);
        morphinfo{1}=mask; %Original mask
        morphinfo{2}=realCen; %Original centroid.
        morphinfo{3}=symAxis; %The vertical symmetric axes
        morphinfo{4}=symOrtho; %The horizontal symmetric axes
        morphinfo{5}=tipPts; %tips of forewings
        morphinfo{6}=refPts; %Important segment points
        morphinfo{7}=wingParts; %right and left X fore and hindwings. Total: 4
        morphinfo{8}=WingAxesSlopes; %Slopes for the front edge of forewings and rear edge of hindwings
        morphinfo{9}=tformWingImg; %Re-posed wings
        morphinfo{10}=shapeImg; %Re-posed fore-hindwings in the same panel and related points in new coordinations.
        morphinfo{11}=allinfo; %Reference points before and after realignment
        morphinfo{12}=scale; %The scale bar (number of pixels = 1 cm)
        morphinfo{13}=segmented_img;
        morphinfo{14}=bodyCharacters; %The body length and width (in cm)
        morphinfo{15}=antennaCharacters; %Antennae length, width, bulb width, degree of curved (all in cm); first row is left one, second is right one.

        matoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{7},[barcode,'_',vdlist{side},flag,'_morph-seg.mat']);
        save(matoutname,'morphinfo'); %save the specimen matrix

        %Save the segmented image
        segvisoutname=fullfile(Result_directory,'Shape_analysis',subFolderList{8},[barcode,'_',vdlist{side},flag,'_wing_segmentation.jpg']);
        figseg=figure('visible', 'off');
        imshow(segmented_img_vis)
        %saveas(figmask, maskoutname);
        export_fig(figseg,segvisoutname, '-jpg','-r150');
        close(figseg);
        disp('An image showing img segmentation has been saved.');

        %Move those images to a subdirectory after they are analyzed
        finishedDir='done';
        if ~exist(fullfile(morph_mat_directory,finishedDir), 'dir')
            mkdir(fullfile(morph_mat_directory,finishedDir));
        end
        movefile(matin,fullfile(morph_mat_directory,finishedDir));
        disp('Analyzed images have been moved to "done" directory.');
    catch
        disp(['STOP analyzing specimen: ', barcode,'_',vdlist{side},flag]);
    end
    
    if bflag==1 %force stop the script
        break
    end
end
end