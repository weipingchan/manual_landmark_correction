function manuallyCalculateMorph3(morph_mat_directory,Code_directory,spp_json_directory,Result_directory,antTrimPx, jumbo, SphingidaeOrNot)
% Code_directory='D:\Milk desk\Dropbox\Harvard\Coloration_research\Multi_spectra_processing/shape_analysis_v1';
% morph_mat_directory='D:\Milk desk\Dropbox\Harvard\Coloration_research\Drawer_result\Shape_analysis\wing_shape_matrices';
% Result_directory='D:\Milk desk\Dropbox\Harvard\Coloration_research\Drawer_result\Shape_analysis\wing_shape_matrices';
addpath(genpath(Code_directory)) %Add the library to the path

% Sphingidae=SphingidaeOrNot;
Sphingidae=0; %if the pre-posed-wing image is desired desperatedly for Sphingidae-like specimens, manually turn this to 1, but there may have many errors in the following process.
% Turn off this warning "Warning: Image is too big to fit on screen; displaying at 33% "
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
    
    try
            matin=fullfile(morph_mat_directory,matinname);
            %read original matrice
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
                stat=regionprops(ref,'Centroid','BoundingBox');
                tipPts=[stat.BoundingBox(1),stat.BoundingBox(2) ; stat.BoundingBox(1)+stat.BoundingBox(3),stat.BoundingBox(2)];
                refPts=zeros(6,2)+stat.Centroid;
                mask=imclearborder(sppmat{end-1});
                scale=sppmat{end};
            end
            ptNameList= {'L-F&H','L-F&B','R-F&B','R-F&H','R-H&B','L-H&B'};
            tipList={'tip-LF','tip-RF'};
        %     {'Left corner dividing Fore & Hing Wings','Left corner dividing Fore Wing & Body','Right corner dividing  Fore Wing & Body',...
        %         'Right corner dividing Fore & Hing Wings','Right corner dividing Hind Wing & Body','Left corner dividing Hind Wing & Body'};

        %Manually define those key points
        [newTipList,newRefList]=manuallyDefineKeyRefPts2(ref,tipPts,tipList, refPts(1:6,:), ptNameList);
        disp(['Key reference points in Img No. ', num2str(matinID),' out of  ',num2str(length(img_listing)),' have been manually defined.']);

        %Adjusted analyzing process
        [realCen, symAxis, symOrtho, boundingBox,  tipPts, refPts, wingParts, WingAxesSlopes, tformWingImg, shapeImg, allComImg, allinfo]=wingShapePreprocessingManualRef2(barcode,flag,mask,spp_json_directory,Result_directory,subFolderList,side,Sphingidae,newTipList,newRefList);
        
        bodyTrimPx=5;
%         jumbo=0;
%         antTrimPx=5; %Adjust this to extart bold antenna if necessary (20-30 for big antenna; 5-10 for thin antenna)
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
        if side==1 %In dorsal side, fore wings are dominant
            segmented_img(wingParts{2}==1)=3;
            segmented_img(wingParts{4}==1)=4;
            segmented_img(wingParts{1}==1)=1;
            segmented_img(wingParts{3}==1)=2;
        elseif side==2 %In ventral side, hind wings are dominant
            segmented_img(wingParts{1}==1)=2;
            segmented_img(wingParts{3}==1)=1;
            segmented_img(wingParts{2}==1)=4;
            segmented_img(wingParts{4}==1)=3;
        end 
        %0 is background
        %1 is left fore wing from dorsal side
        %2 is right fore wing from dorsal side
        %3 is left hind wing from dorsal side
        %4 is right hind wing from dorsal side
        %5 is body
        %6 is antenna
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
        morphinfo{3}=symAxis; %The vertical symatric axes
        morphinfo{4}=symOrtho; %The horizontal symatric axes
        morphinfo{5}=tipPts; %tips of forewings
        morphinfo{6}=refPts; %Important segment points
        morphinfo{7}=wingParts; %right and left X fore and hind wings. Total: 4
        morphinfo{8}=WingAxesSlopes; %Slopes for the front edge of fore wings and rare edge of hind wings
        morphinfo{9}=tformWingImg; %Re-posed wings
        morphinfo{10}=shapeImg; %Re-posed fore-hind wings in the same panel and related points in new corrdinatios.
        morphinfo{11}=allinfo; %Reference points before and after realignment
        morphinfo{12}=scale; %The scale bar (number of pixels = 1 cm)
        morphinfo{13}=segmented_img;
        morphinfo{14}=bodyCharacters; %The body length and width (in cm)
        morphinfo{15}=antennaCharacters; %Antennae length, width, bolb width, degree of curved (all in cm); first row is left one, second is right one.

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

        %Move those images having been analyzed to a subdirectory
        finishedDir='done';
        if ~exist(fullfile(morph_mat_directory,finishedDir), 'dir')
            mkdir(fullfile(morph_mat_directory,finishedDir));
        end
        movefile(matin,fullfile(morph_mat_directory,finishedDir));
        disp('Images analyzed have been moved to done directory.');
    catch
        disp(['STOP analyzing specimen: ', barcode,'_',vdlist{side},flag]);
    end
end
end