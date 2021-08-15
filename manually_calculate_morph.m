Code_directory='E:\WP_work\Dropbox\Harvard\Coloration_research\Multi_spectra_processing\Method_summary\matlab_scripts_organized\Image_Segmentation\manual_landmark_correction';
morph_mat_directory='E:\WP_work\Dropbox\Harvard\Coloration_research\Multi_spectra_processing\Method_summary\Examplar_imgs\morphology_analysis\manual_correction\wing_shape_matrices-reSeg'; %It can also be the dir of allBandsMask.mat
spp_json_directory='E:\WP_work\Dropbox\Harvard\Coloration_research\Multi_spectra_processing\Method_summary\Examplar_imgs\spp_beautiful_RGB_Imgs\segmented';
Result_directory='E:\WP_work\Dropbox\Harvard\Coloration_research\Multi_spectra_processing\Method_summary\Examplar_imgs\morphology_analysis\manual_correction';
addpath(genpath(Code_directory)) %Add the library to the path
SphingidaeOrNot=0;
antTrimPx=10; %Adjust this to extart bold antenna if necessary (20-30 for big antenna; 5-10 for thin antenna; <5 for tiny specimens)
jumboAntenna=0; % set to 1 for the jumbo antenna which have multiple forks (jumbo mode will dilate a little to remove those forks)

% Turn off this warning "Warning: Image is too big to fit on screen; displaying at 33% "
% To set the warning state, you must first know the message identifier for the one warning you want to enable. 
warning('off', 'Images:initSize:adjustingMag');
manuallyCalculateMorph3(morph_mat_directory,Code_directory,spp_json_directory,Result_directory,antTrimPx, jumboAntenna, SphingidaeOrNot);
%Plot order is as follows:
%Left corner dividing Fore & Hing Wings
%Left corner dividing Fore Wing & Body
%Right corner dividing  Fore Wing & Body
%Right corner dividing Fore & Hing Wings
%Right corner dividing Hind Wing & Body
%Left corner dividing Hind Wing & Body
%Left fore wing tip
%right fore wing tip


%Press Enter to confirm using the point showing on the screen
%Press 'f' to force exit the entire manual process
%Press 'r' to redo the manual process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%