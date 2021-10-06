Code_directory=''; %Directory in which MATLAB morphometric segmentation scripts live
morph_mat_directory=''; %Directory of morph-seg.mat files pulled for problematic images; It can also be the dir of allBandsMask.mat
spp_json_directory=''; %Directory of JSON files created during manual wing segmentation process
Result_directory=''; %Directory to which you would like to write shape analysis results
addpath(genpath(Code_directory)) %Add the library to the path
SphingidaeOrNot=0; %If the pre-posed wing image is desired desperately for Sphingidae-like specimens, manually change this to 1. Please note that the resulting Sphingidae shape analyses may be very error-prone.
antTrimPx=10; %Adjust this to extract bold antennae if necessary (20-30 for big antennae; 5-10 for thin antennae; <5 for tiny specimens)
jumboAntenna=0; % set to 1 for the jumbo antennae which have multiple forks (jumbo mode will dilate a little to remove those forks)

% Turn off this warning "Warning: Image is too big to fit on screen; displaying at 33% "
% To set the warning state, you must first know the message identifier for the one warning you want to enable. 
warning('off', 'Images:initSize:adjustingMag');
manuallyCalculateMorph3(morph_mat_directory,Code_directory,spp_json_directory,Result_directory,antTrimPx, jumboAntenna, SphingidaeOrNot);
%Plot order is as follows:
%Left corner dividing fore & hindwings
%Left corner dividing forewing & body
%Right corner dividing forewing & body
%Right corner dividing fore & hindwings
%Right corner dividing hindwing & body
%Left corner dividing hindwing & body
%Left forewing tip
%Right forewing tip


%Press Enter to confirm using the point showing on the screen
%Press 'r' to redo the manual process
%Press 'f' to force-exit the entire manual process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%