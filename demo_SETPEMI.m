%% Overall Control for SETP-EMI DS Preprocessing
%###########################################
% Step1: initialize parameters
% Step2: Check historical variables for Sequential processing
% Step3: SETP-EMI processing (Sequential Read SAR data, Sequential Identify SHP,Sequential TP-EMI)
% Step4: Save final results and compressed SLC for Future processing
%############################################

%   =================================================================
% created by Yian Wang,20230303
% modified by Yian Wang,20231008   Re-organize the data structure, modularized
% modified by Yian Wang,20231015   add parallel processing 
% modified by Yian Wang,20231023   add Spatial block processing
%   =================================================================
clear;clc;
%% parms setting ----------------------------------------------------
disp(['Step1: initialize SETP-EMI parameters']);
SETP_EMI_settings; %
disp(['Step1: Done!']);

%% Environmental parameter settings and check------------------------
disp(['Step2: Check historical variables for Sequential processing']);
SETP_EMI_check;
disp(['Step2: Done!']);

%% SETP-EMI processing ----------------------------------------------
tic
disp(['Step3: SETP-EMI processing']);
SETP_EMI_process;
disp(['Step3: Done!']);
timecost=toc;
%  processed_num=0;stamps_save('processed_num.mat',processed_num);clear;clc;
%% Merge patch -----------------------------------------------
disp(['Step4: merge patch']);
merge_patch;
disp(['Step4: Done!']);

%% Save final results -----------------------------------------------
disp(['Step5: Save final results and compressed SLC for Future processing']);
SETP_EMI_save;
disp(['Step5: Done!']);

%% check result (Optional)
