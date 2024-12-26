%% Overall Control for SETP-EMI DS Preprocessing
%###########################################
% Step1: initialize parameters
% Step2: Check historical variables for Sequential processing
% Step3: SETP-EMI processing 
% Step4: Merge patch 
% Step5: Save final results and compressed SLC for Future processing
%############################################
%   =================================================================
% created by Yian Wang,20230303
% modified by Yian Wang,20231008   Re-organize the data structure, modularized
% modified by Yian Wang,20231015   add parallel processing 
% modified by Yian Wang,20231023   add Spatial block processing
%   =================================================================
clear;clc;
%% parms setting -----------------------------------------------------
disp(['Step1: initialize SETP-EMI parameters']);
SETP_EMI_settings; % The user needs to manually modify this file
disp(['Step1: Done!']);
%% Environmental parameter settings and check-- ----------------------
disp(['Step2: Check historical variables for Sequential processing']);
SETP_EMI_check;
disp(['Step2: Done!']);
%% SETP-EMI processing -----------------------------------------------
disp(['Step3: SETP-EMI processing']);
SETP_EMI_process;
disp(['Step3: Done!']);
%% Merge patch -------------------------------------------------------
disp(['Step4: merge patch']);
merge_patch;
disp(['Step4: Done!']);
%% Save final results ------------------------------------------------
disp(['Step5: Save final results and compressed SLC for Future processing']);
SETP_EMI_save;
disp(['Step5: Done!']);
%% --------------------------------------------------------------------