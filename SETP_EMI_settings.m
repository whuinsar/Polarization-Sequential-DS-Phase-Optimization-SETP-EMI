%% parms setting #######################################################
% The files and structures that users need to prepare shold 
% according to the SETP-EMI Program Guide Book

%% work path ------------------------------------------------------------
workpath=['E:\代码及资料交接\Polarization-Sequential-DS-Phase-Optimization-SETP-EMI-main']; 
% Find the explanation in SETP-EMI Program Guide Book
%% Set the number of CPU parallel cores ---------------------------------
cores=8; %Number of CPU parallel processing threads
%% channels--------------------------------------------------------------
% Polarization SAR data by the user: channels=["vv";"hh";"vh"]; 
% or channels=["vv"; "vh"]; or channels=["vv"];
channels=["vv"; "vh"];
%% number of lines in slc -----------------------------------------------
% (The size of SLC, can be found in slc.par file)
nlines=301;
%% Reference image ------------------------------------------------------
masterID=1;
%% Crop range -----------------------------------------------------------
crop_flag=0; % Crop: 1  / No crop: 0
r0=350;  % start number on rows
rN=650;  % end number on rows
c0=1200; % start number on cols
cN=1700; % end number on cols

%% The size of spatial patches ------------------------------------------
% The number of blocks in the azimuth and range directions
paz=2; %azimuth
prg=2; %range
% Subset size (numbers of SLCs in subset) -------------------------------
interval=10;

%% SHP window settings (for HTCI test) ----------------------------------
hW_w=10; %half window size from range direction
hW_l=5;  %half window size from azimuth direction
Alpha=0.05; % significance level
% Minimum threshold for DS point identification -------------------------
% The number of SHPs greater than this threshold is considered a DS point
minshp=25;
% ----------------------------------------------