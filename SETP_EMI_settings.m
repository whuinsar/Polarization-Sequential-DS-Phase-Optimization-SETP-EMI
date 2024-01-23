%% parms setting ------------------------------------------------
% work path ----------------------------------------------------
workpath=['E:\project\simu_exp'];
% --------------------------------------------------------------------
% number of lines in slc (The size of SLC, can be found in slc.par file)
nlines=301;

% channels:----------------------------------------------------
%channels=["vv";"hh";"vh"]; or channels=["vv"]; or channels=["vv"; "vh"];
channels=["vv"; "vh"];

% Crop range ----------------------------------------------------
crop_flag=0; % Crop: 1  / no crop: 0
r0=400; %start number on rows
rN=700; %end number on rows
c0=1300; %start number on cols
cN=1700; %end number on cols

% The size of spatial patches
paz=1;
prg=1;
% The size of temporal mini-subsets(numbers of SLCs in subset) -----------------------------
interval=10;

% SHP window settings (for HTCI test) ----------------------------------
hW_w=5; %half window size from range direction
hW_l=5;  %half window size from azimuth direction
Alpha=0.05; % significance level

%Reference image ----------------------------------------------------
masterID=1;

% Minimum threshold for DS point identification ------------------
%(the number of SHPs greater than this threshold is considered a DS point)
minshp=20;

%Set the number of CPU parallel cores
cores=12;


% ----------------------------------------------