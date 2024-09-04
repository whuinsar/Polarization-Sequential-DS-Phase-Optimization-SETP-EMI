%% parms setting -----------------------------------------------------

% work path -----------------------------------------------------------
workpath=['E:/project/xueyou'];
% channels:------------------------------------------------------------
%channels=["vv";"hh";"vh"]; or channels=["vv"; "vh"]; or channels=["vv"];
channels=["vv"; "vh"];
% --------------------------------------------------------------------
% number of lines in slc (The size of SLC, can be found in slc.par file)
nlines=901;
%Reference image ------------------------------------------------------
masterID=58;

% Crop range ----------------------------------------------------------
crop_flag=1; % Crop: 1  / no crop: 0
r0=350; %start number on rows
rN=650; %end number on rows
c0=1200; %start number on cols
cN=1700; %end number on cols

% The size of spatial patches ------------------------------------------
paz=1;
prg=1;
%subset size (numbers of SLCs in subset) -------------------------------
interval=10;

% SHP window settings (for HTCI test) ----------------------------------
hW_w=10; %half window size from range direction
hW_l=5;  %half window size from azimuth direction
Alpha=0.05; % significance level
% Minimum threshold for DS point identification ------------------------
%(the number of SHPs greater than this threshold is considered a DS point)
minshp=25;

%Set the number of CPU parallel cores ----------------------------------
cores=4;



% ----------------------------------------------