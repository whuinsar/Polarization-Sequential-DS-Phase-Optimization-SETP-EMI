%% Environmental parameter settings and check ----------------------
mkdir('opt_diff'); % File path to store the optimized Interferograms
mkdir('com_slc'); % File path to store the compressed SLCs

CalWin=[2*hW_l+1 2*hW_w+1];

% start parallel if it not running------------
if isempty(gcp('nocreate'))
    mypool=parpool(cores); 
end
% ----------------------------------------------
%% check the total number of images
imgpath=[workpath '/' char(channels(1)) '/rmli/'];
tag_files = dir([imgpath,'*','.rmli']);
total_num=size(tag_files,1); % Total number of SLC images
npart=ceil(total_num/interval); % The number of mini-time series images generated 
% by segmentation in the time dimension
% ----------------------------------------------
if crop_flag~=1
   r0=1;c0=1;
   [rN,cN]=size(freadbkj([imgpath,tag_files(1).name],nlines,'float32'));
   rows=rN-r0+1;  
   cols=cN-c0+1;
end
rows=rN-r0+1;  
cols=cN-c0+1;

%% generate spatial patches
% overlap between adjacent patches
overlap_az=2*hW_l+1;
overlap_rg=2*hW_w+1;
% create patch list
[patchlist, patch_noover_local, patch_over_local, patch_noover, patch_over, patchsize]=...
    Createpatch(r0, c0, rN, cN,paz,  prg,  overlap_az, overlap_rg);
save('patch.mat','nlines','hW_l','hW_w','rN','cN', 'overlap_az','overlap_rg','paz','prg',...
    'patchlist','patch_noover_local','patch_over_local','patch_noover','patch_over','patchsize'); 
%%
%---check mask ------------------
% if mask_flag==1
%    load('mask.mat');
% end
% -------------------------------
%%






