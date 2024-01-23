%% Environmental parameter settings and check ----------------------
mkdir('opt_diff'); % File path to store the optimized Interferograms
mkdir('com_slc'); % File path to store the compressed SLCs

CalWin=[2*hW_l+1 2*hW_w+1];

% start parallel ------------
mypool=parpool(cores);
% ----------------------------------------------
%% check the total number of images
imgpath=[workpath '/' char(channels(1)) '/rmli/'];
tag_files = dir([imgpath,'*','.rmli']);
total_num=size(tag_files,1);
npart=ceil(total_num/interval); % Split into npart subsets
% ----------------------------------------------
if crop_flag~=1
   r0=1;c0=1;
   [rN,cN]=size(freadbkj([imgpath,tag_files(1).name],nlines,'float32'));
   rows=rN-r0+1;  
   cols=cN-c0+1;
end
rows=rN-r0+1;  
cols=cN-c0+1;

% tic
%% generate spatial patches
% overlap between adjacent patches
overlap_az=2*hW_l+1;
overlap_rg=2*hW_w+1;
% create
[patchlist, patch_noover_local, patch_over_local, patch_noover, patch_over, patchsize]=...
    Createpatch(r0, c0, rN, cN,paz,  prg,  overlap_az, overlap_rg);
save('patch.mat','nlines','hW_l','hW_w','rN','cN', 'overlap_az','overlap_rg','paz','prg',...
    'patchlist','patch_noover_local','patch_over_local','patch_noover','patch_over','patchsize'); 

%%






