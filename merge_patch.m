%---  merge each patch optimal phase  ---%

disp('% -----  merge phase  -----% ');
load('patch.mat');     
npatch=size(patchlist,1);
% check Number of images processed -----------------------
if exist ('PATCH_001/oriprocessed_num.mat','file')
    load('PATCH_001/oriprocessed_num.mat', 'oriprocessed_num');
else
    oriprocessed_num=0;
end
patch_ph=zeros(rows,cols,total_num-oriprocessed_num);
% patch_SHPrecord=zeros(rows,cols,npart);
% patch_gamma=zeros(rows,cols,npart);
% patch_meancoh=zeros(rows,cols,npart);

clear comslc
switch length(channels)
    case 1
        comslc.com_channel1=zeros(rows,cols,no_subset);
    case 2
        comslc.com_channel1=zeros(rows,cols,no_subset);
        comslc.com_channel2=zeros(rows,cols,no_subset);
    case 3
        comslc.com_channel1=zeros(rows,cols,no_subset);
        comslc.com_channel2=zeros(rows,cols,no_subset);
        comslc.com_channel3=zeros(rows,cols,no_subset);
end

for i=1:npatch
    patchname=patchlist(i,:);
    if exist(patchname,'dir')    
        cd(patchname);
    else 
        error('Error, NO such patch: %s \n', patchname);
    end
    fprintf('Processing PATCH: %s \n',patchname);

    r1=patch_noover_local(i,1)-patch_over_local(i,1)+1 ;
    r2=patch_noover_local(i,1)-patch_over_local(i,1)+1+ patch_noover_local(i,3)-patch_noover_local(i,1);
    c1=patch_noover_local(i,2)-patch_over_local(i,2)+1 ;
    c2=patch_noover_local(i,2)-patch_over_local(i,2)+1+ patch_noover_local(i,4)-patch_noover_local(i,2);
    % merge SHPrecord, gamma and meancoh -----------------------
%     SHPrecord=stamps_read('SHPrecord', (patchsize(i,1))*(patchsize(i,2)), 'float32');
%     SHPrecord=reshape(SHPrecord,(patchsize(i,1)),(patchsize(i,2)),size(SHPrecord,2));
%     patch_SHPrecord(patch_noover_local(i,1):patch_noover_local(i,3),patch_noover_local(i,2):patch_noover_local(i,4),:) = SHPrecord(r1:r2,c1:c2,:);
     
    % merge patch phase -------------------------------------
    opt_phase=stamps_read('opt_phase', (patchsize(i,1))*(patchsize(i,2)), 'float32');
    patch_phase=reshape(opt_phase,(patchsize(i,1)),(patchsize(i,2)),size(opt_phase,2));
    patch_ph(patch_noover_local(i,1):patch_noover_local(i,3),patch_noover_local(i,2):patch_noover_local(i,4),:) = patch_phase(r1:r2,c1:c2,:);
 

    cd ..;   
end

