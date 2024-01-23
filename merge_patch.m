%---  merge each patch optimal phase  ---%

disp('% -----  merge phase  -----% ');
load('patch.mat');     
npatch=size(patchlist,1);
patch_ph=zeros(rows,cols,total_num);
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
    

    % merge patch phase -------------------------------------
    opt_phase=stamps_read('opt_phase', (patchsize(i,1))*(patchsize(i,2)), 'float32');
    patch_phase=reshape(opt_phase,(patchsize(i,1)),(patchsize(i,2)),size(opt_phase,2));
   
    patch_ph(patch_noover_local(i,1):patch_noover_local(i,3),patch_noover_local(i,2):patch_noover_local(i,4),:) = patch_phase(r1:r2,c1:c2,:);
 
    % merge compressed SLC -------------------------------------
    switch length(channels)
        case 1
            com_channel1=stamps_read('com_channel1',(patchsize(i,1))*(patchsize(i,2)), 'cpxfloat32');
            com_channel1=reshape(com_channel1,(patchsize(i,1)),(patchsize(i,2)),size(com_channel1,2));
            comslc.com_channel1(patch_noover_local(i,1):patch_noover_local(i,3),patch_noover_local(i,2):patch_noover_local(i,4),:)=com_channel1(r1:r2,c1:c2,:);

        case 2
            com_channel1=stamps_read('com_channel1',(patchsize(i,1))*(patchsize(i,2)), 'cpxfloat32');
            com_channel2=stamps_read('com_channel2',(patchsize(i,1))*(patchsize(i,2)), 'cpxfloat32');
            com_channel1=reshape(com_channel1,(patchsize(i,1)),(patchsize(i,2)),size(com_channel1,2));
            com_channel2=reshape(com_channel2,(patchsize(i,1)),(patchsize(i,2)),size(com_channel2,2));
            comslc.com_channel1(patch_noover_local(i,1):patch_noover_local(i,3),patch_noover_local(i,2):patch_noover_local(i,4),:)=com_channel1(r1:r2,c1:c2,:);
            comslc.com_channel2(patch_noover_local(i,1):patch_noover_local(i,3),patch_noover_local(i,2):patch_noover_local(i,4),:)=com_channel2(r1:r2,c1:c2,:);
        case 3
            com_channel1=stamps_read('com_channel1',(patchsize(i,1))*(patchsize(i,2)), 'cpxfloat32');
            com_channel2=stamps_read('com_channel2',(patchsize(i,1))*(patchsize(i,2)), 'cpxfloat32');
            com_channel3=stamps_read('com_channel3',(patchsize(i,1))*(patchsize(i,2)), 'cpxfloat32');
            com_channel1=reshape(com_channel1,(patchsize(i,1)),(patchsize(i,2)),size(com_channel1,2));
            com_channel2=reshape(com_channel2,(patchsize(i,1)),(patchsize(i,2)),size(com_channel2,2));
            com_channel3=reshape(com_channel3,(patchsize(i,1)),(patchsize(i,2)),size(com_channel3,2));
            comslc.com_channel1(patch_noover_local(i,1):patch_noover_local(i,3),patch_noover_local(i,2):patch_noover_local(i,4),:)=com_channel1(r1:r2,c1:c2,:);
            comslc.com_channel2(patch_noover_local(i,1):patch_noover_local(i,3),patch_noover_local(i,2):patch_noover_local(i,4),:)=com_channel2(r1:r2,c1:c2,:);
            comslc.com_channel3(patch_noover_local(i,1):patch_noover_local(i,3),patch_noover_local(i,2):patch_noover_local(i,4),:)=com_channel3(r1:r2,c1:c2,:);
    end

    cd ..;   
end

