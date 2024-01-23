%% SETP-EMI procesing----------------------------------------
%
patch=load('patch.mat');
patchlist=patch.patchlist;
npatch=size(patchlist,1);
% process in each patch -------------------------------------------
for k=1:npatch
    patchname=patchlist(k,:);
    if exist(patchname,'dir')
        cd(patchname);
    else
        error('Error, NO such patch: %s \n', patchname);
    end
    fprintf('Processing PATCH: %s \n',patchname);
    fprintf('current time: %s \n',datestr(now));
    % check SHPrecord file-----------------------
    if exist ('SHPrecord','file')
        SHPrecord=stamps_read('SHPrecord', patchsize(k,1)*patchsize(k,2), 'float32');
    end
    % check compressed SLC file-----------------------
    if exist('com_channel1','file')
        comslc.com_channel1=stamps_read('com_channel1', patchsize(k,1)*patchsize(k,2), 'cpxfloat32');
    end
    if exist('com_channel2','file')
        comslc.com_channel2=stamps_read('com_channel2', patchsize(k,1)*patchsize(k,2), 'cpxfloat32');
    end
    if exist('com_channel3','file')
        comslc.com_channel3=stamps_read('com_channel3', patchsize(k,1)*patchsize(k,2), 'cpxfloat32');
    end
    % check Number of images processed -----------------------
    if exist ('processed_num.mat','file')
        load('processed_num.mat', 'processed_num');
    else
        processed_num=0;
    end
    % process in each subset SLCs -------------------------
    temp_ph=[];
    while processed_num<total_num
        no_subset=floor(processed_num/interval)+1;
        disp(['########## progress:' patchname, ': '   ' The ', num2str(no_subset),'th of ' num2str(npart)  ' part ##########']);
        % -- Read Data--------------------------------------------------
        [Mslc,M_MLI]=SAR_data_input(workpath,channels,nlines,processed_num,interval,masterID,[patch_over(k,1), patch_over(k,3), patch_over(k,2), patch_over(k,4)]);
        disp(['Read data: ',patchname, ': ' num2str(no_subset),'th of ' num2str(npart)  ' Done!']);
        %  Homogeneous pixel selection (SHP) using 'HTCI' ---------------
%----with shp update ----------------------------------
        if processed_num+interval > total_num
            load('SHP.mat');
        else
            [SHP]=SHP_SelPoint_HTCI(M_MLI,CalWin,Alpha);
        end
%----NO shp update ----------------------------------
%         if processed_num==0
%             [SHP]=SHP_SelPoint_HTCI(M_MLI(:,:,1:10),CalWin,Alpha);
%         else
%             load('SHP.mat');
%         end
%--------------------------------------

        disp(['SHP SelPoint: ',patchname, ': ' num2str(no_subset),'th of ' num2str(npart)  ' Done!']);
        % phase linking with SETP-EMI -----------------------------------
        %     tic
        if processed_num==0
            [opt_ph,SHPrecord,comslc,~]=SETPEMI_parallel(Mslc,SHP,minshp,hW_l,hW_w);
        else
            if no_subset~=npart
              [opt_ph,SHPrecord,comslc,~]=SETPEMI_parallel(Mslc,SHP,minshp,hW_l,hW_w,comslc,SHPrecord);
            else
              [opt_ph,SHPrecord,comslc,adddatum]=SETPEMI_parallel(Mslc,SHP,minshp,hW_l,hW_w,comslc,SHPrecord,1);
            end
        end
        % toc;
        % Store optimized time series interferograms ----------------------
        temp_ph=[temp_ph,opt_ph];

        % save temp file-----------------------------------------------------
        if size(opt_ph,2)< interval
            add_num=0;
        else
            add_num=interval;
        end
        processed_num=processed_num+add_num;
        stamps_save('processed_num.mat',processed_num);

% If the number of the last datasets is insufficient for SHP identification,
% the SHP identified by the previous subset will be enabled
        stamps_save('SHP.mat',SHP);
        if no_subset==npart-1
          stamps_save('SHP.mat',SHP);
        end
        disp([patchname ': The ' num2str(no_subset) 'th of ' num2str(npart) ' subset Done!']);
        if no_subset==npart
            break
        end
    end
    % Reference phase connection --------------------------------------
    if interval~=total_num
        addph = kron(adddatum, ones(1, interval));
        addph(:,size(temp_ph,2)+1:end)=[];
        opt_phase=temp_ph+addph;
        opt_phase = angle(exp(1j*opt_phase));
    else
        opt_phase=temp_ph;
    end

    % save each patch phase and compressed SLC files -------------------------
    % 
    stamps_write(SHPrecord, 'SHPrecord', 'float32');
    % save phase 
    stamps_write(opt_phase, 'opt_phase', 'float32');
    % save compressed SLC
    switch length(channels)
        case 1
            stamps_write(comslc.com_channel1, 'com_channel1', 'cpxfloat32');
        case 2
            stamps_write(comslc.com_channel1, 'com_channel1', 'cpxfloat32');
            stamps_write(comslc.com_channel2, 'com_channel2', 'cpxfloat32');
        case 3
            stamps_write(comslc.com_channel1, 'com_channel1', 'cpxfloat32');
            stamps_write(comslc.com_channel2, 'com_channel2', 'cpxfloat32');
            stamps_write(comslc.com_channel3, 'com_channel3', 'cpxfloat32');
    end
    % -------------------------
    cd ..;
    processed_num=0;
end
%%



