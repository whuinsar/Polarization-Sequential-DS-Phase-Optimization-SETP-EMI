%% SETP-EMI procesing ####################################################
patch=load('patch.mat');
patchlist=patch.patchlist;
npatch=size(patchlist,1);
%% process in each Spatial patch -------------------------------------------
for k=1:npatch
%    for k=2
    patchname=patchlist(k,:);
    if exist(patchname,'dir')
        cd(patchname);
    else
        error('Error, NO such patch: %s \n', patchname);
    end
    fprintf('Processing PATCH: %s \n',patchname);
    fprintf('current time: %s \n',datestr(now));
    % check SHPrecord file------------------------------------------------
    if exist ('./SHPrecord','file')
        SHPrecord=stamps_read('SHPrecord', patchsize(k,1)*patchsize(k,2), 'float32');
    else
        if ~exist ('SHPrecord','var')
            SHPrecord=[];
        end
    end
    % check compressed SLC file-------------------------------------------
    if exist('./com_channel1','file')
        comslc.com_channel1=stamps_read('com_channel1', patchsize(k,1)*patchsize(k,2), 'cpxfloat32');
        num_cslc=size(comslc.com_channel1,2);
    end
    if exist('./com_channel2','file')
        comslc.com_channel2=stamps_read('com_channel2', patchsize(k,1)*patchsize(k,2), 'cpxfloat32');
        num_cslc=size(comslc.com_channel2,2);
    end
    if exist('./com_channel3','file')
        comslc.com_channel3=stamps_read('com_channel3', patchsize(k,1)*patchsize(k,2), 'cpxfloat32');
        num_cslc=size(comslc.com_channel3,2);
    end
    if ~exist('comslc','var')
        comslc=[];
        num_cslc=0;
    end
    % check Number of images processed -----------------------------------
    if exist ('./processed_num.mat','file')
        load('processed_num.mat', 'processed_num');
        oriprocessed_num=processed_num;
    else
        oriprocessed_num=0;
        processed_num=0;
        comslc=[];
        num_cslc=0;
    end
    % check Number of compressed SLC with processed_num-------------------
    if num_cslc>(processed_num/interval)
        SHPrecord(:,end)=[];
        if exist('./com_channel1','file')
            comslc.com_channel1(:,end)=[];
        end
        if exist('./com_channel2','file')
            comslc.com_channel2(:,end)=[];
        end
        if exist('./com_channel3','file')
            comslc.com_channel3(:,end)=[];
        end
    end
    %% process in each mini time seires SLC subset  ----------------------
    temp_ph=[];
    while processed_num<total_num
        no_subset=floor(processed_num/interval)+1;
        disp(['########## progress:' patchname, ': '   ' The ', ...
            num2str(no_subset),'th of ' num2str(npart)  ' part ##########']);
        % -- Read Data--------------------------------------------------
        [Mslc,M_MLI]=SAR_data_input(workpath,channels,nlines,processed_num,...
            interval,masterID,[patch_over(k,1), patch_over(k,3), patch_over(k,2), patch_over(k,4)]);
        disp(['Read data: ',patchname, ': ' num2str(no_subset),'th of ' num2str(npart)  ' Done!']);
        % ---------------------------------------------------------------
        %  Homogeneous pixel selection (SHP) using 'HTCI' ---------------
        %----with shp update 
        if processed_num+interval > total_num
            load('SHP.mat');
        else
            % using FaSHPs identification -------------------------------
            [SHP]=SHP_SelPoint_HTCI(M_MLI,CalWin,Alpha,'HTCI');
        end
        disp(['SHP SelPoint: ',patchname, ': ' num2str(no_subset),'th of ' num2str(npart)  ' Done!']);
        % phase linking with SETP-EMI -----------------------------------
        if processed_num==0
            [opt_ph,SHPrecord,comslc]=SETPEMI_parallel(Mslc,SHP,minshp,hW_l,hW_w,comslc,SHPrecord);
        else
            if no_subset~=npart
                [opt_ph,SHPrecord,comslc,~]=SETPEMI_parallel(Mslc,SHP,minshp,hW_l,hW_w,comslc,SHPrecord,0);
            else
                [opt_ph,SHPrecord,comslc,adddatum]=SETPEMI_parallel(Mslc,SHP,minshp,hW_l,hW_w,comslc,SHPrecord,1);
            end
        end
        disp(['phase linking with SETP-EMI: ',patchname, ': ' num2str(no_subset),'th of ' num2str(npart)  ' Done!']);
        % Store optimized time series interferograms ----------------------
        temp_ph=[temp_ph,opt_ph];

        % save temp file-----------------------------------------------------
        if size(opt_ph,2)< interval
            add_num=0;
        else
            add_num=interval;
        end
        processed_num=processed_num+add_num;
        % If the number of the last datasets is insufficient for SHP identification,
        % the SHP identified by the previous subset will be enabled
        % stamps_save('SHP.mat',SHP);
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
        %addph(:,size(temp_ph,2)+1:end)=[];
        addph=addph(:,oriprocessed_num+1:oriprocessed_num+size(temp_ph,2));
        opt_phase=angle(exp(1j*temp_ph).*exp(1j*addph));
        %opt_phase = angle(exp(1j*opt_phase));
    else
        opt_phase=temp_ph;
    end

    % save each patch phase and compressed SLC files -------------------------
    stamps_write(SHPrecord, 'SHPrecord', 'float32');
    stamps_write(adddatum, 'adddatum', 'float32');
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

    stamps_save('processed_num.mat',processed_num);
    stamps_save('oriprocessed_num.mat',oriprocessed_num);
    % -------------------------
    cd ..;
    processed_num=0;
end
%%



