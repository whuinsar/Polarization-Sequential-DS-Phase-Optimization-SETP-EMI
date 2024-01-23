%% Save final results ----------------------------------------------------------------------------------------------
% stamps_save('SHP.mat',SHP);
% stamps_write(opt_phase, 'opt_diff/opt_phase', 'float32');

%read name  -----------------------------------------------
imgpath=[workpath '/'  char(channels(1)) '/diff/'];
tag_files = dir([imgpath,'*','.diff']);
filename=cell2mat({tag_files.name}');


% save opt diff  -----------------------------------------------
% opt_ph=reshape(opt_phase,cols,rows,size(opt_phase,2));
% opt_ph=permute(opt_ph,[2,1,3]);
patch_ph=(exp(1j*patch_ph));
for i=1: length(filename)
    % stamps_write(opt_phase(:,i), ['opt_diff/' filename(i,:)], 'float32');
    if ~exist(['opt_diff/' filename(i,:)],'file')
        %stamps_write(patch_ph(:,:,i), ['opt_diff/' filename(i,:)],'cpxfloat32');
        fwritebkj(patch_ph(:,:,i), ['opt_diff/' filename(i,:)],'cpxfloat32','b');
    end
end

% save compressed_slc  -----------------------------------------------
logit('Writing compressed_slc into binary files')

if npart*interval~=total_num
    switch length(channels)
        case 1
            stamps_write((comslc.com_channel1(:,1:npart-1)), 'com_slc/com_channel1', 'cpxfloat32');
        case 2
            stamps_write((comslc.com_channel1(:,1:npart-1)), 'com_slc/com_channel1', 'cpxfloat32');
            stamps_write((comslc.com_channel2(:,1:npart-1)), 'com_slc/com_channel2', 'cpxfloat32');
        case 3
            stamps_write((comslc.com_channel1(:,1:npart-1)), 'com_slc/com_channel1', 'cpxfloat32');
            stamps_write((comslc.com_channel2(:,1:npart-1)), 'com_slc/com_channel2', 'cpxfloat32');
            stamps_write((comslc.com_channel3(:,1:npart-1)), 'com_slc/com_channel3', 'cpxfloat32');
    end
    stamps_write(SHPrecord(:,1:npart-1), 'SHPrecord', 'float32');
else
   switch length(channels)
        case 1
            stamps_write(comslc.com_channel1, 'com_slc/com_channel1', 'cpxfloat32');
        case 2
            stamps_write(comslc.com_channel1, 'com_slc/com_channel1', 'cpxfloat32');
            stamps_write(comslc.com_channel2, 'com_slc/com_channel2', 'cpxfloat32');
        case 3
            stamps_write(comslc.com_channel1, 'com_slc/com_channel1', 'cpxfloat32');
            stamps_write(comslc.com_channel2, 'com_slc/com_channel2', 'cpxfloat32');
            stamps_write(comslc.com_channel3, 'com_slc/com_channel3', 'cpxfloat32');
    end

    stamps_write(SHPrecord, 'SHPrecord', 'float32');
end

% -----------------------------------------------