%% Save final results ----------------------------------------------------------------------------------------------
% stamps_save('SHP.mat',SHP);
% stamps_write(opt_phase, 'opt_diff/opt_phase', 'float32');

%read name  -----------------------------------------------
imgpath=[workpath '/'  char(channels(1)) '/diff/'];
tag_files = dir([imgpath,'*','.diff']);
filename=cell2mat({tag_files.name}');

% check Number of images processed -----------------------
if exist ('PATCH_001/oriprocessed_num.mat','file')
    load('PATCH_001/oriprocessed_num.mat', 'oriprocessed_num');
else
    oriprocessed_num=0;
end
filename=filename(oriprocessed_num+1:end,:);
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
% save SHPrecord, gamma and meancoh --------------------------------
% stamps_write(patch_SHPrecord, 'SHPrecord', 'float32');
% stamps_write(patch_gamma, 'gamma', 'float32');
% stamps_write(patch_meancoh, 'meancoh', 'float32');



% -----------------------------------------------