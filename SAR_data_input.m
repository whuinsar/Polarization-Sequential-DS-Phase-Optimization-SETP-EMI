function [M_SLC,M_MLI]=SAR_data_input(workpath,channels,nlines,processed_num,interval,masterID,Blk)
%% introduction ---------------------------------
% This is a program that reads multipolarized (or single-polarized) SLC,
% differential interferograms, and intensity maps. The output is the SLC
% data set with topographic phase removed and the time series intensity map.
%----------------------------------------------------------
%% usage ---------------------------------
% The "workpath" should contain a folder named after the polarization channel,
% and each subfolder should store three folders: diff, slc_f, and rmli.
% 1.The diff folder stores time series differential interferograms based on
% a single master image, with the suffix '.diff'
% 2.The slc folder stores the time series SLC, with the suffix '.slc_f'
% 3.The rmli folder stores time series intensity images, with the suffix '.rmli'

% channels : "vv" or "hh" or "vh" in a string matrix
% nlines : number of lines in slc (The size of SLC, can be found in slc.par file)
% processed_num : Number of images processed,the start number to read data
% interval : Number of images need to read.
% masterID : the ID of master image
% Blk : [r0, rN, c0, cN] : Crop range

%% example --------------------------------
% [M_SLC,M_MLI]=SAR_data_input(workpath,["vv"; "vh"],901,0,10,58,[0, 300, 0, 1000]);

%% Revision history
% created by Yian Wang, 10.10.2023

%% processing---------------------------------------
if (sum(ismember(channels, 'vv')))
    %% for vv
    % Read mli vv data stepwise--------------------
    temp=ImgRead_stepwise([workpath '/vv/rmli'] ,'.rmli',nlines,processed_num,interval,'float32','b',Blk);
    mlivv=temp.datastack;
    %-------------------------------------
    % Read diff vv data stepwise -------------------
    temp=ImgRead_stepwise([workpath '/vv/diff'] ,'.diff',nlines,processed_num,interval,'cpxfloat32','b',Blk);
    cint_minrefdemvv=temp.datastack;
    % Read slc vv data stepwise --------------------
    temp=ImgRead_stepwise([workpath '/vv/slc'] ,'.rslc_f',nlines,processed_num,interval,'cpxfloat32','b',Blk);
    sarcpxvv=temp.datastack;
    % Read master ID data --------------------
    % vv diff master
    temp=ImgRead_one([workpath '/vv/diff'] ,['.diff'],nlines,'cpxfloat32','b',Blk,masterID);
    master_diff=temp.datastack;
    % vv slc master
    temp=ImgRead_one([workpath '/vv/slc'] ,['.rslc_f'],nlines,'cpxfloat32','b',Blk,masterID);
    master_slc=temp.datastack;
    % Calculate SLC with terrain phase removed --------------------
    cintvv= master_slc.*conj(sarcpxvv);
    if processed_num+interval>=masterID & processed_num < masterID
        cint_minrefdemvv(:,:,masterID-processed_num)=master_diff.*0+(1+0i);
    end
    cintvv=exp(1i.*angle(cintvv));
    sarcpx_minrefvv=sarcpxvv.*(cintvv.*conj(cint_minrefdemvv));
    clear cint_minrefdemvv sarcpxvv cintvv
end
if (sum(ismember(channels, 'vh')))
    %% for vh
    % Read mli vh data stepwise--------------------
    temp=ImgRead_stepwise([workpath '/vh/rmli'] ,'.rmli',nlines,processed_num,interval,'float32','b',Blk);
    mlivh=temp.datastack;
    % Read diff vh data stepwise ----------------------------------
    temp=ImgRead_stepwise([workpath '/vh/diff'] ,'.diff',nlines,processed_num,interval,'cpxfloat32','b',Blk);
    cint_minrefdemvh=temp.datastack;
    % Read slc vh data stepwise ------------------------------------
    temp=ImgRead_stepwise([workpath '/vh/slc'] ,'.rslc_f',nlines,processed_num,interval,'cpxfloat32','b',Blk);
    sarcpxvh=temp.datastack;
    % Read master ID data  -----------------------------------------
    % vv diff master
    temp=ImgRead_one([workpath '/vh/diff'] ,['.diff'],nlines,'cpxfloat32','b',Blk,masterID);
    master_diff=temp.datastack;
    % vv slc master
    temp=ImgRead_one([workpath '/vh/slc'] ,['.rslc_f'],nlines,'cpxfloat32','b',Blk,masterID);
    master_slc=temp.datastack;

    % Calculate SLC with terrain phase removed --------------------
    cintvh= master_slc.*conj(sarcpxvh);
    if processed_num+interval>=masterID & processed_num<masterID
        cint_minrefdemvh(:,:,masterID-processed_num)=master_diff.*0+(1+0i);
    end

    cintvh=exp(1i.*angle(cintvh));
    sarcpx_minrefvh=sarcpxvh.*(cintvh.*conj(cint_minrefdemvh));
    clear cintvh sarcpxvh cint_minrefdemvh
end
if (sum(ismember(channels, 'hh')))
        %% for hh
        % Read mli hh data stepwise--------------------
    temp=ImgRead_stepwise([workpath '/hh/rmli'] ,'.rmli',nlines,processed_num,interval,'float32','b',Blk);
    mlihh=temp.datastack;
    % Read diff hh data stepwise -------------------
    temp=ImgRead_stepwise([workpath '/hh/diff'] ,'.diff',nlines,processed_num,interval,'cpxfloat32','b',Blk);
    cint_minrefdemhh=temp.datastack;
    % Read slc hh data stepwise --------------------
    temp=ImgRead_stepwise([workpath '/hh/slc'] ,'.rslc_f',nlines,processed_num,interval,'cpxfloat32','b',Blk);
    sarcpxhh=temp.datastack;
    % Read master ID data  --------------------
    % hh diff master
    temp=ImgRead_one([workpath '/hh/diff'] ,['.diff'],nlines,'cpxfloat32','b',Blk,masterID);
    master_diff=temp.datastack;
    % hh slc master
    temp=ImgRead_one([workpath '/hh/slc'] ,['.rslc_f'],nlines,'cpxfloat32','b',Blk,masterID);
    master_slc=temp.datastack;

    % Calculate SLC with terrain phase removed --------------------
    cinthh= master_slc.*conj(sarcpxhh);
    if processed_num+interval>=masterID & processed_num<masterID
        cint_minrefdemhh(:,:,masterID-processed_num)=master_diff.*0+(1+0i);
    end

    cinthh=exp(1i.*angle(cinthh));
    sarcpx_minrefhh=sarcpxhh.*(cinthh.*conj(cint_minrefdemhh));
    clear cinthh sarcpxhh cint_minrefdemhh
end
if exist('sarcpx_minrefvv') 
    M_SLC.slcexpVV = sarcpx_minrefvv;
end
if  exist('sarcpx_minrefhh')
    M_SLC.slcexpHH = sarcpx_minrefhh;
end
if  exist('sarcpx_minrefvh')
    M_SLC.slcexpVH = sarcpx_minrefvh;
end
if exist('mlivv','var')
    M_MLI=mlivv;
else 
    if exist('mlihh','var')
        M_MLI=mlihh;
    else 
        if exist('mlivv','var')
        M_MLI=mlivh;
        end
    end
end

end


