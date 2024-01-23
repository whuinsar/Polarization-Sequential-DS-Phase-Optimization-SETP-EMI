function [patchlist, patch_noover_local, patch_over_local, patch_noover, patch_over,patchsize]=generatepatch(M_start, N_start, M, N, paz, prg, overlap_az, overlap_rg)
%function [patchlist, patch_noover_local, patch_over_local, patch_noover, patch_over,patchsize]=generatepatch(M_start, N_start, M, N, paz, prg, overlap_az, overlap_rg)
% Usage: generate several patches for patch prcossing, output the start and
% end pixel indexes for each patch and make corresponding patch folders
%%--------input--------%
% M_start, N_start:  topleft pixel of being extracted image
% M, N: size of being extracted image
% paz, prg: patches along azimuth and range
% overlap_az, overlap_rg: overlap between adjacent patches
%%-------output-------%
% patchlist(paz*prg, :): list of all pathes, patch along row first
% patch_noover_local(paz*prg, 4): nooverlap local pathes, each row: start_az1,end_az1,start_rg1,end_rg1
% patch_over_local(paz*prg, 4): nooverlap local pathes, each row: start_az2,end_az2,start_rg2,end_rg2
% patch_noover(paz*prg, 4): nooverlap pathes with pixel coords in the whole image , each row: start_az3,end_az3,start_rg3,end_rg3
% patch_over(paz*prg, 4): overlap pathes, each row: start_az4,end_az4,start_rg4,end_rg4
% patchsize(paz*prg, 2): size for each patch
%
%   =================================================================
%   18-Aug-2015 17:18:53   DJ:  create
%   20-Nov-2015 17:11:27   DJ:  modify
%   =================================================================

% size of each patch
length_p=ceil((M-M_start+1)/paz); width_p=ceil((N-N_start+1)/prg);   %????

irg=0;   iaz=0;  ip=0;
patch_noover_local=zeros(paz*prg, 4); % each row: start_az1,start_rg1,end_az1,end_rg1
patch_over_local=zeros(paz*prg, 4); % each row: start_az2,start_rg2,end_az2,end_rg2   added by DJ 20151120
patch_noover=zeros(paz*prg, 4); % each row: start_az3,start_rg3,end_az3,end_rg3
patch_over=zeros(paz*prg, 4); % each row: start_az4,start_rg4,end_az4,end_rg4
patchsize=zeros(paz*prg, 2);
while iaz < paz
    iaz=iaz+1;
    while irg < prg
        irg=irg+1;
        ip=ip+1;
        
        % azimuth 
        start_az1=length_p*(iaz-1)+1;  % no overlap 
        end_az1=length_p*iaz;
        start_az2=start_az1-overlap_az; % overlap 
        end_az2=end_az1+overlap_az;
        if start_az2 < 1    start_az2=1;   end
        if end_az2 > M   end_az2=M;   end
        if start_az1 < 1    start_az1=1;   end
        if end_az1 > M   end_az1=M;   end        
        
        % range
        start_rg1=width_p*(irg-1)+1;
        end_rg1=width_p*irg;
        start_rg2=start_rg1-overlap_rg; 
        end_rg2=end_rg1+overlap_rg;
        if start_rg2 <1    start_rg2=1;  end
        if end_rg2 > N    end_rg2=N;   end
        if start_rg1 <1    start_rg1=1;  end
        if end_rg1 > N    end_rg1=N;   end
        
        % convert to pixel index in the whole image
        start_az3=start_az1+(M_start-1); 
        end_az3=end_az1+(M_start-1);
        start_az4=start_az2+(M_start-1); 
        end_az4=end_az2+(M_start-1);
        
        start_rg3=start_rg1+(N_start-1); 
        end_rg3=end_rg1+(N_start-1);
        start_rg4=start_rg2+(N_start-1); 
        end_rg4=end_rg2+(N_start-1);
        
        
        patch_noover_local(ip,:)=[start_az1,start_rg1,end_az1,end_rg1];
        patch_over_local(ip,:)=[start_az2,start_rg2,end_az2,end_rg2];
        patch_noover(ip,:)=[start_az3,start_rg3,end_az3,end_rg3];
        patch_over(ip,:)=[start_az4,start_rg4,end_az4,end_rg4];
        
        patchsize(ip,:)=[end_az2-start_az2+1,end_rg2-start_rg2+1];
           
        if ip < 10
            ip_str=['00',num2str(ip)];
        elseif ip <100  & ip >= 10
            ip_str=['0',num2str(ip)];
        else 
            ip_str=num2str(ip); 
        end
        
        patchname=strcat('PATCH_',num2str(ip_str));
        patchlist(ip,:)=patchname;
        if ~exist(patchname,'dir')
            mkdir(patchname);
        end
    end
    irg=0;
end
