function [opt_ph,SHPrecord,comslc,adddatum]=SETPEMI_parallel(sarin,SHP,minshp,hW_l,hW_w,comslc,SHPrecord,end_flag)
%%-----------------------input---------------------%
%{%}
% sarin :  SLC data set
% SHP : Statistically homogeneous pixel sets
% minshp : Minimum threshold for DS point identification
% hW_l: half window size from range direction of SHP window
% hW_w: half window size from azimuth direction of SHP window
% comslc : Compressed SLC sets
% SHPrecord : SHP number record from historical data
% end_flag: Check if it is the last sequential subset
%%----------------------output---------------------%
%{%}
% opt_ph : Optimized phase
% SHPrecord : SHP number record
% comslc :  Compressed SLC sets
% adddatum : Phase reference for sequential subsets

%   =================================================================
%   02/2023  Yian Wang:  create
%   09/2023  Yian Wang:  modified. Change the data structure to reduce data I/O time
%   10/2023  Yian Wang:  modified. Single or full-polarization data processing supported
%   =================================================================

% function handle: get the remainder
% rem(x,y): sign of remainder is same to dividend(x)
MOD=@(a,b) rem(b+rem(a,b),b)+1;

if nargin < 6
    first_subset=1;
    comslc=[];
else
    first_subset=0;
end
if nargin < 8
    end_flag=0;
end
%% Read the SLC data set and determine the polarization channel
pol_names = char(fieldnames(sarin));
m_channel=length (unique (fieldnames (sarin)));

for i=1:m_channel
    if pol_names(i,:) == ['slcexpVV']
        polvv=sarin.slcexpVV;
        p_M=size(sarin.slcexpVV,1); p_N=size(sarin.slcexpVV,2);  D=size(sarin.slcexpVV,3);
    end
    if pol_names(i,:) == ['slcexpHH']
        polhh=sarin.slcexpHH;
        p_M=size(sarin.slcexpHH,1); p_N=size(sarin.slcexpHH,2);  D=size(sarin.slcexpHH,3);
    end
    if pol_names(i,:) == ['slcexpVH']
        polvh=sarin.slcexpVH;
        p_M=size(sarin.slcexpVH,1); p_N=size(sarin.slcexpVH,2);  D=size(sarin.slcexpVH,3);
    end
    if pol_names(i,:) == ['slcexpHV']
        polhv=sarin.slcexpHV;
        p_M=size(sarin.slcexpHV,1); p_N=size(sarin.slcexpHV,2);  D=size(sarin.slcexpHV,3);
    end
end

%% ----Pauli basis scattering vector generation----------------------
if m_channel==1
    if exist('polvv','var')
        channels=['vv'];
        pol1=sarin.slcexpVV;
    end
    if exist('polhh','var')
        channels=['hh'];
        pol1=sarin.slcexpHH;
    end
end

if m_channel==2
    if exist('polhh','var') && exist('polhv','var')
        channels=['hh+hv'];
        pol1=sarin.slcexpHH;
        pol2=2*(sarin.slcexpHV);
    end
    if exist('polhh','var') && exist('polvh','var')
        channels=['hh+vh'];
        pol1=sarin.slcexpHH;
        pol2=2*(sarin.slcexpVH);
    end
    if exist('polvv','var') && exist('polvh','var')
        channels=['vv+vh'];
        pol1=sarin.slcexpVV;
        pol2=2*(sarin.slcexpVH);
    end
    if exist('polvv','var') && exist('polhh','var')
        channels=['vv+hh'];
        pol1=(sarin.slcexpHH+sarin.slcexpVV)/sqrt(2);
        pol2=(sarin.slcexpHH-sarin.slcexpVV)/sqrt(2);
    end
    if exist('polvv','var') && exist('polhh','var')
        channels=['vv+hh'];
        pol1=(sarin.slcexpHH+sarin.slcexpVV)/sqrt(2);
        pol2=(sarin.slcexpHH-sarin.slcexpVV)/sqrt(2);
    end
end

if m_channel==3
    if exist('polvv','var') && exist('polvh','var') && exist('polhh','var')
        channels=['vv+hh+vh'];
        pol1=(sarin.slcexpVV+sarin.slcexpHH)/sqrt(2);
        pol2=(sarin.slcexpHH-sarin.slcexpVV)/sqrt(2);
        pol3=(sarin.slcexpVH)*sqrt(2);
    end
end

%% Predefine Matrix to reduce time cost
SHP_check0=zeros(p_M*p_N,1);
opt_ph=zeros(p_M*p_N,D);
if first_subset==1
    npart=1;
    com_channel1=zeros(p_M*p_N,(2*hW_l+1)*(2*hW_w+1));
    com_channel2=zeros(p_M*p_N,(2*hW_l+1)*(2*hW_w+1));
    com_channel3=zeros(p_M*p_N,(2*hW_l+1)*(2*hW_w+1));
    temp_com_channel1=zeros(p_M*p_N,npart);
    temp_com_channel2=zeros(p_M*p_N,npart);
    temp_com_channel3=zeros(p_M*p_N,npart);
    com_channel10=com_channel1;
    com_channel20=com_channel2;
    com_channel30=com_channel3;
    adddatum=0;
else
    npart=size(SHPrecord,2)+1;
    com_channel1=[];
    com_channel2=[];
    com_channel3=[];
    if isfield(comslc,'com_channel1')
        com_channel1=reshape(comslc.com_channel1,p_M,p_N,npart-1) ;
        %com_channel1=permute(com_channel1,[2,1,3]);
        temp_com_channel1=zeros(p_M*p_N,npart);
        com_channel10=comslc.com_channel1;
        com_channel20=zeros(size(comslc.com_channel1));
        com_channel30=zeros(size(comslc.com_channel1));
    end
    if isfield(comslc,'com_channel2')
        com_channel2=reshape(comslc.com_channel2,p_M,p_N,npart-1) ;
        %com_channel2=permute(com_channel2,[2,1,3]);
        temp_com_channel2=zeros(p_M*p_N,npart);
        com_channel20=comslc.com_channel2;
        com_channel30=zeros(size(comslc.com_channel2));
    end
    if isfield(comslc,'com_channel3')
        com_channel3=reshape(comslc.com_channel3,p_M,p_N,npart-1) ;
        %com_channel3=permute(com_channel3,[2,1,3]);
        temp_com_channel3=zeros(p_M*p_N,npart);
        com_channel30=comslc.com_channel3;
    end
    adddatum=zeros(p_M*p_N,npart);

end

% mask cohmatrix (test)
% cohmask=1;

%%
all_step=p_M*p_N;
p=1;num=1;
%% process on each pixels
%temp_mslc = cell(p_M, p_N);
%idd=1;
SHP0=SHP;
clear SHP;
pol10=0;
pol20=0;
pol30=0;
%comslc0=comslc;
%clear comslc;

switch m_channel
    case 1
        pol10=pol1;clear pol1;
    case 2
        pol10=pol1;clear pol1;
        pol20=pol2;clear pol2;
    case 3
        pol10=pol1;clear pol1;
        pol20=pol2;clear pol2;
        pol30=pol3;clear pol3;
end

%tic;
coh_all=0;
coh_count=0;
% parfor idd=1:all_step
for idd=1:all_step
    %comslc=comslc0;
    pol1=0;
    pol2=0;
    pol3=0;
    SHP=SHP0;
    p_N0=p_N;
    p_M0=p_M;
    switch m_channel
        case 1
            pol1=pol10;
        case 2
            pol1=pol10;
            pol2=pol20;
        case 3
            pol1=pol10;
            pol2=pol20;
            pol3=pol30;
    end
    % -------------------------------------------------------
    %     j=ceil(idd/p_M);
    %     i=mod(idd,p_M);
    %     if i==0
    %         i=p_M;
    %     end
    %     if i==100 && j==40   
    %         a=1;
    %     end
    % extract SHP ---------------------------------------------------
    [SHP_check,SHPpatch,id_center,ix_l,ix_p]=extract_shp(SHP,idd,hW_l,hW_w,p_N0,p_M0,minshp);
    % extract SHP in SLCs---------------------------------------
    switch m_channel
        case 1
            sarblock1 = extract_shp_slc( pol1,ix_l,ix_p,SHPpatch,hW_l,hW_w,D);
        case 2
            sarblock1 = extract_shp_slc( pol1,ix_l,ix_p,SHPpatch,hW_l,hW_w,D);
            sarblock2 = extract_shp_slc( pol2,ix_l,ix_p,SHPpatch,hW_l,hW_w,D);
        case 3
            sarblock1 = extract_shp_slc( pol1,ix_l,ix_p,SHPpatch,hW_l,hW_w,D);
            sarblock2 = extract_shp_slc( pol2,ix_l,ix_p,SHPpatch,hW_l,hW_w,D);
            sarblock3 = extract_shp_slc( pol3,ix_l,ix_p,SHPpatch,hW_l,hW_w,D);
    end
    %-----------------------------
%         if  size(sarblock1,2) < minshp
%             continue
%         end
    test_nonan = size(sarblock1,2);
    if test_nonan<id_center
        id_center=test_nonan;
    end
    %-----------------------------
    if (first_subset==0)
        switch m_channel
            case 1
                cslc1 = extract_shp_slc( com_channel1,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart-1);
            case 2
                cslc1 = extract_shp_slc( com_channel1,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart-1);
                cslc2 = extract_shp_slc( com_channel2,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart-1);
            case 3
                cslc1 = extract_shp_slc( com_channel1,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart-1);
                cslc2 = extract_shp_slc( com_channel2,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart-1);
                cslc3 = extract_shp_slc( com_channel3,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart-1);
        end
    end
    %-----------------------------

    %% phase linking with TP-EMI
    if (first_subset)
        % for the first sequential SLC subset
        switch m_channel
            case 1
                % Generate coherence matrix
                cohmat1 = CovEst(sarblock1,SHP_check);
                % coherence matrix Stacking and EMI phase optimization
                [slcopt,phtemp] = optph_emi(cohmat1,SHP_check);
                opt_ph(idd,:)=phtemp';
                %Linear transformation for compresses original SLC
                slcopt1=slcopt;
                slcopt1 = slcopt1/norm(slcopt1);
                tempPCA = slcopt1'*sarblock1;
                temp_com_channel1(idd,:) = tempPCA(id_center);

            case 2
                % Generate coherence matrix
                cohmat1 = CovEst(sarblock1,SHP_check);
                cohmat2 = CovEst(sarblock2,SHP_check);
                % coherence matrix Stacking and EMI phase optimization
                [slcopt,phtemp] = optph_emi(cohmat1+cohmat2,2*SHP_check);
                opt_ph(idd,:)=phtemp';
                %Linear transformation for compresses original SLC
                slcopt1=slcopt;
                slcopt1 = slcopt1/norm(slcopt1);
                tempPCA = slcopt1'*sarblock1;
                tempPCA2 = slcopt1'*sarblock2;
                temp_com_channel1(idd,:) = tempPCA(id_center);
                temp_com_channel2(idd,:) = tempPCA2(id_center);     
                
            case 3
                % Generate coherence matrix
                cohmat1 = CovEst(sarblock1,SHP_check);
                cohmat2 = CovEst(sarblock2,SHP_check);
                cohmat3 = CovEst(sarblock3,SHP_check);
                % coherence matrix Stacking and EMI phase optimization
                [slcopt,phtemp] = optph_emi(cohmat1+cohmat2+cohmat3,3*SHP_check);
                opt_ph(idd,:)=phtemp';
                %Linear transformation for compresses original SLC
                slcopt1=slcopt;
                slcopt1 = slcopt1/norm(slcopt1);
                tempPCA = slcopt1'*sarblock1;
                tempPCA2 = slcopt1'*sarblock2;
                tempPCA3 = slcopt1'*sarblock3;

                temp_com_channel1(idd,:) = tempPCA(id_center);
                temp_com_channel2(idd,:) = tempPCA2(id_center);
                temp_com_channel3(idd,:) = tempPCA3(id_center);

        end
    else
        % for the sequential subset slc
        switch m_channel
            case 1
                Initslc1 = vertcat(cslc1,sarblock1);
                % Generate coherence matrix ----------------------------
                cohmat1 = CovEst(Initslc1,SHP_check);
                % Coherence matrix Stacking and EMI phase optimization-
                [slcopt,phtemp] = optph_emi(cohmat1,SHP_check);
                % Calculate the optimal phase of a subset --------------
                phDec = phtemp - phtemp(npart);
                opt_ph(idd,:) = phDec(npart:length(slcopt));
                % Linear transformation for compresses original SLC ----
                slcopt1=slcopt.*conj(slcopt(npart));
                slcopt1=slcopt1(npart:length(slcopt));
                slcopt1 = slcopt1/norm(slcopt1);
                tempPCA = slcopt1'*sarblock1;
                temp_com_channel1(idd,:)= [com_channel10(idd,:) ,tempPCA(id_center)];

            case 2
                Initslc1 = vertcat(cslc1,sarblock1);
                Initslc2 = vertcat(cslc2,sarblock2);
                % Generate coherence matrix ----------------------------
                cohmat1 = CovEst(Initslc1,SHP_check);
                cohmat2 = CovEst(Initslc2,SHP_check);
                % Coherence matrix Stacking and EMI phase optimization-
                [slcopt,phtemp] = optph_emi(cohmat1+cohmat2,2*SHP_check);
                % Calculate the optimal phase of a subset ------------
                phDec = phtemp - phtemp(npart);
                opt_ph(idd,:) = phDec(npart:length(slcopt));
                % Linear transformation for compresses original SLC ---
                slcopt1=slcopt.*conj(slcopt(npart));
                slcopt1=slcopt1(npart:length(slcopt));
                slcopt1 = slcopt1/norm(slcopt1);
                tempPCA = slcopt1'*sarblock1;
                tempPCA2 = slcopt1'*sarblock2;
                temp_com_channel1(idd,:)= [com_channel10(idd,:) ,tempPCA(id_center)];
                temp_com_channel2(idd,:)= [com_channel20(idd,:) ,tempPCA2(id_center)];
            case 3
                Initslc1 = vertcat(cslc1,sarblock1);
                Initslc2 = vertcat(cslc2,sarblock2);
                Initslc3 = vertcat(cslc3,sarblock3);
                % Generate coherence matrix ------------------------
                cohmat1 = CovEst(Initslc1,SHP_check);
                cohmat2 = CovEst(Initslc2,SHP_check);
                cohmat3 = CovEst(Initslc3,SHP_check);
                % Coherence matrix Stacking and EMI phase optimization-
                [slcopt,phtemp] = optph_emi(cohmat1+cohmat2+cohmat3,3*SHP_check);
                % Calculate the optimal phase of a subset ------------
                phDec = phtemp - phtemp(npart);
                opt_ph(idd,:) = phDec(npart:length(slcopt));
                % Linear transformation for compresses original SLC ---
                slcopt1=slcopt.*conj(slcopt(npart));
                slcopt1=slcopt1(npart:length(slcopt));
                slcopt1 = slcopt1/norm(slcopt1);
                tempPCA = slcopt1'*sarblock1;
                tempPCA2 = slcopt1'*sarblock2;
                tempPCA3 = slcopt1'*sarblock3;

                temp_com_channel1(idd,:)= [com_channel10(idd,:) ,tempPCA(id_center)];
                temp_com_channel2(idd,:)= [com_channel20(idd,:) ,tempPCA2(id_center)];
                temp_com_channel3(idd,:)= [com_channel30(idd,:) ,tempPCA3(id_center)];

        end
        SHP_check0(idd)=SHP_check;

    end

    %% Triangular Phase Closure Coherence Calculation
    % gamma(i,j)=estgamma(reshape(optph(i,j,:),D,1),cohmat);   %temporal coherence
    % meancoh(i,j)=sum(sum(abs(triu(cohmat,1)),1),2)/(D*(D-1)/2); % mean coherence
    % counting ----------------------------------------------------------
    %         idd=idd+1;
    %             num=num+1;
    %             if num == floor(all_step * 0.1*p);
    %                 disp(['progress: ', num2str(1*10*p),'%']);
    %                 p = p+1;
    %             end

end

%% Link to phase reference ----------------------------

if end_flag ==1
    switch m_channel
        case 1
            link_slc1=reshape(temp_com_channel1,p_M,p_N,npart) ;
            %link_slc1=permute(link_slc1,[2,1,3]);
            link_slc2=link_slc1;
            link_slc3=link_slc1;
        case 2
            link_slc1=reshape(temp_com_channel1,p_M,p_N,npart) ;
            %link_slc1=permute(link_slc1,[2,1,3]);
            link_slc2=reshape(temp_com_channel2,p_M,p_N,npart) ;
            %link_slc2=permute(link_slc2,[2,1,3]);
            link_slc3=link_slc1;
        case 3
            link_slc1=reshape(temp_com_channel1,p_M,p_N,npart) ;
            %link_slc1=permute(link_slc1,[2,1,3]);
            link_slc2=reshape(temp_com_channel2,p_M,p_N,npart) ;
            %link_slc2=permute(link_slc2,[2,1,3]);
            link_slc3=reshape(temp_com_channel3,p_M,p_N,npart) ;
            %link_slc3=permute(link_slc3,[2,1,3]);
    end
    parfor idd=1:all_step
   %  for idd=1:all_step
        SHP=SHP0;
        % -------------------------------------------------------
        % extract SHP ---------------------------------------------
        [SHP_check,SHPpatch,~,ix_l,ix_p]=extract_shp(SHP,idd,hW_l,hW_w,p_N,p_M,minshp);
        % extract SHP in SLCs---------------------------------------
        switch m_channel
            case 1
                sarblock1 = extract_shp_slc( link_slc1,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart);
                cohmat1 = CovEst(sarblock1,SHP_check);
                [~,phtemp] = optph_emi(cohmat1,SHP_check);
                adddatum(idd,:)=phtemp';
            case 2
                sarblock1 = extract_shp_slc( link_slc1,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart);
                sarblock2 = extract_shp_slc( link_slc2,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart);
                % Generate coherence matrix
                cohmat1 = CovEst(sarblock1,SHP_check);
                cohmat2 = CovEst(sarblock2,SHP_check);
                % coherence matrix Stacking and EMI phase optimization
                [~,phtemp] = optph_emi(cohmat1+cohmat2,2*SHP_check);
                adddatum(idd,:)=phtemp';
            case 3
                sarblock1 = extract_shp_slc( link_slc1,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart);
                sarblock2 = extract_shp_slc( link_slc2,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart);
                sarblock3 = extract_shp_slc( link_slc3,ix_l,ix_p,SHPpatch,hW_l,hW_w,npart);
                % Generate coherence matrix
                cohmat1 = CovEst(sarblock1,SHP_check);
                cohmat2 = CovEst(sarblock2,SHP_check);
                cohmat3 = CovEst(sarblock3,SHP_check);
                % coherence matrix Stacking and EMI phase optimization
                [~,phtemp] = optph_emi(cohmat1+cohmat2+cohmat3,3*SHP_check);
                adddatum(idd,:)=phtemp';
        end
        %----------------------
    end
end
%%
%toc;
switch m_channel
    case 1
        com_channel1=temp_com_channel1;
        comslc= struct('com_channel1', com_channel1);
    case 2
        com_channel1=temp_com_channel1;
        com_channel2=temp_com_channel2;
        comslc= struct('com_channel1', com_channel1, 'com_channel2', com_channel2);
    case 3
        com_channel1=temp_com_channel1;
        com_channel2=temp_com_channel2;
        com_channel3=temp_com_channel3;
        comslc= struct('com_channel1', com_channel1, 'com_channel2', com_channel2,'com_channel3', com_channel3);
end

if npart==1
    SHPrecord=SHP_check0;
else
    SHPrecord=[SHPrecord,SHP_check0];
end

if(1)
    temppp=abs(coh_all)/(coh_count);
    save('coh.mat','temppp');
end
% close(h1);
end

function covMT = CovEst( slcstack,num_shp)
%tempslc = exp(1j*angle(slcstack));
covMT = slcstack*slcstack';
covMT = covMT/num_shp;
abscov=sqrt(diag(covMT)*diag(covMT)');
abscov(abscov==0)=1;
covMT=covMT./abscov;
end
function sarblock = extract_shp_slc( slcs,ix_l,ix_p,SHPpatch,hW_l,hW_w,D)
sarblock=slcs(ix_l(:),ix_p(:),:);
sarblock=SHPpatch.*sarblock;
sarblock=reshape(sarblock,(2*hW_l+1)*(2*hW_w+1),D);
nanRows = any(isnan(sarblock), 2);sarblock = sarblock(~nanRows, :);
sarblock=sarblock.';
end


