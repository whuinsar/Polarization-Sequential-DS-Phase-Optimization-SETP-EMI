function [SHP]=SHP_SelPoint_HTCI(mlistack,CalWin,Alpha,EstAgr,prefilter)
%This function is used to select homogeneous pixels (SHP) on SAR intensity stack with single look 
%   Usage:
%       [SHP]=SHP_SelPoint(mlistack,CalWin,Alpha,EstAgr)
%   Inputs:
%   - mlistack: A height by width by page matrix
%   - CalWin:   Fixed boxcar window size
%   - Alpha:    A value between 0 and 1 specifying the
%               significance level. Default is 0.05 for 5% significance.
%   - EstAgr:   Pixel selection algorithm: 1) HTCI; 2) BWS (support Alpha =.05 and .01 only)
%   Outputs:
%   - SHP.PixelInd: A CalWin(1)*CalWin(2) by size(mlistack,1)*size(mlistack,2) array with elements of type logical, containing a SHPs set per pixel 
%   - SHP.BroNum:   The SHP number per pixel (reference pixel is not included) 
%   - SHP.CalWin:   Fixed boxcar window size
if nargin < 5
    prefilter = 0;
end
guass_win=2;
guass_sigma=1;

if nargin < 4
    EstAgr='HTCI';
end

if nargin < 3
    Alpha = 0.05;
end

if nargin < 2
    CalWin = [15 15];
end

if nargin < 1
    help SHP_SelPoint
    return;
end

%tic;

if length(size(mlistack))~=3
    error('Please input 3D matrix...');
end

[nlines,nwidths,npages] = size(mlistack);
mlistack=single(mlistack);

% guass prefilter for slc
if prefilter == 1
    mlistack=guass_prefilter(mlistack,guass_win,guass_sigma);
end
%

%Parameter prepare:
RadiusRow=(CalWin(1)-1)/2;
RadiusCol=(CalWin(2)-1)/2;
InitRow=(CalWin(1)+1)/2; % InitRow is CenterRow
InitCol=(CalWin(2)+1)/2; % InitCol is CenterCol
LRT_nl = 3;
LRT_nw = 3;
if RadiusRow<LRT_nl
    LRT_nl=1;
end
if RadiusCol<LRT_nw
    LRT_nw=1;
end
% Select pixels randomly
% LRT_index = randsample(CalWin(1)*CalWin(1),(2*LRT_nl+1)*(2*LRT_nw+1)); 
[line,pixel] = size(mlistack(:,:,1));
%Statistical threshold:
CR_lo = finv(Alpha/2,2*npages,2*npages);
CR_up = finv(1-Alpha/2,2*npages,2*npages);
Galpha_L = gaminv(Alpha/2,npages,1);
Galpha_U = gaminv(1-Alpha/2,npages,1);
%Edeg mirror-image
mlistack = padarray(mlistack,[RadiusRow RadiusCol],'symmetric');
[pixelind,lineind] = meshgrid(1:pixel,1:line);
lineindpad = padarray(lineind,[RadiusRow RadiusCol],'symmetric');
pixelindpad = padarray(pixelind,[RadiusRow RadiusCol],'symmetric');
indmatrix = sub2ind([line,pixel],lineindpad,pixelindpad);
meanmli = mean(mlistack,3);
[nlines_EP,nwidths_EP]= size(meanmli);
% SHP.PixelInd=false(CalWin(1)*CalWin(2),nlines*nwidths);
SHP.PixelSub=zeros(CalWin(1)*CalWin(2),nlines*nwidths);

%estimate SHPs
num=1;
p=1;
all = nlines*nwidths;
all_step = floor(all/100);

if strcmpi(EstAgr,'HTCI') 
    for kk=InitCol:nwidths_EP-RadiusCol
        for ll=InitRow:nlines_EP-RadiusRow       
            %Initial estimation (Likelihood-ratio test)
            temp = meanmli(ll-LRT_nl:ll+LRT_nl,kk-LRT_nw:kk+LRT_nw);
            tempind = indmatrix(ll-RadiusRow:ll+RadiusRow,kk-RadiusCol:kk+RadiusCol);
            T = meanmli(ll,kk)./temp;
            T = T>CR_lo&T<CR_up;
            SeedPoint = mean(temp(T));
            %iteration (Gamma Confidence interval)
            MeanMatrix = meanmli(ll-RadiusRow:ll+RadiusRow,kk-RadiusCol:kk+RadiusCol);
            SeedPoint = MeanMatrix>Galpha_L*SeedPoint/npages&MeanMatrix<Galpha_U*SeedPoint/npages; %check membership
            SeedPoint(InitRow,InitCol)=true;
            %connection
            LL = bwlabel(SeedPoint); %double
            flag = LL(:)==LL(InitRow,InitCol);
            SHP.PixelSub(:,num)= flag.*tempind(:);  
            num=num+1;
%             if num == all_step * p;
%                 disp(['progress: ', num2str(1*p),'%']);
%                 p = p+1;
%             end
        end
    end
else %Baumgartner-weiB-Schindler test (How about KS test?)
    for kk=InitCol:nwidths_EP-RadiusCol
        for ll=InitRow:nlines_EP-RadiusRow
            Matrix = mlistack(ll-RadiusRow:ll+RadiusRow,kk-RadiusCol:kk+RadiusCol,:);
            Ref = Matrix(InitRow,InitCol,:);
            T = BWStest(repmat(Ref(:),[1,CalWin(1)*CalWin(2)])...
                ,reshape(Matrix,[CalWin(1)*CalWin(2),npages])',Alpha);   
            temp=reshape(~T,[CalWin(1),CalWin(2)]);
            %connection
            LL = bwlabel(temp);
            SHP.PixelInd(:,num)=LL(:)==LL(InitRow,InitCol);       
            num=num+1;
            if num == all_step * p;
                disp(['progress: ', num2str(5*p),'%']);
                p = p+1;
            end              
        end
    end
end
%SHPs map            
SHP.BroNum = sum(SHP.PixelSub~=0,1);
SHP.BroNum = reshape(SHP.BroNum(:),[nlines,nwidths]);
SHP.BroNum = single((SHP.BroNum-1));          
SHP.CalWin = CalWin;            
%toc;
% figure;imagesc(SHP.BroNum);axis image off
% ti=title ('Homogeneous Pixel Number');
% set(ti,'fontweight','bold');
%disp(['SHP_SelPoint operation completed in ',int2str(t),' second(s).']);
%disp('Done!');            
