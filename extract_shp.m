% extract SHP groups phase in SLCs---------------------------------------------------
function [SHP_check,SHPpatch,id_center,ix_l,ix_p]=extract_shp(SHP,idd,hW_l,hW_w,p_N,p_M,minshp)

j=ceil(idd/p_M);
i=mod(idd,p_M);
if i==0
    i=p_M;
end
SHPpatch=logical(reshape(SHP.PixelSub(:,idd),2*hW_l+1,2*hW_w+1));
%SHPpatch=squeeze(SHP(i,j,:,:));
SHPpatch=bwlabel(logical(SHPpatch),4); % Connectivity check: 8 connectivity, or 4
SHPpatch(SHPpatch==SHPpatch(hW_l+1,hW_w+1))=1000;
SHPpatch(SHPpatch<1000)=0;
SHPpatch(SHPpatch==1000)=1;

% if sum(SHPpatch(:)) < minshp    % set the minimum SHP number: 25
%     %SHPpatch(:,:)=0;
%     SHPpatch(hW_l+1-2:hW_l+1+2,hW_w+1-2:hW_w+1+2)=1;
% end
SHP_check= sum(sum(SHPpatch(:,:,1)));

ix_l=[i-hW_l:i+hW_l];        ix_p=[j-hW_w:j+hW_w];
ix_l=MOD(ix_l-1,p_M);        ix_p=MOD(ix_p-1,p_N);
SHPpatch(SHPpatch==0)=nan;
%
SHPpatch0=SHPpatch;
SHPpatch0(hW_l+1,hW_w+1)=1000;
SHPpatch0=reshape(SHPpatch0,(2*hW_l+1)*(2*hW_w+1),1);
SHPpatch0(isnan(SHPpatch0))=[];
id_center=find(SHPpatch0==1000);
%% ------------------------------

end
function [c]=MOD(a,b)
c=rem(b+rem(a,b),b)+1;
end