function fibrous_cap=TCFA_boundary_extra_signelRegion(Img_gw_remove,cordX_ori,lumen_index,band_width,inten_width)
M=size(Img_gw_remove,1);N=size(Img_gw_remove,2);
add_lines=12;

r_add_u=max(cordX_ori(1)-add_lines,1);r_add_d=max(cordX_ori(1)-1,1);
add_ind=cordX_ori(1)-r_add_u;
cordX=[r_add_u:r_add_d cordX_ori min(M,cordX_ori(end)+1):min(M,cordX_ori(end)+add_lines)];
%inten_dif=zeros(length(cordX),band_width);
Img_gw_ROI=zeros(length(cordX),band_width+inten_width);
for i=1:length(cordX)
    Img_gw_ROI(i,:)=Img_gw_remove(cordX(i),lumen_index(cordX(i)):lumen_index...
        (cordX(i))+band_width+inten_width-1);
end

inten_dif=TCFA_energy(Img_gw_ROI',7,15,inten_width,7); % The partial cost function
% for r=1:size(Img_gw_ROI,1)
%     inten_dif=diff(Img_gw_ROI(1:end));
% end

connectivity=7;
[border,countour_index]=DP_fast(inten_dif,connectivity);
TCFA_index=double(countour_index)+double(lumen_index(cordX));

GaussianDieOff = .0001;
pw = 1:30;
sig=3;
ssq = sig*sig;
width = max(find(exp(-(pw.*pw)/(2*ssq))>GaussianDieOff));
if isempty(width)
    width = 1;
end
t = (-width:width);
gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);
gau=gau/sum(gau);
W=width;
L=length(TCFA_index);
if L>W
    x1=TCFA_index;
    ad=round(length(gau)/2)-1;
    xx=conv(x1,gau);
    TCFA_index_smooth=xx(ad+1:end-ad);
    TCFA_index_smooth=round(TCFA_index_smooth);
else
    TCFA_index_smooth=TCFA_index;
end

fibrous_cap=zeros(M,N);
for i=1:size(cordX_ori,2)
    fibrous_cap(cordX_ori(i),lumen_index(cordX_ori(i)):TCFA_index_smooth(i+add_ind))=1;
end