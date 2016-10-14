function [border,index]=lumen_segmentation_DP(im,cathe_pos,connectivity,band_width)
if nargin<4
    band_width=min(30,cathe_pos-1);
end
if nargin<3
    connectivity=11;
end
if nargin<2
    cathe_pos=80;
end

add_lines=20; % To guarantee the contour is closed, there may be better ways to do it
im=[im(end-add_lines+1:end,:);im;im(1:add_lines,:)];

M=size(im,1);N=size(im,2);

f=zeros(M,N);
for j=cathe_pos:N-band_width
    f(:,j)=mean(im(:,j:min(j+band_width,N)),2)-mean(im(:,j-band_width:j),2);
end
% for j=cathe_pos:N-band_width
%     f(:,j)=2*mean(im(:,j:min(j+band_width,N)),2)-mean(im(:,j-band_width:j),2)-mean(im(:,min(j+10,N):min(j+band_width+60,N)),2);
% end

%[border,index]=Dynamic_Programming(f,connectivity);
[border,index]=DP_fast(f',connectivity);
border=border';

border=border(add_lines+1:end-add_lines,:,:);
index=index(add_lines+1:end-add_lines);

border=connect_border(border);