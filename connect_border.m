function [cont_border,index]=connect_border(dis_border)
% dis_border should be binary image
% index stores all the column position of the stent strut at each row
border=dis_border;
M=size(border,1);N=size(border,2);
index=zeros(M,1);
for m = 1:M
    for n = 1:N
        if border(m,n) == 1
            index(m)= n;
            break;
        end
    end
end
% avoid jump
for m=1:M
    if index(m)~=0
        if m>2&&m<M
            if index(m-1)~=0 && index (m+1)~=0
                if abs(index(m)-index(m-1))>10&&abs(index(m)-index(m+1))>10
                    index(m)=round((index(m+1)+index(m-1))/2);
                    border(m,:)=border(m,:)*0;
                    border(m,index(m))=1;
                end
            elseif index(m-1)~=0 && index (m+1)==0
                if abs(index(m)-index(m-1))>10
                    index(m)=index(m-1)+sign(index(m)-index(m-1));
                    border(m,:)=border(m,:)*0;
                    border(m,index(m))=1;
                end
            elseif index(m-1)==0 && index (m+1)~=0
                if abs(index(m)-index(m+1))>10
                    index(m)=index(m+1)-sign(index(m+1)-index(m));
                    border(m,:)=border(m,:)*0;
                    border(m,index(m))=1;
                end
            end
        end
    end
end

for m = 1:M
    if index(m)==0
        p=m;q=m;
        wrap_p = 0;wrap_q = 0;
        % find the next point, save the row position into p
        while 1
            p=p+1;
            if p == M+1
                p = 1;
                wrap_p=M;
            end
            if index (p)~=0
                break;
            end
        end
        % find the last point, save the row position into p
        while 1
            q=q-1;
            if q == 0
                q = M;
                wrap_q=M;
            end
            if index (q)~=0
                break;
            end
        end
        index(m)=round(linear_interpolation(p+wrap_p,q-wrap_q,index(p),index(q),m));
    end
end

for m =1:M
    if m == M
        border(m,index(m))=1;
        break;
    end
    if index(m+1)== index(m) ||abs(index(m+1)-index(m)) == 1
        border(m,index(m))=1;
        continue;
    elseif index(m+1)>index(m)||index(m+1)==index(m)
        border(m,index(m):index(m+1)-1)=1;
    else
        border(m,index(m+1):index(m)-1)=1;
    end
end
cont_border=border;
%--------------------------------------------------------------------------