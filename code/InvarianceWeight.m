%**********************************
% global exposure level measurement
% parameters:
%           img:input sequence
%           l:block size
%**********************************
function Su=InvarianceWeight(img,l)%l=9
%% 初始化
[r1,c1,h,n]=size(img);
gauss = fspecial('gaussian', [1 9], 1.2);
M=zeros(r1,c1,n);
GX=zeros(r1,c1,n);
GY=zeros(r1,c1,n);
G=zeros(r1,c1,n);
ANG=zeros(r1,c1,n);
T=zeros(r1,c1,h,n);
I=zeros(r1,c1,n);
%% 只用v通道
for j=1:n
    T(:,:,:,j)=rgb2hsv(img(:,:,:,j));
    I(:,:,j)=T(:,:,3,j);
end
%% 求梯度方向
for i=1:n  
    extend = pad_reflect(abs(I(:,:,i)), 4.5);
    Te=conv2(gauss,gauss',extend, 'valid');%2d gaussian
    G(:,:,i) = Te;
    [GX(:,:,i),GY(:,:,i)]=gradient(G(:,:,i));%一阶导数
    M(:,:,i)=sqrt(GX(:,:,i).*GX(:,:,i)+GY(:,:,i).*GY(:,:,i));%大小
    GX(:,:,i)=GX(:,:,i)+1e-25;
    %ANG(:,:,i)为梯度方向
    tmp=atan2((GY(:,:,i)),(GX(:,:,i)));    
    tmp(tmp(:,:)<0) = 2*pi + tmp(tmp(:,:)<0);  %[0, 2pi];
    ANG(:,:,i)=tmp;
    
end
clear  Te extend G GX GY tmp T I gauss
factor=1;%调整mask区域大小
%% 利用均值滤波器来算窗口的均值,很快啊
Avg=cell(n,n);
for j=1:n
    for i=1:n
        Avg{j,i}=zeros(r1,c1);
    end
end

avg=fspecial('average',[2*l+1 2*l+1]);
for i=1:n-1
    for j=i+1:n
        tt=min(2*pi-abs(ANG(:,:,i)-ANG(:,:,j)),abs(ANG(:,:,i)-ANG(:,:,j)));
        tt=sin(tt);
        Sine=M(:,:,i).*tt;
        Avg{j,i}=imfilter(Sine,avg,'replicate');
        tem=1-Avg{j,i}*factor;
        tem(tem<0)=0;
        Avg{j,i}=tem;
        
        Sine=M(:,:,j).*tt;
        Avg{i,j}=imfilter(Sine,avg,'replicate');
        tem=1-Avg{i,j}*factor;
        tem(tem<0)=0;
        Avg{i,j}=tem;
    end
end
%% 累积
clear avg tt Sine M
Su=zeros(r1,c1,n);
for i=1:n
    for j=1:n
        if i~=j%忽略和自己的差，Su(:,:,i)范围在[0,n-1]
           tem=Su(:,:,i)+Avg{j,i};
           tem(tem>0.01)=1;
%            tem=bwmorph(tem,'erode');%some situation need this step to eliminate noise
%            tem=bwmorph(tem,'dilate');
           Su(:,:,i)=tem;
        end
    end
end
% 避免任意两幅图都不一样的情况
tem=sum(Su,3);
t=Su(:,:,1);
t(tem==0)=1;
Su(:,:,1)=t;
end