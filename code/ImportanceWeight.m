%**********************************
% saliency measurement JND
% parameters:
%           img:input sequence
%**********************************
% Matlab codes for pixel-based JND (Just-Noticeable Difference) model
% This program implements the pixel-based JND (Just-Noticeable Difference) model described in the following papers: 
% X.K. Yang, W.S. Lin, Z.K. Lu, E. Ong and S. Yao, Motion-compensated Residue Pre-processing 
% in Video Coding Based on Just-noticeable-distortion Profile, IEEE Trans. Circuits and Systems for Video Technology, vol.15(6),
% pp.742-750, June, 2005.
function W=ImportanceWeight(img)
I=uint8(img*255);
N = size(I,4);
W = zeros(size(I,1),size(I,2),N);
for i=1:N
    gray=rgb2gray(I(:,:,:,i));
    W(:,:,i)=getJND(gray)*10/255;
end
end

function JND=getJND(I,type)
% JND=JND_img(I,type)

%Initialize
if nargin<1
    close all;clc;
    I =imread('demo.bmp');
    I =rgb2gray(I);
%     Test=1;
end

if nargin<2
    type='Yang';
end

I=double(I);
[H,W]=size(I);
JNDl=zeros(H,W); % JND: Luminance component
% JNDt=zeros(H,W); % JND: Texture component
T=zeros(H,W,2);  % Temp matrix
bg = func_bg(I); % Average background luminance
Gm = func_Gm(I,type); % Maximal weighted average of gradients (with and without identify edge)
%%
% Calculating JNDl
T0 = 17;
GAMMA = 3/128;
for i = 1:H
   for j = 1:W
      if bg(i,j) <= 127
         JNDl(i,j) = T0*(1-sqrt(bg(i,j)/127))+3;
      else
         JNDl(i,j) = GAMMA*(bg(i,j)-127)+3;
      end
   end
end
%%
% Calculating JNDt
LANDA = 1/2;
alpha = 0.0001*bg+0.115;
belta = LANDA-0.01*bg;
JNDt = Gm.*alpha + belta;
%%
% Calculating final JND
T(:,:,1)=JNDl;
T(:,:,2)=JNDt;
C_TG = 0.3;
JND = sum(T,3)-C_TG*min(T,[],3);
%%
end



%%
function output = func_bg(input)
Mask=[1 1 1 1 1
      1 2 2 2 1
      1 2 0 2 1
      1 2 2 2 1
      1 1 1 1 1];
output = filter2(Mask,input)/32;
end
%%
function output = func_Gm(input,type)
G1=[0  0  0  0  0
    1  3  8  3  1
    0  0  0  0  0
   -1 -3 -8 -3 -1
    0  0  0  0  0];

G2=[0 0  1  0  0
    0 8  3  0  0
    1 3  0 -3 -1
    0 0 -3 -8  0
    0 0 -1  0  0];

G3=[0  0  1 0 0
    0  0  3 8 0
   -1 -3  0 3 1
    0 -8 -3 0 0
    0  0 -1 0 0];

G4=[0 1 0 -1 0
    0 3 0 -3 0
    0 8 0 -8 0
    0 3 0 -3 0
    0 1 0 -1 0];
[H,W]=size(input);
grad=zeros(H,W,4);
grad(:,:,1) = filter2(G1,input)/16;
grad(:,:,2) = filter2(G2,input)/16;
grad(:,:,3) = filter2(G3,input)/16;
grad(:,:,4) = filter2(G4,input)/16;
Gm = max(abs(grad),[],3);

edge_threshold = 0.5;
img_edge = edge(input,'canny',edge_threshold);
se = strel('disk',6,6);
img_edge = imdilate(img_edge,se);
img_supedge = 1-0.95*double(img_edge);
gaussian_kernal = fspecial('gaussian',7,0.8);
img_supedge = filter2(gaussian_kernal,img_supedge);

switch type 
case 'Chou'
    output = Gm;
case 'Yang'
    output = Gm.*img_supedge;
end

end