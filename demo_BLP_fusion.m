function demo_BLP_fusion

clear all
close all
addpath('code');
I= load_images('images\land');
%%
% recommend small size for testing
I = min(1,max(0, imresize(I,1/16) )); 
fprintf('processing %d images sequence:land\n',size(I,4));
% nlev: pyramid level
nlev=2;
R = exposure_fusion(I,nlev);
% figure('Name','Result'); 
% imshow(R); 
imwrite(R,'result\land_01.bmp');

