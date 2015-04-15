%*************************************************************
% Implementation of Exposure Fusion
%
% Usage:
%   result = exposure_fusion(I,nlev);
%   Arguments:
%     'I': represents a stack of N color images (at double
%       precision). Dimensions are (height x width x 3 x N).
%     'nlev': decomposition level of BLP pyramid.
%*************************************************************

function R = exposure_fusion(I,nlev)

r = size(I,1);
c = size(I,2);
N = size(I,4);

%% get weight maps
% ExposureWeight: local cue for exposure level measurement
EW=ExposureWeight(I);
% ImportanceWeight: saliency cue using JND
IW=ImportanceWeight(I);
% compute three cues and combines them into a weight map
% InvarianceWeight: global cue for exposure level assessment, esp. moving objects detection for dynamic scene
W=EW.*InvarianceWeight(I,9).*IW;

% EWN: exposure guidance
EWN=EW+1e-12;
EWN=EWN./repmat(sum(EWN,3),[1 1 N]);
% GE:exposure guidance matrix
GE=EWN>0.001;
clear EW EWN

% IWN: importance guidance
IWN=IW+1e-12;
IWN=IWN./repmat(sum(IWN,3),[1 1 N]);
%importance guidance matrix
GJ=IWN.*GE;
clear IW IWN
%normalize weights: make sure that weights sum to one for each pixel
W = W + 1e-12; %avoids division by zero
W = W./repmat(sum(W,3),[1 1 N]);
%% fusion
% create empty pyramid
pyr = gaussian_pyramid(zeros(r,c,3),nlev);

% initialize parameters
sigma_r = 0.4;
%每个像素放大程度不同，用权值控制
alpha = 0.25.^W;
% beta=0,BLP;beta=1,standard Laplacian
beta = 0;
colorRemapping = 'rgb';
domain = 'lin';
% gamma: for base boosting
gamma=[0.8 1 1.2];

% multiresolution blending
fprintf('building BLP...\n');
for i = 1:N
    % construct pyramid from each input image
	pyrW = gaussian_pyramid(W(:,:,i));
    fprintf('image %d-%d:\n',N,i);
	pyrI = lapfilter(I(:,:,:,i),sigma_r,alpha(:,:,i),beta,colorRemapping,domain,GE(:,:,i),GJ(:,:,i),gamma(i),nlev);
    % blend
    for l = 1:nlev
        w = repmat(pyrW{l},[1 1 3]);
        pyr{l} = pyr{l} + w.*pyrI{l};
    end
end
%% reconstruct
R = reconstruct_laplacian_pyramid(pyr);
R=postprocessing(R,domain,beta,colorRemapping);
