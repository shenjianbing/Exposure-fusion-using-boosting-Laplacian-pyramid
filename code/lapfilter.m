% Laplacian Filtering 
%   - public Matlab implementation for reproducibility
%   - about 30x slower than our single-thread C++ version
%
% This script implements edge-aware detail and tone manipulation as 
% described in Paris, Hasinoff, and Kautz, "Laplacian Filters: Edge-aware 
% Image Processing with a Laplacian Pyramid", ACM Transactions on 
% Graphics (Proc. SIGGRAPH 2011), 30(4), 2011.
%
% This is a wrapper around the core algorithm (see lapfilter_core.m).
% It defines the remapping function, the color treatment, the processing 
% domain (linear or log), and it implements simple postprocessing
% for our tone mapping results.
%
% Arguments:
%   image 'I'
%   pixel-wise remapping parameters 'sigma_r', 'alpha', 'beta'
%   remapping method for color images 'colorRemapping' ['rgb' or 'lum']
%   processing domain 'domain' ['lin' or 'log']
%
% sam.hasinoff@gmail.com, April 2011
%

function [R,Iratio] = lapfilter(I,sigma_r,alpha,beta,colorRemapping,domain,GE,GJ,gamma,nlev)

% interpret the input arguments
if size(I,3)==1
    I = repmat(I,[1 1 3]); 
end
if strcmp(domain,'log')
    sigma_r = log(sigma_r);
end

% detail remapping function
noise_level = 0.01;

function out = fd(I,g0,sigma_r,alpha_sub,gj)
    % 只对rgb向量大小增强
    dnrm=sqrt(sum((I-g0).^2,3));
    % gj: importance guidance
    n=0.2.^gj;%n在[0.7,1],n越小out越大，对比度越高
    theta=0.5;%theta越小颜色越艳
    out=dnrm.^n./(dnrm.^n+theta.^n);
    % 保持向量方向
    unit = (I-g0)./repmat(eps + dnrm,[1 1 3]);
    idx=find(alpha_sub<1);
    tau = smooth_step(noise_level,2*noise_level,dnrm(idx)*sigma_r);%平坦区域（1-tau）不增强
    out(idx) = tau.*out(idx) + (1-tau).*dnrm(idx);
    out=unit.*repmat(out,[1 1 3]); 
end

% define the overall pixel-wise remapping function r, using
% the threshold sigma_r for edge-detail separation
switch colorRemapping
    case 'rgb'
        % process pixels as vectors in RGB color space
        r = @(i,g0,alpha_sub,GJ)(r_color(i,g0,alpha_sub,sigma_r,GJ,@fd));%匿名函数r,调用时用Iremap=r(参数列表);参数传递给(i,g0)
    case 'lum'
        % save RGB color ratios for later, process the luminance
        IY = luminance(I);
        Iratio = I ./ repmat(IY+eps,[1 1 3]);
        I = IY;
        r = @(i,g0,alpha_sub)(r_gray(i,g0,alpha_sub,sigma_r,@fd,@fe));
    otherwise
        error('invalid color remapping');
end

% define the processing domain
switch domain
    case 'lin',
        to_domain   = @(I) I;

    case 'log', 
        to_domain   = @(I) log(I + eps);

    otherwise
        error('invalid domain');
end

% % call the core Laplacian filtering algorithm
if  beta==1  
    R = laplacian_pyramid(I);
else
    I = to_domain(I);
    R = lapfilter_core(I,r,GE,GJ,alpha,gamma,nlev);   %return new Laplacian pyr
%     R = from_domain(R);
end
%% helper functions

% smooth step edge between (xmin,0) and (xmax,1)
function y = smooth_step(xmin,xmax,x)
    y = (x - xmin)/(xmax - xmin);
    y = max(0,min(1,y));
    y = y.^2.*(y-2).^2;
end

% convert RGB to grayscale intensity
function Y = luminance(I)
    switch size(I,3),
        case 1, Y = I;
        case 3, Y = (20*I(:,:,1) + 40*I(:,:,2) + I(:,:,3))/61;
    end
end

% color remapping function
function inew = r_color(patchI,g0,alpha_sub,sigma_r,gj,fd)
    %复制高斯滤波系数g0成patch
    g0 = repmat(g0,[size(patchI,1) size(patchI,2) 1]);
    %调用fd增强高频成分
    inew=g0+fd(patchI,g0,sigma_r,alpha_sub,gj);
end

end
