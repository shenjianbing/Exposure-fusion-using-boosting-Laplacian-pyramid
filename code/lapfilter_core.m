% Laplacian Filtering 
%   - public Matlab implementation for reproducibility
%   - about 30x slower than our single-thread C++ version
%
% This script implements the core image processing algorithm 
% described in Paris, Hasinoff, and Kautz, "Laplacian Filters: Edge-aware 
% Image Processing with a Laplacian Pyramid", ACM Transactions on 
% Graphics (Proc. SIGGRAPH 2011), 30(4), 2011.
%
% Processes an input image using a general pointwise remapping function 
% r(I,g0) in a pyramid-based way. Its running time is O(N log N), where 
% N is the number of pixels.
%
% Most of the code is bookkeeping to isolate the subpyramid contributing 
% to a particular Laplacian coefficient. See below for a 14-line naive 
% O(N^2) implementation that gives identical results.
%
% Arguments:
%   image 'I'
%   pixel-wise remapping function 'r', expects arguments r(I,g0)
%
% sam.hasinoff@gmail.com, March 2011
%

function L = lapfilter_core(I,r,GE,GJ,alpha,gamma0,nlev)%return new Laplacian pyr L,cell type

G = gaussian_pyramid(I,nlev); 
GJ1=GJ;
% build up the result, one Laplacian coefficient at a time
L = laplacian_pyramid(zeros(size(I)),nlev);  % allocate space for result
tic;
for lev0 = 1:length(L)-1
    hw = 3*2^lev0 - 2;  % half-width of full-res footprint (conservative)
    fprintf('level %d (%dx%d), footprint %dx%d ...   0%',lev0,size(G{lev0},1),size(G{lev0},2),min(2*hw+1,size(I,1)),min(2*hw+1,size(I,2)));
    for y0 = 1:size(G{lev0},1)
        for x0 = 1:size(G{lev0},2)
            % exposure guidance 
            if GE(y0,x0)==0
                continue;
            end
            % coords in full-res image corresponding to (lev0,y0,x0)
            yf = (y0-1)*2^(lev0-1) + 1;
            xf = (x0-1)*2^(lev0-1) + 1;
            % subwindow in full-res image needed to evaluate (lev0,y0,x0) in result
            % patch大小根据需要分解的子金字塔层数决定
            yrng = [max(1,yf-hw) min(size(I,1),yf+hw)];
            xrng = [max(1,xf-hw) min(size(I,2),xf+hw)];
            Isub = I(yrng(1):yrng(2),xrng(1):xrng(2),:);
            alpha_sub=alpha(yrng(1):yrng(2),xrng(1):xrng(2),:);
            GJ_sub=GJ(yrng(1):yrng(2),xrng(1):xrng(2),:);
            % use the corresponding Gaussian pyramid coefficient to remap
            % the full-res subwindow
            g0 = G{lev0}(y0,x0,:);
            Iremap = r(Isub,g0,alpha_sub,GJ_sub);
            % compute Laplacian pyramid for remapped subwindow
            Lremap = laplacian_pyramid(Iremap,lev0+1,[yrng xrng]);
            
            % bookkeeping to compute index of (lev0,y0,x0) within the
            % subwindow, at full-res and at current pyramid level
            yfc = yf - yrng(1) + 1;
            xfc = xf - xrng(1) + 1;
            yfclev0 = floor((yfc-1)/2^(lev0-1)) + 1;
            xfclev0 = floor((xfc-1)/2^(lev0-1)) + 1;
            
            % set coefficient in result based on the corresponding
            % coefficient in the remapped pyramid
            L{lev0}(y0,x0,:) = Lremap{lev0}(yfclev0,xfclev0,:);            
        end
        fprintf('\b\b\b\b%3d%%',floor(y0/size(G{lev0},1)*100));
    end
    fprintf('\n');
    GE= imresize(GE,1/2);
    GJ1=imresize(GJ1,0.5);
end
%% base boosting   
gj=log(GJ1+1e-10);
nbr=sum(GE(:));
sgj=sum(gj,3)/nbr;
nbr=size(gj,1)*size(gj,2);
t=gamma0*exp(sgj);
gamma=sum(t(:))/nbr;
L{end} = G{end}.^gamma;  % residual not affected

toc;

end

