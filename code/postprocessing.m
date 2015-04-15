% postprocessing
function R=postprocessing(R,domain,beta,colorRemapping)

switch domain
    case 'lin',
        from_domain = @(R) R;
    case 'log', 
        from_domain = @(R) exp(R) - eps;
    otherwise
        error('invalid domain');
end
R = from_domain(R);
if strcmp(domain,'log') && beta<=1
    % for tone mapping, remap middle 99% of intensities to
    % fixed dynamic range using a gamma curve
    DR_desired = 100;
    prc_clip = 0.5;
    RY = luminance(R);
    Rmax_clip = prctile(RY(:),100-prc_clip);
    Rmin_clip = prctile(RY(:),prc_clip);
    DR_clip = Rmax_clip/Rmin_clip;
    exponent = log(DR_desired)/log(DR_clip);
    R = max(0,R/Rmax_clip) .^ exponent;
end
if strcmp(colorRemapping,'lum')
    % if working with luminance, reintroduce color ratios 
    R = repmat(R,[1 1 3]) .* Iratio;
end
% clip out of bounds intensities
R = max(0,R);
if beta<=1
    R = min(1,R);  
end
if strcmp(domain,'log') && beta<=1
    % for tone mapping, gamma correct linear intensities for display
    gamma_val = 2.2;
    R = R.^(1/gamma_val);
end