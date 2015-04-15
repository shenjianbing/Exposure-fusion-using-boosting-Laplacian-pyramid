%lowpass,weightmap额外再下采样两级
function I=extDownsample(lowpass,filter)
I=lowpass;
for l = 1:2
    I = downsample(I,filter);
end

end