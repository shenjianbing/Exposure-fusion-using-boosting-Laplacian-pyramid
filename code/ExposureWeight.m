%**********************************
% local exposure level measurement
% parameters:
%           I:input sequence
%**********************************
function W=ExposureWeight(I)

a=3.2;%[-3.71,-2.69]
b=-1.3;%[-1.71,-0.88]
c=1/0.4;%[0.30,0.49]
N = size(I,4);

W = zeros(size(I,1),size(I,2),N);
for i=1:N
    I(:,:,:,i)=(log(I(:,:,:,i).^c./(1-I(:,:,:,i).^c))-b)/a;
    W(:,:,i)=rgb2gray(1-(abs(I(:,:,:,i))));
end

end