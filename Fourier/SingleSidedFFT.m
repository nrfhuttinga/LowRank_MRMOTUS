function [P1,f]=SingleSidedFFT(X,Fs,d)
% Perform single sided FFT on 'X', which has sample frequency 'Fs' in Hz
%
% Copyright UMC Utrecht, 2020. Written by Niek Huttinga, 2020. For academic purpose only.

sx = size(X);
if numel(sx)>2
    error('Wrong input dimension');
end

if nargin < 2 || isempty(d)
    d = 1;
end

% other dimension than d
d_2 = find(size(SliceData(X,d,1))~=1);

L = sx(d);
for i=1:sx(d_2)
    Y = fft(SliceData(X,d_2,i),[],d);
    P2 = abs(Y/L);
    P1(i,:) = P2(1:round(L/2)+1);
    P1(i,2:end-1) = 2*P1(i,2:end-1);
    f = Fs*(0:round(L/2))/L;
end

if d==1
    P1 = permute(P1,[2 1]);
end