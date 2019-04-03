function y = linsurr_ind(x)
%Creates surrogate multichannel data, by adding random numbers to phase component
%of all channel data.

% Email from Mike Breakspear
%Using independent phase randomization tests the null that the apparent
%effects are due to finite sample length and/or auto-correlations within
%each channel; (this code)
%Using uniform phase randomization tests the above null and/or the presence
%of cross-correlations between the time series; (lin_surr)
%Using uniform phase randomization tests the above and/or the presence of a
%static nonlinear transformation of the time series amplitude. (lin_surr)

[r,c] = size(x);
% if r < c % modified by EJC, don't transpose if c > r just give it time-by-channel
% 	x = x.';   % make each column a timeseries
% end;
[n,cc] = size(x);
m = 2^nextpow2(n);
y = fft(real(x),m);
phsrnd=zeros(m,cc);

for i=1:cc 
    phsrnd(2:m/2,i)=rand(m/2-1,1)*2*pi; 
    phsrnd(m/2+2:m,i)=-phsrnd(m/2:-1:2,i);
end
y = y.*exp(phsrnd*sqrt(-1));
y = ifft(y,m);
y = y(1:n,:);
% if r < c % modified by EJC, don't transpose if c > r just give it time-by-channel
%    y = y.';
% end
y=real(y);    %small imag. component created by rounding error
