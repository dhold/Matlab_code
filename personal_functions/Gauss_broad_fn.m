function [broaden_sig,renorm_fct] = Gauss_broad_fn(sigma,omega,signal,dim_to_conv)

static_broad_fn = exp(-(omega-median(omega)).^2/2/sigma^2)/sqrt(2*pi)/sigma;
static_ft = fft(static_broad_fn);
sz = size(signal);
if nargin == 3
    dim_to_conv  = find(sz==length(omega),1,'first');
end
sz(dim_to_conv) = 1; lg = ones(size(sz)); lg(dim_to_conv) = length(omega);
static_ft = reshape(static_ft,lg);

broaden_sig = fftshift(ifft(fft(signal,[],dim_to_conv)...
                .*repmat(static_ft,sz),[],dim_to_conv),dim_to_conv);
%next factor normalises the signal to have the same amplitude as before
renorm_fct = trapz(omega,broaden_sig,dim_to_conv)./trapz(omega,signal,dim_to_conv);