function ac = autocorrelationForCircularVars(S,params)
% NB: this fxn assumes that Ndims = 1

% params
[Ntraj,T] = size(S);
N = params.N;
NSmin = params.smin(:,strcmp(params.NS,params.mods));
NSmax = params.smax(:,strcmp(params.NS,params.mods));
trajmin = NSmin;
trajmax = N/(N-1)*(NSmax - NSmin) + NSmin;
tauMax = T-1;





% autocorrelation of the data
if 0
    acx = 0;
    acy = 0;
    for iTraj = 1:Ntraj
        s = scalefxn(squeeze(S(iTraj,:)),trajmin,trajmax,0,2*pi);
        acx = acx + xcorr(cos(s),cos(s),tauMax,'none');
        acy = acy + xcorr(sin(s),sin(s),tauMax,'none');
        %%% using the mean of these is identical to using the real part of:
        %
        %     % autocorrelation
        %     Nfft = 2^nextpow2(2*T-1);
        %     z = exp(1i*s);
        %     r = ifft( fft(z,Nfft) .* conj(fft(z,Nfft)) );
        %
        %     % rearrange and keep values corresponding to lags: -(len-1):+(len-1)
        %     r = [r(end-T+2:end) r(1:T)];
        %
        %%% It's also quite close to using the magnitude of this, especially in
        %%% the region of interest
        
    end
    ac = [acx; acy]/Ntraj;
    
else
    ac = zeros(1,tauMax);
    for iTraj = 1:Ntraj
        s = scalefxn(squeeze(S(iTraj,:)),trajmin,trajmax,0,2*pi);
        for tau = 0:(tauMax-1)
            a = s(1:(end-tau))';
            b = s((1+tau):end)';
            n = length(a);
            
            num3 = (cos(a)'*cos(b))*(sin(a)'*sin(b)) - (cos(a)'*sin(b))*(sin(a)'*cos(b));
            den3 = sqrt((n^2 - sum(cos(2*a))^2 - sum(sin(2*a))^2)*...
                (n^2 - sum(cos(2*b))^2 - sum(sin(2*b))^2))/4;
            ac(tau+1) = ac(tau+1) + num3/den3;
        end
    end
    ac = ac/Ntraj;
    ac = [fliplr(ac) ac(2:end)]; %%% symmetrize
end

end













