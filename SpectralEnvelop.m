% Calculate spectral envelop by cepstrum method
clear

StiPath = ['a1.wav'];
[Sti,fs] = audioread(StiPath);
if mod(length(Sti),2), Sti(end) = []; end
PreWeighted = filter([1 -.99],1,Sti);
winlen = length(PreWeighted);
CepstLen = 45;% window width
winlen2 = winlen/2;
fre = (0:winlen2-1)*fs/winlen;
PreWeighted2 = PreWeighted.*hamming(winlen);
PreWeightedFFT = fft(PreWeighted2);
PreWeightedFFT_abs = log(abs(PreWeightedFFT(1:winlen2)));
Cepst = ifft(PreWeightedFFT_abs);
cepst = zeros(1,winlen2);
cepst(1:CepstLen) = Cepst(1:CepstLen);
cepst(end-CepstLen+2:end) = Cepst(end-CepstLen+2:end);
SpectralEnvelope_0 = real(fft(cepst));        
SpectralEnvelope_temp = resample(SpectralEnvelope_0,ceil(length(SpectralEnvelope_0)/fs*512),length(SpectralEnvelope_0));% downsample to 64 Hz
SpectralEnvelope((i1-1)*4+i2,:) = SpectralEnvelope_temp;
        
figure
plot(fre,PreWeightedFFT_abs,'k:')
hold on
plot(fre(1:length(SpectralEnvelope_0)/length(SpectralEnvelope_temp):length(SpectralEnvelope_0)),SpectralEnvelope_temp,'b','linewidth',2)
axis([1 8000 -10 4])        
        
