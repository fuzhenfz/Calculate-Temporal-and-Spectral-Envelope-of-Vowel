%% calculate temporal envelope, Hilbert + 1 bands + 8 Hz LP
clear

LP = 8; % Hz
window = 4;% hamming
order = 100;

StiPath = ['a1.wav'];
[Sti,fs] = audioread(StiPath);
Envelope_temp = abs(hilbert(Sti))';
hn = usefir1(2,order,LP,[],window,fs);%2 = LP
EnvelopeFilt = conv(hn,Envelope_temp);
Envelope_temp = resample(EnvelopeFilt,ceil(length(EnvelopeFilt)/fs*64),length(EnvelopeFilt));% downsample to 64 Hz
Envelope = Envelope_temp;

figure
plot(1:length(Sti),Sti)
hold on
plot(1:length(Sti)/size(Envelope_temp,2):length(Sti),Envelope_temp)
