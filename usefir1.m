function [h] = usefir1(mode,n,fp,fs,window,sample,r)
% mode:1--HP; 2--LP; 3--BP; 4--BF
% n:order

if window == 1
    w = boxcar(n+1);
end
if window == 2
    w = triang(n+1);
end
if window == 3
    w = bartlett(n+1);
end
if window == 4
    w = hamming(n+1);
end
if window == 5
    w = hanning(n+1);
end
if window == 6
    w = blackman(n+1);
end
if window == 7
    w = kaiser(n+1,r);
end
if window == 8
    w = chebwin(n+1,r);
end
wp = 2*fp/sample;
ws = 2*fs/sample;
if mode == 1
    h = fir1(n,wp,'high',w);
end
if mode == 2
    h = fir1(n,wp,'low',w);
end
if mode == 3
    h = fir1(n,[wp,ws],w);
end
if mode == 4
    h = fir1(n,[wp,ws],'stop',w);
end
