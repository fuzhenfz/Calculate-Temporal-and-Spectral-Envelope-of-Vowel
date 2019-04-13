% extract spectral envelop by cepstrum method
clear

Vowel = ['a','e','i','o','u','v'];
Tone = ['1','2','3','4'];
for i1=1:length(Vowel)
    for i2 = 1:length(Tone)
        StiPath = [cd '\record\' Vowel(i1) Tone(i2) '_p.wav'];
        [Sti,fs] = audioread(StiPath);
        if mod(length(Sti),2), Sti(end) = []; end
        PreWeighted = filter([1 -.99],1,Sti);                          % 预加重
        % u=x;
        winlen = length(PreWeighted);                                  % 帧长
        CepstLen = 45;                                        % 倒频率上窗函数的宽度
        winlen2 = winlen/2;
        fre = (0:winlen2-1)*fs/winlen;                        % 计算频域的频率刻度
        PreWeighted2 = PreWeighted.*hamming(winlen);		                     % 信号加窗函数
        PreWeightedFFT = fft(PreWeighted2);                                       % 按式(9-2-1)计算
        PreWeightedFFT_abs = log(abs(PreWeightedFFT(1:winlen2)));                      % 按式(9-2-2)计算
        Cepst = ifft(PreWeightedFFT_abs);                               % 按式(9-2-3)计算
        cepst = zeros(1,winlen2);
        cepst(1:CepstLen) = Cepst(1:CepstLen);                 % 按式(9-2-5)计算
        cepst(end-CepstLen+2:end) = Cepst(end-CepstLen+2:end);
        SpectralEnvelope_0 = real(fft(cepst));                          % 按式(9-2-6)计算
        % downsample
        SpectralEnvelope_temp = resample(SpectralEnvelope_0,ceil(length(SpectralEnvelope_0)/fs*512),length(SpectralEnvelope_0));% downsample to 64 Hz
        SpectralEnvelope((i1-1)*4+i2,:) = SpectralEnvelope_temp;
        figure
        plot(fre,PreWeightedFFT_abs,'k:')
        hold on
        plot(fre(1:length(SpectralEnvelope_0)/length(SpectralEnvelope_temp):length(SpectralEnvelope_0)),SpectralEnvelope_temp,'b','linewidth',2)
        axis([1 8000 -10 4])
        titleName = [Vowel(i1) Tone(i2)];
        title(titleName)
        FileName = [Vowel(i1) Tone(i2) '.jpg'];
        myFileName = [cd '\SpectralEnvelope\figure\' FileName];
        print('-djpeg',myFileName,'-r100');
        close
        
    end
end
% save([cd '\SpectralEnvelope\cepstrum.mat'],'SpectralEnvelope')

%% calculate similarity across all Tone combination
clear
load([cd '\SpectralEnvelope\cepstrum.mat']);% downsample to 64 Hz
SpectralEnvelope = SpectralEnvelope(:,1:22);% <4 kHz
ToneXY = [1 2;1 3;1 4;2 3;2 4;3 4;1 1;2 2;3 3;4 4];

for i0=1:size(ToneXY,1)
    clear dis_mse dis_fd dis_corr dis_Eucli
    ToneX = ToneXY(i0,1);
    ToneY = ToneXY(i0,2);
    i3 = 1;
    if ToneX~=ToneY
        for i1 = ToneX:4:24
            for i2=ToneY:4:24
%                 disp([num2str(i1),'-',num2str(i2)])
                dis_mse(i3) = mse(SpectralEnvelope(i1,:),SpectralEnvelope(i2,:));
                dis_fd(i3) = frechet((1:size(SpectralEnvelope,2))',SpectralEnvelope(i1,:)',(1:size(SpectralEnvelope,2))',SpectralEnvelope(i2,:)');
                dis_corr(i3) = corr(SpectralEnvelope(i1,:)',SpectralEnvelope(i2,:)');
                dis_Eucli(i3) = mean((SpectralEnvelope(i1,:) - SpectralEnvelope(i2,:)));
                i3 = i3+1;
            end
        end
    else % ToneX=ToneY
        for i1 = ToneX:4:24
            for i2=i1+4:4:24
                dis_mse(i3) = mse(SpectralEnvelope(i1,:),SpectralEnvelope(i2,:));
                dis_fd(i3) = frechet((1:size(SpectralEnvelope,2))',SpectralEnvelope(i1,:)',(1:size(SpectralEnvelope,2))',SpectralEnvelope(i2,:)');
                dis_corr(i3) = corr(SpectralEnvelope(i1,:)',SpectralEnvelope(i2,:)');
                dis_Eucli(i3) = mean((SpectralEnvelope(i1,:) - SpectralEnvelope(i2,:)));
                i3 = i3+1;
            end
        end
    end
    dist_mse(i0) = mean(dis_mse);
    dist_fd(i0) = mean(dis_fd);
    dist_corr(i0) = mean(dis_corr);
    dist_Eucli(i0) = mean(dis_Eucli);
end

%% calculate similarity across all Vowel combination
clear
Vowel = ['a','e','i','o','u','v'];
load([cd '\SpectralEnvelope\cepstrum.mat']);% downsample to 64 Hz
SpectralEnvelope = SpectralEnvelope(:,1:22);% <4 kHz
VowelXY = {
    'a','e';'a','i';'a','o';'a','u';'a','v';...
    'e','i';'e','o';'e','u';'e','v';...
    'i','o';'i','u';'i','v';...
    'o','u';'o','v';...
    'u','v';...
    'a','a';'e','e';'i','i';'o','o';'u','u';'v','v'};

for i0=1:size(VowelXY,1)
    clear dis_mse dis_fd dis_corr dis_Eucli
    VowelX = VowelXY{i0,1};
    VowelY = VowelXY{i0,2};
    indexVX = strfind(Vowel,VowelX);
    indexVY = strfind(Vowel,VowelY);
    i3 = 1;
    if indexVX~=indexVY
        for i1 = ((indexVX-1)*4+1):(indexVX*4)
            for i2 = ((indexVY-1)*4+1):(indexVY*4)
%                 disp([num2str(i1),'-',num2str(i2)])
                dis_mse(i3) = mse(SpectralEnvelope(i1,:),SpectralEnvelope(i2,:));
                dis_fd(i3) = frechet((1:size(SpectralEnvelope,2))',SpectralEnvelope(i1,:)',(1:size(SpectralEnvelope,2))',SpectralEnvelope(i2,:)');
                dis_corr(i3) = corr(SpectralEnvelope(i1,:)',SpectralEnvelope(i2,:)');
                dis_Eucli(i3) = mean((SpectralEnvelope(i1,:) - SpectralEnvelope(i2,:)));
                i3 = i3+1;
            end
        end
    else % VowelX=VowelY
        for i1 = ((indexVX-1)*4+1):(indexVX*4)
            for i2=i1+1:(indexVX*4)
%                 disp([num2str(i1),'-',num2str(i2)])
                dis_mse(i3) = mse(SpectralEnvelope(i1,:),SpectralEnvelope(i2,:));
                dis_fd(i3) = frechet((1:size(SpectralEnvelope,2))',SpectralEnvelope(i1,:)',(1:size(SpectralEnvelope,2))',SpectralEnvelope(i2,:)');
                dis_corr(i3) = corr(SpectralEnvelope(i1,:)',SpectralEnvelope(i2,:)');
                dis_Eucli(i3) = mean((SpectralEnvelope(i1,:) - SpectralEnvelope(i2,:)));
                i3 = i3+1;
            end
        end
    end
    dist_mse(i0) = mean(dis_mse);
    dist_fd(i0) = mean(dis_fd);
    dist_corr(i0) = mean(dis_corr);
    dist_Eucli(i0) = mean(dis_Eucli);
end





