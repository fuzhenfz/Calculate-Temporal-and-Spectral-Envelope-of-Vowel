%% extract amplitude envelope
% Hilbert + 8 bands
clear

Vowel = ['a','e','i','o','u','v'];
Tone = ['1','2','3','4'];
cutoffFre = [80 205 405 724 1236 2055 3366 5463 8820];% 8 bands
for i1=1:length(Vowel)
    for i2 = 1:length(Tone)
        StiPath = [cd '\record\' Vowel(i1) Tone(i2) '_p.wav'];
        [Sti,fs] = audioread(StiPath);
        for i3=1:length(cutoffFre)-1
            b = fir1(256,[cutoffFre(i3)/(fs/2),cutoffFre(i3+1)/(fs/2)], hamming(256+1));
            Signal_SubBand = filter(b,1,Sti);
            Envelope_SubBand_temp = abs(hilbert(Signal_SubBand))';
            Envelope_SubBand(i3,:) = resample(Envelope_SubBand_temp,42,length(Envelope_SubBand_temp));% downsample to 80 samples
        end
        Envelope((i1-1)*4+i2,:) = sum(Envelope_SubBand);
        
        %         figure
        %         plot(1:length(Sti),Sti)
        %         hold on,
        %         plot(1:length(Sti)/size(Envelope,2):length(Sti),Envelope((i1-1)*4+i2,:));
        %         titleName = [Vowel(i1) Tone(i2)];
        %         title(titleName)
    end
end
save([cd '\Envelope\Envelope_8band_Hilbert.mat'],'Envelope')

%% Hilbert + 1 bands
clear
Vowel = ['a','e','i','o','u','v'];
Tone = ['1','2','3','4'];
for i1=1:length(Vowel)
    for i2 = 1:length(Tone)
        StiPath = [cd '\record\' Vowel(i1) Tone(i2) '_p.wav'];
        [Sti,fs] = audioread(StiPath);
        Envelope_temp = abs(hilbert(Sti))';
        Envelope_temp = resample(Envelope_temp,42,length(Envelope_temp));% downsample to 80 samples
        output1 = smoothts(Envelope_temp,'b',2);
        Envelope((i1-1)*4+i2,:) = Envelope_temp;
        figure
        plot(1:length(Sti),Sti)
        hold on,
        plot(1:length(Sti)/size(Envelope,2):length(Sti),Envelope((i1-1)*4+i2,:));
        plot(1:length(Sti)/size(Envelope,2):length(Sti),output1)
        titleName = [Vowel(i1) Tone(i2)];
        title(titleName)
    end
end
save([cd '\Envelope\Envelope_1band_Hilbert.mat'],'Envelope')

%% Hilbert + 1 bands + 8 Hz LP
clear
Vowel = ['a','e','i','o','u','v'];
Tone = ['1','2','3','4'];
LP = 8; % Hz
window = 4;% hamming
order = 100;
for i1=1:length(Vowel)
    for i2 = 1:length(Tone)
        StiPath = [cd '\record\' Vowel(i1) Tone(i2) '_p.wav'];
        [Sti,fs] = audioread(StiPath);
        Envelope_temp = abs(hilbert(Sti))';
        hn = usefir1(2,order,LP,[],window,fs);%2 = LP
        EnvelopeFilt = conv(hn,Envelope_temp);
        Envelope_temp = resample(EnvelopeFilt,ceil(length(EnvelopeFilt)/fs*64),length(EnvelopeFilt));% downsample to 64 Hz
        Envelope((i1-1)*4+i2,:) = Envelope_temp;
        figure
        plot(1:length(Sti),Sti)
        hold on
        plot(1:length(Sti)/size(Envelope_temp,2):length(Sti),Envelope_temp)
        titleName = [Vowel(i1) Tone(i2)];
        title(titleName)
        FileName = [Vowel(i1) Tone(i2) '.jpg'];
        myFileName = [cd '\Envelope\figure\' FileName];
        print('-djpeg',myFileName,'-r100');
        close
    end
end
save([cd '\Envelope\Envelope_1band_Hilbert_8HzLP_DS64Hz.mat'],'Envelope')

%%
% Low pass
clear

Vowel = ['a','e','i','o','u','v'];
Tone = ['1','2','3','4'];
order = 200;
window = 4;% hamming
LP = 25;% Hz
clear Envelope
for i1=1:length(Vowel)
    for i2 = 1:length(Tone)
        StiPath = [cd '\record\' Vowel(i1) Tone(i2) '_p.wav'];
        [Sti,fs] = audioread(StiPath);
        hn = usefir1(2,order,LP,[],window,fs);%2 = LP
        StiFilt = conv(hn,Sti);
        Envelope_temp = abs(hilbert(StiFilt));
        Envelope((i1-1)*4+i2,:) = resample(Envelope_temp,80,length(Envelope_temp));% downsample to 80 samples
        %         figure
        %         plot(1:length(Sti),Sti)
        %         hold on
        %         plot(1:size(Envelope,2),Envelope((i1-1)*4+i2,:));
        %         titleName = [Vowel(i1) Tone(i2)];
        %         title(titleName)
    end
end
save([cd '\Envelope\Envelope_25Hz_LP.mat'],'Envelope')

%% calculate similarity across all Tones
clear
% load([cd '\Envelope\Envelope_25Hz_LP.mat']);
% load([cd '\Envelope\Envelope_25Hz_LP.mat']);
% load([cd '\Envelope\Envelope_25Hz_LP.mat']);
% load([cd '\Envelope\Envelope_8band_Hilbert.mat']);% downsample to 80 samples
% load([cd '\Envelope\Envelope_1band_Hilbert.mat']);% downsample to 80 samples
load([cd '\Envelope\Envelope_1band_Hilbert_8HzLP_DS64Hz.mat']);% downsample to 64 Hz
Envelope = Envelope*100;% 原始数值太小

ToneXY = [1 2;1 3;1 4;2 3;2 4;3 4;1 1;2 2;3 3;4 4];
for i0=1:size(ToneXY,1)
    clear dis_mse dis_fd
    ToneX = ToneXY(i0,1);
    ToneY = ToneXY(i0,2);
    i3 = 1;
    if ToneX~=ToneY
        for i1 = ToneX:4:24
            for i2=ToneY:4:24
                dis_mse(i3) = mse(Envelope(i1,:),Envelope(i2,:));
                dis_fd(i3) = frechet((1:size(Envelope,2))',Envelope(i1,:)',(1:size(Envelope,2))',Envelope(i2,:)');
                dis_corr(i3) = corr(Envelope(i1,:)',Envelope(i2,:)');
                dis_Eucli(i3) = mean((Envelope(i1,:) - Envelope(i2,:)));
                i3 = i3+1;
            end
        end
    else % ToneX=ToneY
        for i1 = ToneX:4:24
            for i2=i1+4:4:24
                dis_mse(i3) = mse(Envelope(i1,:),Envelope(i2,:));
                dis_fd(i3) = frechet((1:size(Envelope,2))',Envelope(i1,:)',(1:size(Envelope,2))',Envelope(i2,:)');
                dis_corr(i3) = corr(Envelope(i1,:)',Envelope(i2,:)');
                dis_Eucli(i3) = mean((Envelope(i1,:) - Envelope(i2,:)));
                i3 = i3+1;
            end
        end
    end
    dist_mse(i0) = mean(dis_mse);
    dist_fd(i0) = mean(dis_fd);
    dist_corr(i0) = mean(dis_corr);
    dist_Eucli(i0) = mean(dis_Eucli);
end

