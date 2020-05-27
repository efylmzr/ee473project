[y, Fs] = audioread('ilknota.wav');
L0 = length(y);
y = lpfilter(y, Fs/32, Fs*1.5/32, Fs, 1000);
y = downsample(y,16);

Fsn = Fs/16;
F0 = 103.83;
Nhar = 10;
L = length(y);
Ts = (1/Fsn)*(0:L-1);
fftL = 2^ceil(log2(L));

harmonics = zeros(1, Nhar);
[Y,fs] = freqz(y, 1, 2^ceil(log2(L)+1), Fsn);
[pks,locs] = findpeaks(abs(Y),fs);
locs0 = round(locs/F0);
for n=1:Nhar
    if isempty(pks(locs0 == n))
        harmonics(1,n) = n*F0;
    else
        [~,I] = max(pks(locs0 == n));
        temp = locs(locs0 == n);
        harmonics(1,n) = temp(I);
    end
end

avals = zeros(Nhar, L);
for n = 1:Nhar
    ys = cos(2*pi*harmonics(1,n)*Ts) .* (y');
    out = lpfilter(ys,F0/2,F0*1.6/2,Fsn,1000);
    avals(n,:) = abs(out);
end

figure, plot(Ts, 10*avals');
title('Demodulation Method/Guitar, f_{note} = 103.83 Hz (G^#_2)');
xlabel('t (s)');
ylabel('a_n(t)');
labels = {};
for n = 1:Nhar
    labels{end+1} = ['n = ', num2str(n)];
end
legend(labels);

audio = zeros(1,L);
for n = 1:Nhar
    audio = audio + cos(2*pi*n*F0*Ts) .* avals(n,:);
end
sound(12*audio, Fsn);

%yout = interp1(Ts, audio, (0:L0-1)*(1/Fs), 'linear');
% yout = upsample(audio, 16);
% yout = lpfilter(yout, Fsn/2, Fsn*1.5/2, Fs, 1000);
% sound(128*yout, Fs);

function out = lpfilter(in, Fc1, Fc2, Fs, N)
    f = [0 (2*Fc1/Fs) (2*Fc2/Fs) 1];
    a = [1 1 0 0];
    b = firpm(N,f,a);
    %[Yf,ff] = freqz(b,1,2^ceil(log2(length(in))), Fs);
    %figure, plot(ff, abs(Yf));
    %title('some lowpass');
    out = filter(b,1,in);
end