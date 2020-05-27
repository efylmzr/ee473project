[y, Fs] = audioread('ilknota.wav');
chsz = 2048;
chcnt = ceil(length(y)/chsz);
F0 = 103.83;

vals = zeros(10, chcnt);
for k = 0:chcnt-1
    if k == chcnt-1
        ch = y((chsz*k+1):end);
    else
        ch = y((chsz*k+1):chsz*(k+1));
    end
    vals(:,k+1) = get_harmonics(ch, Fs, F0, 16*chsz, 10);
end
figure, plot((1/Fs)*chsz*(0:chcnt-1), (vals')/350);
title('Chunk Method/Guitar, f_{note} = 103.83 Hz (G^#_2)');
xlabel('t (s)');
ylabel('a_n(t)');
labels = {};
for n = 1:10
    labels{end+1} = ['n = ', num2str(n)];
end
legend(labels);

audio = zeros(1, chcnt*chsz);
Ts = (1/Fs)*(0:length(audio)-1);
for n=1:10
    snwave = sin(2*pi*Ts*n*F0);
    pp = interp1((1/Fs)*chsz*(0:chcnt-1), vals(n,:), Ts, 'linear');
    %figure, plot((1/Fs)*chsz*(0:chcnt-1), vals(n,:), 'b-', Ts, pp, 'r--');
    audio = audio + snwave.*pp;
end
figure, plot(Ts, audio/500);
figure, plot((1/Fs)*(0:length(y)-1), y);
sound(audio/125, Fs);
%audiowrite('ilknota_5.wav', audio/500, Fs);

function [harmonics,phase] = get_harmonics(chunk, Fs, F0, snum, hnum)
    harmonics = zeros(1,hnum);
    phase = zeros(1,hnum);
    [Y,f] = freqz(chunk, 1, snum, Fs);
    [pks,locs] = findpeaks(abs(Y),f);
    locs0 = round(locs/F0);
    for n=1:hnum
        if isempty(pks(locs0 == n))
            harmonics(1,n) = 0;
        else
            harmonics(1,n) = max(pks(locs0 == n));
        end
    end
end
