clear;
clc;
close all;
rng(1);

%% ========================================================================
%  GPS L1 C/A - Acquisition -> DLL/FLL/PLL -> PRN bazli korelasyon
%              -> MUSIC -> guc siniflandirma -> MVDR
%
%  Ana fikir:
%    1) PRN kodu alici tarafinda yerel olarak uretilir.
%    2) Acquisition kod fazi + Doppler adaylarini bulur.
%    3) Her aday AYRI DLL/FLL/PLL kanalinda takip edilir.
%    4) Referans antendeki takip tahminleri tum 4 antene ORTAK uygulanir.
%       Boylece antenler arasi uzaysal faz farki korunur.
%    5) Her takip dali MUSIC'e tek-kaynak olarak verilir (NumSignals = 1).
%    6) PRN1'in iki dali guce gore siralanir.
%       Dusuk guclu dal "authentic", yuksek guclu dal "spoofer" kabul edilir.
%    7) Authentic MUSIC AoA -> MVDR Direction.
%       Spoofer dali -> MVDR training kovaryansi (null olusmasini kolaylastirir).
%% ========================================================================

%% Temel parametreler
numbits = 25;                % 25 x 20 ms = 500 ms

% Her kaynak farkli navigation data tasir. Boylece MVDR/NLMS sonrasinda
% hangi navigation verisinin geri kazanildigini BER ile dogrulayabiliriz.
navdata_1 = randi([0 1], numbits, 1);
navdata_2 = randi([0 1], numbits, 1);
navdata_spf = randi([0 1], numbits, 1);

% Tesadufen ayni veya tam tersi veri dizisi uretilmesini engelle.
while isequal(navdata_2,navdata_1) || isequal(navdata_2,1-navdata_1)
    navdata_2 = randi([0 1], numbits, 1);
end
while isequal(navdata_spf,navdata_1) || isequal(navdata_spf,1-navdata_1) || ...
      isequal(navdata_spf,navdata_2) || isequal(navdata_spf,1-navdata_2)
    navdata_spf = randi([0 1], numbits, 1);
end

fs     = 16.368e6;
fIF    = 4.092e6;
fL1    = 1575.42e6;
fChip  = 1.023e6;
sigma2 = 1;
c      = physconst("LightSpeed");

N_1ms      = round(fs * 1e-3);
ornPerChip = fs / fChip;     % 16 sample/chip
assert(abs(ornPerChip-round(ornPerChip)) < 1e-12, ...
    'Bu kod fs/fChip tam sayi olacak sekilde yazildi.');

%% Gercek uydu 1 - PRN1
prn1   = 1;
sat1Az = 30;
sat1El = 50;
sat1Fd = 1500;
sat1CN0 = 45;

%% Gercek uydu 2 - PRN2
prn2   = 2;
sat2Az = 120;
sat2El = 70;
sat2Fd = -500;
sat2CN0 = 44;

%% Spoofer - PRN1 ile ayni PRN
spfPRN      = 1;
spfAz       = 45;
spfEl       = 20;
spfFdFark   = 30;            % 1530 Hz toplam Doppler
spfKodKayma = 3;             % +3 chip
spfFaz      = 0;
spfGucFark  = 10;            % gercek PRN1'den +10 dB

%% Dalga bicimleri
wave1 = gpsWaveformGenerator( ...
    SignalType="legacy", PRNID=prn1, EnablePCode=false, SampleRate=fs);
baseband1 = wave1(navdata_1);

wave2 = gpsWaveformGenerator( ...
    SignalType="legacy", PRNID=prn2, EnablePCode=false, SampleRate=fs);
baseband2 = wave2(navdata_2);

waveSpf = gpsWaveformGenerator( ...
    SignalType="legacy", PRNID=spfPRN, EnablePCode=false, SampleRate=fs);
basebandSpf = waveSpf(navdata_spf);

t = (0:length(baseband1)-1).' / fs;
num_epochs = floor(length(baseband1) / N_1ms);

sat1_if = baseband1 .* exp(1j*2*pi*(fIF + sat1Fd)*t);
sat2_if = baseband2 .* exp(1j*2*pi*(fIF + sat2Fd)*t);

kaymaOrnek = round(spfKodKayma * ornPerChip);
spf_bb = circshift(basebandSpf, kaymaOrnek);
spf_if = spf_bb .* exp(1j*(2*pi*(fIF + sat1Fd + spfFdFark)*t + spfFaz));

%% 4 elemanli CRPA
lambda = c/fL1;
dizi = phased.URA(Size=[2 2], ElementSpacing=lambda/2, ArrayNormal="z");

Ps = mean(abs(baseband1).^2);
genlik = @(CN0dB) sqrt(10^(CN0dB/10) * sigma2 / (fs * Ps));

Xsat1 = collectPlaneWave(dizi, genlik(sat1CN0)*sat1_if, [sat1Az; sat1El], fL1);
Xsat2 = collectPlaneWave(dizi, genlik(sat2CN0)*sat2_if, [sat2Az; sat2El], fL1);

gurultu = sqrt(sigma2/2) * ...
    (randn(size(Xsat1)) + 1j*randn(size(Xsat1)));

Xspf = collectPlaneWave(dizi, ...
    genlik(sat1CN0 + spfGucFark)*spf_if, [spfAz; spfEl], fL1);
X = Xsat1 + Xsat2 + Xspf + gurultu;
fprintf('\n===== SENARYO: PRN1 + PRN2 + PRN1 SPOOFER =====\n');

%% Yerel C/A kodlari - alici PRN'i bilir, kod fazini acquisition bulur
ca1 = 1 - 2*double(gnssCACode(prn1,"GPS"));
ca2 = 1 - 2*double(gnssCACode(prn2,"GPS"));

kod1 = repelem(ca1(:), round(ornPerChip));
kod2 = repelem(ca2(:), round(ornPerChip));
kod1 = kod1(1:N_1ms);
kod2 = kod2(1:N_1ms);

%% Acquisition ayarlari
% 1 ms koherent + nNC epoch non-coherent toplama.
% PRN1 icin iki AYRI kod-fazi/Doppler adayi araniyor: gercek + spoofer.
nNC = 10;
fdAra1 = sat1Fd + (-100:5:100);
fdAra2 = sat2Fd + (-100:5:100);

nAday1 = 2;   % PRN1 icin: gercek uydu + ayni PRN'i kullanan spoofer

% Ayni korelasyon tepesinin komsu hucrelerini ikinci aday diye secmemek icin
% guard bolgesi. Gercek/spoof ayrimi 3 chip oldugu icin 1.25 chip yeterli.
guardChip = 1.25;

% Acquisition referans anten = anten 1
[tau1Acq, fd1Acq, acqPow1, S1] = kodEdinmeCoklu( ...
    X(:,1), kod1, N_1ms, fs, fIF, fdAra1, nNC, nAday1, ...
    guardChip, ornPerChip);

[tau2Acq, fd2Acq, acqPow2, S2] = kodEdinmeCoklu( ...
    X(:,1), kod2, N_1ms, fs, fIF, fdAra2, nNC, 1, ...
    guardChip, ornPerChip);

fprintf('\n--- ACQUISITION ---\n');
for q = 1:nAday1
    fprintf('PRN1 aday %d: kod = %.3f chip, Doppler = %.1f Hz, relatif tepe = %.2f dB\n', ...
        q, tau1Acq(q)/ornPerChip, fd1Acq(q), ...
        10*log10(acqPow1(q)/max(acqPow1)));
end
fprintf('PRN2       : kod = %.3f chip, Doppler = %.1f Hz\n', ...
    tau2Acq(1)/ornPerChip, fd2Acq(1));

%% DLL/FLL/PLL parametreleri
prm.dllBn = 1;     % Hz
prm.pllBn = 30;    % Hz
prm.fllBn = 4;     % Hz

% Her PRN1 adayi AYRI tracker ile izlenir.
takip1 = cell(nAday1,1);
for q = 1:nAday1
    takip1{q} = takipAdayToolbox( ...
        X(:,1), fs, fIF, prn1, N_1ms, num_epochs, ...
        tau1Acq(q), fd1Acq(q), ornPerChip, prm);
end

takip2 = takipAdayToolbox( ...
    X(:,1), fs, fIF, prn2, N_1ms, num_epochs, ...
    tau2Acq(1), fd2Acq(1), ornPerChip, prm);

%% Steady-state kontrolu
steadyMs  = min(200,num_epochs);
steadyBas = max(51, num_epochs-steadyMs+1);   % ilk 50 ms'i en azindan at
steadyIdx = steadyBas:num_epochs;

fprintf('\n--- DLL/FLL/PLL STEADY-STATE ---\n');
for q = 1:nAday1
    fprintf(['PRN1 aday %d: kod = %.3f +/- %.4f chip | ', ...
             'Doppler = %.2f +/- %.2f Hz | lock = %.3f\n'], ...
        q, mean(takip1{q}.tauChip(steadyIdx)), std(takip1{q}.tauChip(steadyIdx)), ...
        mean(takip1{q}.fd(steadyIdx)), std(takip1{q}.fd(steadyIdx)), ...
        mean(takip1{q}.lock(steadyIdx)));
end
fprintf(['PRN2       : kod = %.3f +/- %.4f chip | ', ...
         'Doppler = %.2f +/- %.2f Hz | lock = %.3f\n'], ...
    mean(takip2.tauChip(steadyIdx)), std(takip2.tauChip(steadyIdx)), ...
    mean(takip2.fd(steadyIdx)), std(takip2.fd(steadyIdx)), ...
    mean(takip2.lock(steadyIdx)));

%% Takip bilgisi ile 4 antenin ORTAK despreading/korelasyonu
% Kritik nokta: DLL/PLL sadece referans antende tahmin uretiyor.
% Ayni tau/fd tum antenlere uygulanir; antenler arasi faz farki korunur.
C1 = cell(nAday1,1);
for q = 1:nAday1
    C1{q} = uzaysalPromptKorelasyon( ...
        X, kod1, N_1ms, num_epochs, fs, fIF, ...
        takip1{q}.tauChip, takip1{q}.fd, ornPerChip);
end

C2 = uzaysalPromptKorelasyon( ...
    X, kod2, N_1ms, num_epochs, fs, fIF, ...
    takip2.tauChip, takip2.fd, ornPerChip);

%% MUSIC - her takip dali zaten tek kod/Doppler adayi oldugu icin NumSignals=1
azScan = -180:1:180;
elScan = 0:1:90;

musicEst = phased.MUSICEstimator2D( ...
    SensorArray=dizi, OperatingFrequency=fL1, ...
    AzimuthScanAngles=azScan, ElevationScanAngles=elScan, ...
    DOAOutputPort=true, NumSignalsSource="Property", NumSignals=1);

spec1 = cell(nAday1,1);
doa1  = zeros(2,nAday1);
pow1  = zeros(nAday1,1);

for q = 1:nAday1
    [spec1{q}, a] = musicEst(C1{q}(steadyIdx,:));
    doa1(:,q) = a(:,1);
    pow1(q) = mean(abs(C1{q}(steadyIdx,:)).^2, 'all');
end

[spec2, doa2] = musicEst(C2(steadyIdx,:));
doa2 = doa2(:,1);
pow2 = mean(abs(C2(steadyIdx,:)).^2, 'all');

fprintf('\n--- MUSIC ---\n');
for q = 1:nAday1
    fprintf('PRN1 aday %d: AoA = [%.1f deg, %.1f deg], korelasyon gucu = %.2f dB\n', ...
        q, doa1(1,q), doa1(2,q), 10*log10(pow1(q)));
end
fprintf('PRN2       : AoA = [%.1f deg, %.1f deg], korelasyon gucu = %.2f dB\n', ...
    doa2(1), doa2(2), 10*log10(pow2));

%% PRN1 guc siniflandirma + MVDR
[~, idxWeak]   = min(pow1);
[~, idxStrong] = max(pow1);

doaAuthentic = doa1(:,idxWeak);
doaSpoofer   = doa1(:,idxStrong);

fprintf('\n--- PRN1 SINIFLANDIRMA ---\n');
fprintf('Dusuk guclu aday -> AUTHENTIC: aday %d, AoA=[%.1f, %.1f]\n', ...
    idxWeak, doaAuthentic(1), doaAuthentic(2));
fprintf('Yuksek guclu aday -> SPOOFER  : aday %d, AoA=[%.1f, %.1f]\n', ...
    idxStrong, doaSpoofer(1), doaSpoofer(2));
fprintf('Olculen guc farki (strong/weak) = %.2f dB\n', ...
    10*log10(pow1(idxStrong)/pow1(idxWeak)));

% MVDR Direction = MUSIC'in dusuk guclu/authentic PRN1 AoA'si.
% Training = guclu/spoofer takip dalinin uzaysal snapshot'lari.
Xtrain = C1{idxStrong}(steadyIdx,:);
Xtrain = Xtrain / sqrt(mean(abs(Xtrain).^2,'all'));

mvdr1 = phased.MVDRBeamformer( ...
    SensorArray=dizi, OperatingFrequency=fL1, ...
    Direction=doaAuthentic, ...
    TrainingInputPort=true, WeightsOutputPort=true, ...
    DiagonalLoadingFactor=1e-3);

[yMVDR1, wMVDR1] = mvdr1(X, Xtrain);

% PRN2 icin ayri MVDR agirliklari. Huzme PRN2 MUSIC yonune bakar;
% ayni spoofer takip dali training olarak kullanilarak spoofer yonu bastirilir.
mvdr2 = phased.MVDRBeamformer( ...
    SensorArray=dizi, OperatingFrequency=fL1, ...
    Direction=doa2, ...
    TrainingInputPort=true, WeightsOutputPort=true, ...
    DiagonalLoadingFactor=1e-3);

[yMVDR2, wMVDR2] = mvdr2(X, Xtrain);


%% MVDR SONRASI TEKRAR DLL/FLL/PLL TAKIBI
% Yeni akis:
%   MVDR -> DLL/FLL/PLL -> Prompt korelator -> NLMS -> Navigation bit karari
%
% Boylece NLMS artik ham IF tasiyicisini ve C/A kodunu degistirmez.
% NLMS, DLL/PLL ve despreading sonrasindaki 1 ms'lik kompleks Prompt
% korelator cikislarini isler.

% PRN1 authentic dalinin acquisition baslangici
tauAuthAcq = tau1Acq(idxWeak);
fdAuthAcq  = fd1Acq(idxWeak);

% MVDR cikislarini hedef PRN tracker'larina tekrar ver.
trkMVDR1 = takipAdayToolbox( ...
    yMVDR1, fs, fIF, prn1, N_1ms, num_epochs, ...
    tauAuthAcq, fdAuthAcq, ornPerChip, prm);

trkMVDR2 = takipAdayToolbox( ...
    yMVDR2, fs, fIF, prn2, N_1ms, num_epochs, ...
    tau2Acq(1), fd2Acq(1), ornPerChip, prm);

Pmvdr1 = trkMVDR1.P;   % 1 kompleks Prompt / ms
Pmvdr2 = trkMVDR2.P;

lockMVDR1 = mean(trkMVDR1.lock(steadyIdx));
lockMVDR2 = mean(trkMVDR2.lock(steadyIdx));

fprintf('\n--- MVDR SONRASI DLL/FLL/PLL ---\n');
fprintf('PRN1 tracker lock = %.3f | son Doppler = %.2f Hz | son kod = %.3f chip\n', ...
    lockMVDR1, mean(trkMVDR1.fd(steadyIdx)), mean(trkMVDR1.tauChip(steadyIdx)));
fprintf('PRN2 tracker lock = %.3f | son Doppler = %.2f Hz | son kod = %.3f chip\n', ...
    lockMVDR2, mean(trkMVDR2.fd(steadyIdx)), mean(trkMVDR2.tauChip(steadyIdx)));

%% NLMS - TRACKING / DESPREADING SONRASI PROMPT CIKISLARINA
% Prompt ornekleme hizi 1 kHz'tir (1 Prompt / 1 ms).
% 1 epoch gecikmeli ALE kullanilir. Bu noktada C/A kodu ve tasiyici
% tracking ile giderildigi icin NLMS'nin navigation verisini bozma riski
% ham IF uzerinde calismaya gore daha dusuktur.

nlmsDelayEp = 1;     % 1 ms = 1 Prompt epoch
nlmsOrder   = 16;    % daha kisa ve daha kararlı adaptif filtre

% --- PRN1 Prompt -> NLMS ---
xNLMS1 = [complex(0); Pmvdr1(1:end-nlmsDelayEp)];
nlms1 = dsp.LMSFilter(nlmsOrder, Method="Normalized LMS");
nStep1 = min(length(xNLMS1), 200);
muMax1 = maxstep(nlms1, xNLMS1(1:nStep1));
nlms1.StepSize = muMax1/20;
[Pnlms1, eNLMS1, wNLMS1] = nlms1(xNLMS1, Pmvdr1);

% --- PRN2 Prompt -> NLMS ---
xNLMS2 = [complex(0); Pmvdr2(1:end-nlmsDelayEp)];
nlms2 = dsp.LMSFilter(nlmsOrder, Method="Normalized LMS");
nStep2 = min(length(xNLMS2), 200);
muMax2 = maxstep(nlms2, xNLMS2(1:nStep2));
nlms2.StepSize = muMax2/20;
[Pnlms2, eNLMS2, wNLMS2] = nlms2(xNLMS2, Pmvdr2);

fprintf('\n--- NLMS (TRACKING SONRASI PROMPT DOMAIN) ---\n');
fprintf('PRN1 Prompt -> NLMS: D = %d ms, M = %d, mu = %.4g\n', ...
    nlmsDelayEp, nlmsOrder, nlms1.StepSize);
fprintf('PRN2 Prompt -> NLMS: D = %d ms, M = %d, mu = %.4g\n', ...
    nlmsDelayEp, nlmsOrder, nlms2.StepSize);

%% NAVIGATION DATA DOGRULAMA
% MVDR satirlari tracker Prompt cikisindan;
% MVDR+NLMS satirlari ayni Prompt'larin NLMS sonrasindan decode edilir.
% NLMS'den sonra tekrar PLL calistirilmadigi icin PLL lock degeri,
% NLMS ONCESINDEKI (MVDR cikisindaki) tracker lock degeridir.

steadyBitStart = max(1, ceil(steadyBas/20));

[bitsMVDR1, metricMVDR1, berMVDR1, invMVDR1] = ...
    navBitDogrula(Pmvdr1, navdata_1, steadyBitStart);
[bitsNLMS1, metricNLMS1, berNLMS1, invNLMS1] = ...
    navBitDogrula(Pnlms1, navdata_1, steadyBitStart);

[bitsMVDR2, metricMVDR2, berMVDR2, invMVDR2] = ...
    navBitDogrula(Pmvdr2, navdata_2, steadyBitStart);
[bitsNLMS2, metricNLMS2, berNLMS2, invNLMS2] = ...
    navBitDogrula(Pnlms2, navdata_2, steadyBitStart);

% PRN1 cikisini hem authentic hem spoofer navigation datasiyla karsilastir.
berMVDR1_spf = navBERKarsilastir(bitsMVDR1, navdata_spf, steadyBitStart);
berNLMS1_spf = navBERKarsilastir(bitsNLMS1, navdata_spf, steadyBitStart);

fprintf('\n--- NAVIGATION DATA DOGRULAMA (steady-state bitleri) ---\n');
fprintf('PRN1 MVDR Prompt      : BER(auth)=%.3f | BER(spf)=%.3f | pre-NLMS lock=%.3f | flip=%d\n', ...
    berMVDR1, berMVDR1_spf, lockMVDR1, invMVDR1);
fprintf('PRN1 Prompt + NLMS    : BER(auth)=%.3f | BER(spf)=%.3f | pre-NLMS lock=%.3f | flip=%d\n', ...
    berNLMS1, berNLMS1_spf, lockMVDR1, invNLMS1);
fprintf('PRN2 MVDR Prompt      : BER=%.3f | pre-NLMS lock=%.3f | flip=%d\n', ...
    berMVDR2, lockMVDR2, invMVDR2);
fprintf('PRN2 Prompt + NLMS    : BER=%.3f | pre-NLMS lock=%.3f | flip=%d\n', ...
    berNLMS2, lockMVDR2, invNLMS2);

if lockMVDR1 < 0.3
    fprintf('UYARI: PRN1 MVDR cikisinda PLL lock zayif; NLMS BER sonucu guvenilir olmayabilir.\n');
elseif berNLMS1 < berNLMS1_spf
    fprintf('PRN1 NLMS decoded veri AUTHENTIC navigation datasina daha yakin.\n');
elseif berNLMS1_spf < berNLMS1
    fprintf('PRN1 NLMS decoded veri SPOOFER navigation datasina daha yakin.\n');
else
    fprintf('PRN1 NLMS decoded veri authentic/spoofer arasinda ayirt edici degil.\n');
end

Tnav = table( ...
    ["PRN1 MVDR Prompt";"PRN1 Prompt+NLMS";"PRN2 MVDR Prompt";"PRN2 Prompt+NLMS"], ...
    [berMVDR1;berNLMS1;berMVDR2;berNLMS2], ...
    [berMVDR1_spf;berNLMS1_spf;NaN;NaN], ...
    [lockMVDR1;lockMVDR1;lockMVDR2;lockMVDR2], ...
    [invMVDR1;invNLMS1;invMVDR2;invNLMS2], ...
    'VariableNames',{'Cikis','BER_Authentic','BER_Spoofer','PLL_Lock_PreNLMS','PolarityFlip'});
disp(Tnav);

sv = phased.SteeringVector(SensorArray=dizi, PropagationSpeed=c);
aAuth = sv(fL1, doaAuthentic);
aPrn2 = sv(fL1, doa2);
aSpf  = sv(fL1, doaSpoofer);

gAuth = abs(wMVDR1' * aAuth)^2;
gSpf1 = abs(wMVDR1' * aSpf)^2;
sonum1 = 10*log10(gAuth/gSpf1);

gPrn2 = abs(wMVDR2' * aPrn2)^2;
gSpf2 = abs(wMVDR2' * aSpf)^2;
sonum2 = 10*log10(gPrn2/gSpf2);

fprintf('\n--- MVDR1 ---\n');
fprintf('Direction = [%.1f, %.1f] deg (AUTHENTIC PRN1)\n', ...
    doaAuthentic(1), doaAuthentic(2));
fprintf('Spoofer yonune gore relatif sonum = %.2f dB\n', sonum1);

fprintf('\n--- MVDR2 ---\n');
fprintf('Direction = [%.1f, %.1f] deg (PRN2)\n', doa2(1), doa2(2));
fprintf('Spoofer yonune gore relatif sonum = %.2f dB\n', sonum2);

%% ========================================================================
% PLOTLAR
%% ========================================================================

%% 1) DLL/FLL/PLL takip grafikleri
figure('Name','DLL-FLL-PLL Takip','Position',[80 80 1350 650]);
tl = tiledlayout(2,1,'TileSpacing','compact','Padding','compact');

nexttile;
for q = 1:nAday1
    plot(1:num_epochs, takip1{q}.tauChip, 'LineWidth',1.2, ...
        'DisplayName',sprintf('PRN1 aday %d',q));
    hold on;
end
plot(1:num_epochs, takip2.tauChip, 'LineWidth',1.2, 'DisplayName','PRN2');
yline(0,'--','PRN1 gercek','HandleVisibility','off');
yline(spfKodKayma,'--','PRN1 spoofer','HandleVisibility','off');
xline(steadyBas,':','steady-state','HandleVisibility','off');
grid on; legend('Location','best');
xlabel('Zaman [ms]'); ylabel('Kod fazi [chip]');
title('DLL: Kod fazi takibi');

nexttile;
for q = 1:nAday1
    plot(1:num_epochs, takip1{q}.fd, 'LineWidth',1.2, ...
        'DisplayName',sprintf('PRN1 aday %d',q));
    hold on;
end
plot(1:num_epochs, takip2.fd, 'LineWidth',1.2, 'DisplayName','PRN2');
yline(sat1Fd,'--','PRN1 gercek','HandleVisibility','off');
yline(sat2Fd,'--','PRN2 gercek','HandleVisibility','off');
yline(sat1Fd+spfFdFark,'--','PRN1 spoofer','HandleVisibility','off');
xline(steadyBas,':','steady-state','HandleVisibility','off');
grid on; legend('Location','best');
xlabel('Zaman [ms]'); ylabel('Doppler [Hz]');
title('FLL/PLL: Doppler takibi');
title(tl, 'Tracking - PRN1 + PRN2 + PRN1 Spoofer');

%% 2) MUSIC spektrumlari - TEK BIRLESIK GRAFIK
% Her kod/Doppler dali MUSIC'te ayri ayri islenir (NumSignals=1), ancak
% gorsellestirme icin her dalin kendi tepesine normalize edilmis MUSIC
% spektrumu tek bir haritada birlestirilir. Boylece zayif authentic dal da
% guclu spoofer dali tarafindan gorsel olarak bastirilmaz.

nSpec = nAday1 + 1;
specStack = zeros(numel(elScan), numel(azScan), nSpec);

for q = 1:nAday1
    specStack(:,:,q) = musicDb(spec1{q}, numel(elScan));
end
specStack(:,:,end) = musicDb(spec2, numel(elScan));

% Her dal zaten 0 dB'e normalize. Noktasal maksimum, tum dallardaki
% MUSIC tepelerini tek uzaysal spektrumda korur.
sdbAll = max(specStack, [], 3);

figure('Name','Birlesik MUSIC Spektrumu','Position',[140 100 950 620]);
imagesc(azScan, elScan, sdbAll);
axis xy; colorbar; caxis([-35 0]); hold on;

% MUSIC kestirimleri
for q = 1:nAday1
    if q == idxStrong
        plot(doa1(1,q), doa1(2,q), 'rx', 'MarkerSize',13, 'LineWidth',2.2, ...
            'DisplayName',sprintf('PRN1 aday %d - Spoofer',q));
    else
        plot(doa1(1,q), doa1(2,q), 'gx', 'MarkerSize',13, 'LineWidth',2.2, ...
            'DisplayName',sprintf('PRN1 aday %d - Authentic',q));
    end
end
plot(doa2(1), doa2(2), 'mx', 'MarkerSize',13, 'LineWidth',2.2, ...
    'DisplayName','PRN2 MUSIC');

% Simulasyonda bilinen gercek yonler (yalnizca dogrulama icin)
plot(sat1Az, sat1El, 'ko', 'MarkerSize',9, 'LineWidth',1.4, ...
    'DisplayName','PRN1 gercek yon');
plot(spfAz, spfEl, 'ks', 'MarkerSize',9, 'LineWidth',1.4, ...
    'DisplayName','Spoofer gercek yon');
plot(sat2Az, sat2El, 'kd', 'MarkerSize',9, 'LineWidth',1.4, ...
    'DisplayName','PRN2 gercek yon');

xlabel('Azimuth [deg]');
ylabel('Elevation [deg]');
title('Kod/Doppler Dallarindan Birlesik MUSIC Uzaysal Spektrumu');
legend('Location','best');
grid on;

%% 3) MVDR1 + MVDR2 3B beamforming patternlari - tek figure
PAT1 = pattern(dizi,fL1,azScan,elScan, ...
    PropagationSpeed=c,Weights=wMVDR1,Type="powerdb",Normalize=true);
PAT2 = pattern(dizi,fL1,azScan,elScan, ...
    PropagationSpeed=c,Weights=wMVDR2,Type="powerdb",Normalize=true);

if size(PAT1,1) ~= numel(elScan), PAT1 = PAT1.'; end
if size(PAT2,1) ~= numel(elScan), PAT2 = PAT2.'; end

% Cok derin null degerlerini gorsellestirme icin -50 dB'de kirp.
PAT1 = max(PAT1,-50);
PAT2 = max(PAT2,-50);

[AZ,EL] = meshgrid(azScan,elScan);

% Isaretlenecek yonlerin 3B yuzeydeki guc degerleri.
[~,iAzAuth] = min(abs(azScan-doaAuthentic(1)));
[~,iElAuth] = min(abs(elScan-doaAuthentic(2)));
[~,iAzSpf]  = min(abs(azScan-doaSpoofer(1)));
[~,iElSpf]  = min(abs(elScan-doaSpoofer(2)));
[~,iAz2]    = min(abs(azScan-doa2(1)));
[~,iEl2]    = min(abs(elScan-doa2(2)));

zAuth1 = PAT1(iElAuth,iAzAuth);
zSpf1  = PAT1(iElSpf,iAzSpf);
zPrn2  = PAT2(iEl2,iAz2);
zSpf2  = PAT2(iElSpf,iAzSpf);

figure('Name','MVDR 3D Beamforming Patterns','Position',[60 120 1500 620]);
tlMVDR = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

% --- MVDR1 / PRN1 ---
nexttile;
surf(AZ,EL,PAT1,'EdgeColor','none'); hold on;
plot3(doaAuthentic(1),doaAuthentic(2),zAuth1,'go', ...
    'MarkerSize',10,'LineWidth',2,'MarkerFaceColor','g');
plot3(doaSpoofer(1),doaSpoofer(2),zSpf1,'rs', ...
    'MarkerSize',10,'LineWidth',2,'MarkerFaceColor','r');
view(45,35); grid on; box on;
xlabel('Azimuth [deg]'); ylabel('Elevation [deg]'); zlabel('Normalize Guc [dB]');
xlim([azScan(1) azScan(end)]); ylim([elScan(1) elScan(end)]); zlim([-50 0]);
clim([-50 0]); colorbar;
legend({'MVDR1 pattern','PRN1 Authentic','Spoofer'},'Location','best');
title(sprintf('MVDR1 | PRN1=[%.0f, %.0f] | Spoofer sonum %.1f dB', ...
    doaAuthentic(1),doaAuthentic(2),sonum1));

% --- MVDR2 / PRN2 ---
nexttile;
surf(AZ,EL,PAT2,'EdgeColor','none'); hold on;
plot3(doa2(1),doa2(2),zPrn2,'co', ...
    'MarkerSize',10,'LineWidth',2,'MarkerFaceColor','c');
plot3(doaSpoofer(1),doaSpoofer(2),zSpf2,'rs', ...
    'MarkerSize',10,'LineWidth',2,'MarkerFaceColor','r');
view(45,35); grid on; box on;
xlabel('Azimuth [deg]'); ylabel('Elevation [deg]'); zlabel('Normalize Guc [dB]');
xlim([azScan(1) azScan(end)]); ylim([elScan(1) elScan(end)]); zlim([-50 0]);
clim([-50 0]); colorbar;
legend({'MVDR2 pattern','PRN2','Spoofer'},'Location','best');
title(sprintf('MVDR2 | PRN2=[%.0f, %.0f] | Spoofer sonum %.1f dB', ...
    doa2(1),doa2(2),sonum2));

title(tlMVDR,'MVDR 3B Beamforming Diyagramlari');


%% 4) MVDR -> DLL/FLL/PLL PROMPT ve NLMS cikislari
% NLMS artik ham IF orneklerinde degil, 1 ms'lik Prompt korelator
% cikislarinda calisiyor. Son 200 ms steady-state bolgesi gosterilir.

plotEp = steadyIdx;
tPlotMs = plotEp(:);

figure('Name','Tracking Prompt + NLMS','Position',[100 150 1450 560]);
tlNLMS = tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(tPlotMs, real(Pmvdr1(plotEp)), 'LineWidth',1.0, 'DisplayName','MVDR1 -> DLL/PLL Prompt'); hold on;
plot(tPlotMs, real(Pnlms1(plotEp)), 'LineWidth',1.4, 'DisplayName','Prompt + NLMS');
grid on; legend('Location','best');
xlabel('Zaman [ms]'); ylabel('Prompt I');
title(sprintf('PRN1 | pre-NLMS PLL lock=%.3f',lockMVDR1));

nexttile;
plot(tPlotMs, real(Pmvdr2(plotEp)), 'LineWidth',1.0, 'DisplayName','MVDR2 -> DLL/PLL Prompt'); hold on;
plot(tPlotMs, real(Pnlms2(plotEp)), 'LineWidth',1.4, 'DisplayName','Prompt + NLMS');
grid on; legend('Location','best');
xlabel('Zaman [ms]'); ylabel('Prompt I');
title(sprintf('PRN2 | pre-NLMS PLL lock=%.3f',lockMVDR2));

title(tlNLMS, sprintf('MVDR -> DLL/FLL/PLL -> Prompt -> NLMS | D=%d ms, M=%d', ...
    nlmsDelayEp, nlmsOrder));

%% 5) Navigation bit dogrulama grafigi
% Referans ve NLMS ile elde edilen bitler ayni eksende gosterilir.
% Sadece steady-state bitleri cizilir.
nBitPlot1 = min(numel(bitsNLMS1), numel(navdata_1));
nBitPlot2 = min(numel(bitsNLMS2), numel(navdata_2));
idxBit1 = steadyBitStart:nBitPlot1;
idxBit2 = steadyBitStart:nBitPlot2;

figure('Name','Navigation Data Dogrulama','Position',[80 180 1650 520]);
tlNav = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

nexttile;
stairs(idxBit1, navdata_1(idxBit1), 'LineWidth',1.8, 'DisplayName','PRN1 authentic'); hold on;
stairs(idxBit1, double(bitsNLMS1(idxBit1))+0.04, '--', 'LineWidth',1.5, 'DisplayName','PRN1 NLMS decoded');
ylim([-0.2 1.25]); yticks([0 1]); grid on; legend('Location','best');
xlabel('Navigation bit indeksi'); ylabel('Bit');
title(sprintf('PRN1 authentic | BER=%.3f | pre-NLMS lock=%.3f',berNLMS1,lockMVDR1));

nexttile;
stairs(idxBit1, navdata_spf(idxBit1), 'LineWidth',1.8, 'DisplayName','Spoofer navdata'); hold on;
stairs(idxBit1, double(bitsNLMS1(idxBit1))+0.04, '--', 'LineWidth',1.5, 'DisplayName','PRN1 NLMS decoded');
ylim([-0.2 1.25]); yticks([0 1]); grid on; legend('Location','best');
xlabel('Navigation bit indeksi'); ylabel('Bit');
title(sprintf('PRN1 NLMS vs Spoofer | BER=%.3f',berNLMS1_spf));

nexttile;
stairs(idxBit2, navdata_2(idxBit2), 'LineWidth',1.8, 'DisplayName','PRN2 orijinal'); hold on;
stairs(idxBit2, double(bitsNLMS2(idxBit2))+0.04, '--', 'LineWidth',1.5, 'DisplayName','PRN2 NLMS decoded');
ylim([-0.2 1.25]); yticks([0 1]); grid on; legend('Location','best');
xlabel('Navigation bit indeksi'); ylabel('Bit');
title(sprintf('PRN2 | BER=%.3f | pre-NLMS lock=%.3f',berNLMS2,lockMVDR2));

title(tlNav,'MVDR -> DLL/FLL/PLL -> Prompt -> NLMS Navigation Data Dogrulamasi');

%% ========================================================================
% YEREL FONKSIYONLAR
%% ========================================================================

function [tauList,fdList,powList,S] = kodEdinmeCoklu( ...
    y,kodPer,N,fs,fIF,fdAra,nNC,nPeaks,guardChip,ornPerChip)
% PRN acquisition: Doppler x kod-fazi korelasyon haritasi olusturur ve
% farkli kod fazlarina ait en guclu nPeaks tepeyi dondurur.

    Kf = conj(fft(kodPer));
    S = zeros(numel(fdAra),N);

    for iF = 1:numel(fdAra)
        for ep = 1:nNC
            idx = (ep-1)*N + (1:N).';
            carrier = exp(-1j*2*pi*(fIF+fdAra(iF))*(idx-1)/fs);
            yb = y(idx).*carrier;
            cc = ifft(fft(yb).*Kf);
            S(iF,:) = S(iF,:) + abs(cc.').^2;
        end
    end

    W = S;
    tauList = zeros(nPeaks,1);
    fdList  = zeros(nPeaks,1);
    powList = zeros(nPeaks,1);

    guardSamp = max(1,round(guardChip*ornPerChip));

    for p = 1:nPeaks
        [pk,lin] = max(W(:));
        if ~isfinite(pk)
            error('Yeterli acquisition tepesi bulunamadi.');
        end
        [iF,iT] = ind2sub(size(W),lin);

        tauList(p) = iT-1;
        fdList(p)  = fdAra(iF);
        powList(p) = pk;

        % Ayni fiziksel sinyalin farkli Doppler hucrelerinin ikinci aday
        % olarak secilmesini engelle: bulunan kod-fazi cevresini tum
        % Doppler ekseninde kapat.
        codeBins = mod((iT-1) + (-guardSamp:guardSamp),N) + 1;
        W(:,codeBins) = -Inf;
    end
end

function s = takipAdayToolbox(y,fs,fIF,prnID,N,E,tauAcqSamp,fdAcq,ornPerChip,prm)
% Tek bir acquisition adayi icin gnssSignalTracker.
% Takip yalnizca referans antende yapilir. Cikan tau/fd daha sonra tum CRPA
% elemanlarina ortak uygulanarak uzaysal faz bilgisi korunur.

    L = E*N;

    % CRPA sinyali kompleks analitik IF oldugu icin baseband'e indiriyoruz.
    t = (0:L-1).' / fs;
    yBB = y(1:L).*exp(-1j*2*pi*fIF*t);
    % System object giris karmasikligini (real/complex) sabit tut.
    % Ozellikle NLMS cikisinin ilk 1 ms'i sifir oldugundan MATLAB bu
    % ilk frame'i real algilayabiliyor; sonraki frame complex olunca
    % gnssSignalTracker hata veriyor.
    yBB = complex(real(yBB), imag(yBB));

    chip0 = mod(round(tauAcqSamp/ornPerChip),1023);

    gst = gnssSignalTracker( ...
        GNSSSignalType="GPS C/A", ...
        SampleRate=fs, ...
        PRNID=prnID, ...
        IntermediateFrequency=0, ...
        InitialCodePhaseOffset=chip0, ...
        InitialFrequencyOffset=fdAcq, ...
        PLLNoiseBandwidth=prm.pllBn, ...
        FLLNoiseBandwidth=prm.fllBn, ...
        DLLNoiseBandwidth=prm.dllBn, ...
        IntegrationTime=1e-3);

    s.tauChip = zeros(E,1);
    s.fd      = zeros(E,1);
    s.P       = complex(zeros(E,1));
    s.lock    = zeros(E,1);

    for ep = 1:E
        idx = (ep-1)*N + (1:N);
        frame = yBB(idx);
        frame = complex(real(frame), imag(frame));
        [P,info] = gst(frame);

        % Bu projede onceki testlerde gnssSignalTracker NCO ciktilari,
        % acquisition baslangicina gore duzeltme gibi davrandi.
        % Bu nedenle toplam tahminler acquisition degeriyle yeniden kurulur.
        fdTot  = fdAcq - info.FrequencyEstimate;
        tauTot = chip0 - info.DelayEstimate;

        s.fd(ep)      = fdTot;
        s.tauChip(ep) = mod(tauTot + 511.5,1023) - 511.5;
        s.P(ep)       = P;

        % PLL faz hatasi uzerinden basit lock gostergesi.
        s.lock(ep) = cos(2*info.PhaseError);
    end
end

function C = uzaysalPromptKorelasyon(X,kodPer,N,E,fs,fIF,tauChip,fdHz,ornPerChip)
% Referans tracking kanalinin kod/frekans tahminlerini 4 antene ORTAK uygular.
% Boylece her epoch sonunda 1xM uzaysal snapshot elde edilir ve antenler arasi
% faz farki bozulmadan MUSIC'e tasinir.

    M = size(X,2);
    C = zeros(E,M);
    n = (0:N-1).';

    for ep = 1:E
        idx = (ep-1)*N + (1:N).';
        tauSamp = round(tauChip(ep)*ornPerChip);

        kP = kodPer(mod(n - tauSamp,N)+1);
        carrier = exp(-1j*2*pi*(fIF+fdHz(ep))*(idx-1)/fs);

        Ybb = X(idx,:).*carrier;
        C(ep,:) = sum(Ybb.*conj(kP),1)/N;
    end
end

function [bitsBest,metric,berSteady,inverted] = navBitDogrula(P,refBits,steadyBitStart)
% 1 ms prompt korelator cikislarindan 50 bps GPS navigation bitlerini cikarir.
% Her navigation biti 20 ms oldugu icin 20 prompt epoch koherent toplanir.
% Costas PLL'deki 180 derece faz belirsizligi nedeniyle hem normal hem ters
% polarite denenir ve steady-state bolgesinde daha dusuk BER veren secilir.

    epPerBit = 20;
    nBit = min(numel(refBits), floor(numel(P)/epPerBit));
    P = P(1:nBit*epPerBit);

    % BPSK veri isaretini kaldirarak sabit prompt fazini kestir.
    % P^2 ile 180 derece nav-bit isaret degisimi yok olur.
    phi0 = 0.5*angle(mean(P.^2));
    Prot = P .* exp(-1j*phi0);

    metric = zeros(nBit,1);
    for b = 1:nBit
        ii = (b-1)*epPerBit + (1:epPerBit);
        metric(b) = real(sum(Prot(ii)));
    end

    bits0 = metric < 0;
    ref = logical(refBits(1:nBit));

    k0 = min(max(1,steadyBitStart),nBit);
    ii = k0:nBit;
    ber0 = mean(bits0(ii) ~= ref(ii));
    ber1 = mean(~bits0(ii) ~= ref(ii));

    if ber1 < ber0
        bitsBest = ~bits0;
        berSteady = ber1;
        inverted = true;
    else
        bitsBest = bits0;
        berSteady = ber0;
        inverted = false;
    end
end

function ber = navBERKarsilastir(bits,refBits,steadyBitStart)
% Cozulmus navigation bitlerini verilen referans ile karsilastirir.
% Costas 180 derece polarite belirsizligi icin iki polariteyi de dener.
    nBit = min(numel(bits),numel(refBits));
    k0 = min(max(1,steadyBitStart),nBit);
    ii = k0:nBit;
    b = logical(bits(ii));
    r = logical(refBits(ii));
    ber = min(mean(b ~= r), mean(~b ~= r));
end

function sdb = musicDb(s,nEl)
    sdb = 10*log10(s/max(s(:)));
    if size(sdb,1) ~= nEl
        sdb = sdb.';
    end
end
