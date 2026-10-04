clear;
clc;
close all;
rng(1);   % Tekrarlanabilirlik (navdata + gürültü her çalıştırmada aynı)

%% GPS L1 C/A - SPOOFER SENARYO ANALİZİ (MUSIC + MVDR + NLMS)

numbits = 5;
navdata = randi([0 1], numbits, 1);

fs     = 16.368e6;
fIF    = 4.092e6;
fL1    = 1575.42e6;
fChip  = 1.023e6;
sigma2 = 1;
c      = physconst("LightSpeed");

% --- 1. Gerçek Uydu (PRN 1) ---
prn1     = 1;
sat1Az   = 30;
sat1El   = 50;
sat1Fd   = 1500;
sat1CN0  = 45;

% --- 2. Gerçek Uydu (PRN 2) ---
prn2     = 2;
sat2Az   = 120;
sat2El   = 70;
sat2Fd   = -500;
sat2CN0  = 44;

% --- SPOOFER: 1. Uydu ile AYNI PRN KODUNU (PRN 1) Kullanır ---
spfPRN       = 1;
spfEl        = 20;
spfFdFark    = 30;
spfKodKayma  = 3;
spfFaz       = 0;

% --- SENARYOLAR: [Azimuth (°)   Güç farkı (dB, PRN 1'e göre)] ---
% Azimuth 30 -> 35 -> 40 -> 45 (5'er derece, 3 artış), sonra 45°'de güç 5 -> 7 -> 10 dB
senaryo = [30  5;
           35  5;
           40  5;
           45  5;
           45  7;
           45 10];
nSen = size(senaryo, 1);

% Sönümün ayrıca raporlanacağı sabit nokta (spoofer'ın son konumu)
takipAz = 45;
takipEl = 50;

%% Temel Bant ve IF Sinyal Üretimi (senaryodan bağımsız, bir kez üretilir)

gpswaveobj1 = gpsWaveformGenerator(SignalType="legacy", PRNID=prn1, EnablePCode=false, SampleRate=fs);
baseband1 = gpswaveobj1(navdata);

gpswaveobj2 = gpsWaveformGenerator(SignalType="legacy", PRNID=prn2, EnablePCode=false, SampleRate=fs);
baseband2 = gpswaveobj2(navdata);

gpswaveobjSpf = gpsWaveformGenerator(SignalType="legacy", PRNID=spfPRN, EnablePCode=false, SampleRate=fs);
basebandSpf = gpswaveobjSpf(navdata);

t = (0:length(baseband1)-1).' / fs;

sat1_if_c = baseband1 .* exp(1j*2*pi*(fIF + sat1Fd)*t);
sat2_if_c = baseband2 .* exp(1j*2*pi*(fIF + sat2Fd)*t);

kaymaOrnek = round(spfKodKayma * fs / fChip);
spf_bb     = circshift(basebandSpf, kaymaOrnek);
spf_if_c   = spf_bb .* exp(1j*(2*pi*(fIF + sat1Fd + spfFdFark)*t + spfFaz));

%% 4 Elemanlı CRPA Dizisi

lambda = c / fL1;
dizi = phased.URA(Size=[2 2], ElementSpacing=lambda/2, ArrayNormal="z");

Ps = mean(abs(baseband1).^2);
genlik = @(CN0dB) sqrt(10^(CN0dB/10) * sigma2 / (fs * Ps));

% Gerçek uydular ve gürültü tüm senaryolarda aynı -> fark sadece spoofer'dan gelir
Xsat1 = collectPlaneWave(dizi, genlik(sat1CN0)*sat1_if_c, [sat1Az; sat1El], fL1);
Xsat2 = collectPlaneWave(dizi, genlik(sat2CN0)*sat2_if_c, [sat2Az; sat2El], fL1);
gurultu = sqrt(sigma2/2) * (randn(size(Xsat1)) + 1j*randn(size(Xsat1)));

N_1ms      = round(fs * 1e-3);
num_epochs = floor(size(Xsat1, 1) / N_1ms);

%% MUSIC / Steering vektör nesneleri

azScan = -180:1:180;
elScan = 0:1:90;

musicEst = phased.MUSICEstimator2D(SensorArray=dizi, OperatingFrequency=fL1, ...
    AzimuthScanAngles=azScan, ElevationScanAngles=elScan, ...
    DOAOutputPort=true, NumSignalsSource="Property", NumSignals=1);

sv = phased.SteeringVector(SensorArray=dizi, PropagationSpeed=c);
aTakip = sv(fL1, [takipAz; takipEl]);

%% Kayıt değişkenleri
spekdB     = cell(nSen, 1);
estSat1    = zeros(2, nSen);
estSat2    = zeros(2, nSen);
estSpf     = zeros(2, nSen);
sonumSpf   = zeros(nSen, 1);   % MVDR-1'in spoofer yönündeki sönümü (hüzme yönüne göre)
sonumTakip = zeros(nSen, 1);   % MVDR-1'in (45°,50°) noktasındaki sönümü
ssrAnten   = zeros(nSen, 1);   % Spoofer/PRN1 korelasyon güç oranı - ham anten (eleman 1)
ssrMVDR    = zeros(nSen, 1);   % Spoofer/PRN1 korelasyon güç oranı - MVDR çıkışı

%% ================= SENARYO DÖNGÜSÜ ================= %%
for k = 1:nSen
    spfAz      = senaryo(k, 1);
    spfGucFark = senaryo(k, 2);

    % --- Sinyal karışımı ---
    Xspf = collectPlaneWave(dizi, genlik(sat1CN0 + spfGucFark)*spf_if_c, [spfAz; spfEl], fL1);
    X = Xsat1 + Xsat2 + Xspf + gurultu;

    % --- Korelasyon (despreading) ---
    Cp1  = korele(X, sat1_if_c, N_1ms, num_epochs);
    Cp2  = korele(X, sat2_if_c, N_1ms, num_epochs);
    Cspf = korele(X, spf_if_c,  N_1ms, num_epochs);

    % --- Korelasyon sonrası MUSIC ---
    [s1, a1] = musicEst(Cp1);
    [s2, a2] = musicEst(Cp2);
    [sS, aS] = musicEst(Cspf);
    estSat1(:, k) = a1(:, 1);
    estSat2(:, k) = a2(:, 1);
    estSpf(:, k)  = aS(:, 1);

    spTop = s1 + s2 + sS;
    spdB  = 10*log10(spTop / max(spTop(:)));
    if size(spdB, 1) ~= numel(elScan), spdB = spdB.'; end
    spekdB{k} = spdB;

    % --- Çoklu hüzme MVDR ---
    mvdr1 = phased.MVDRBeamformer(SensorArray=dizi, OperatingFrequency=fL1, ...
        Direction=a1(:, 1), WeightsOutputPort=true);
    [yMVDR1, wMVDR1] = mvdr1(X);

    mvdr2 = phased.MVDRBeamformer(SensorArray=dizi, OperatingFrequency=fL1, ...
        Direction=a2(:, 1), WeightsOutputPort=true);
    [yMVDR2, wMVDR2] = mvdr2(X);

    % --- MVDR-1 hüzme diyagramından sönüm (hüzme yönüne göre, dB) ---
    gBak   = abs(wMVDR1' * sv(fL1, a1(:, 1)))^2;
    gSpf   = abs(wMVDR1' * sv(fL1, [spfAz; spfEl]))^2;
    gTakip = abs(wMVDR1' * aTakip)^2;
    sonumSpf(k)   = 10*log10(gBak / gSpf);
    sonumTakip(k) = 10*log10(gBak / gTakip);

    % --- Ölçülen spoofer/uydu güç oranı: ham anten vs MVDR çıkışı ---
    Psat_ant = mean(abs(Cp1(:, 1)).^2);
    Pspf_ant = mean(abs(Cspf(:, 1)).^2);
    Psat_mv  = mean(abs(korele(yMVDR1, sat1_if_c, N_1ms, num_epochs)).^2);
    Pspf_mv  = mean(abs(korele(yMVDR1, spf_if_c,  N_1ms, num_epochs)).^2);
    ssrAnten(k) = 10*log10(Pspf_ant / Psat_ant);
    ssrMVDR(k)  = 10*log10(Pspf_mv  / Psat_mv);

    fprintf('Durum %d tamam: Spoofer Az=%d°, El=%d°, dP=%d dB | Sönüm(spf)=%.2f dB\n', ...
        k, spfAz, spfEl, spfGucFark, sonumSpf(k));
end
% Döngü sonunda X, yMVDR1, wMVDR1, wMVDR2 -> SON DURUM (45°, 20°, 10 dB)

%% ZAMANSAL FİLTRELEME: SON DURUM MVDR ÇIKIŞINA NLMS

mvdrSinyalReal = real(yMVDR1);
hamAnten1      = real(X(:, 1));

D = round(1e-3 * fs);
M = 64;

xin   = [zeros(D, 1); mvdrSinyalReal(1:end-D)];
d_sig = mvdrSinyalReal;

nlms = dsp.LMSFilter(M, Method="Normalized LMS");
mumax = maxstep(nlms, xin);
nlms.StepSize = mumax / 5;

[temizSinyal, ~, ~] = nlms(xin, d_sig);

%% ---------------- SÖNÜM TABLOSU ---------------- %%

T = table((1:nSen).', senaryo(:, 1), repmat(spfEl, nSen, 1), senaryo(:, 2), ...
    round(estSpf(1, :).', 1), round(estSpf(2, :).', 1), ...
    round(sonumSpf, 2), round(sonumTakip, 2), ...
    round(ssrAnten, 2), round(ssrMVDR, 2), round(ssrAnten - ssrMVDR, 2), ...
    'VariableNames', {'Durum', 'Spf_Az_deg', 'Spf_El_deg', 'GucFark_dB', ...
                      'MUSIC_SpfAz', 'MUSIC_SpfEl', ...
                      'Sonum_SpfYonu_dB', 'Sonum_45_50_dB', ...
                      'SSR_Anten_dB', 'SSR_MVDR_dB', 'SSR_Iyilesme_dB'});

fprintf('\n================ MVDR (Hüzme 1 - PRN 1) SÖNÜM TABLOSU ================\n');
disp(T);

tabloFig = uifigure('Name', 'MVDR Sönüm Tablosu', 'Position', [150 150 1150 260]);
uitable(tabloFig, 'Data', T, 'Position', [10 10 1130 240]);

%% ---------------- GRAFİKLER ---------------- %%

% 1. GRAFİK: 6 senaryo için MUSIC heatmap
figure('Name', 'MUSIC Heatmap - 6 Senaryo', 'Position', [40, 60, 1500, 800]);
tl = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
for k = 1:nSen
    nexttile;
    imagesc(azScan, elScan, spekdB{k});
    axis xy; colorbar; hold on;
    h1 = plot(sat1Az, sat1El, "wo", 'MarkerSize', 12, 'LineWidth', 2);
    h2 = plot(sat2Az, sat2El, "co", 'MarkerSize', 12, 'LineWidth', 2);
    h3 = plot(senaryo(k, 1), spfEl, "ms", 'MarkerSize', 12, 'LineWidth', 2);
    h4 = plot(estSat1(1, k), estSat1(2, k), "rx", 'MarkerSize', 14, 'LineWidth', 2);
    plot(estSat2(1, k), estSat2(2, k), "rx", 'MarkerSize', 14, 'LineWidth', 2);
    plot(estSpf(1, k),  estSpf(2, k),  "rx", 'MarkerSize', 14, 'LineWidth', 2);
    xlabel("Azimuth [°]"); ylabel("Elevation [°]");
    title(sprintf('Durum %d: Spoofer Az=%d°, El=%d°, \\DeltaP=%d dB', ...
        k, senaryo(k, 1), spfEl, senaryo(k, 2)));
    if k == 1
        legend([h1 h2 h3 h4], {'PRN 1 (Uydu)', 'PRN 2 (Uydu)', 'PRN 1 (Spoofer)', 'MUSIC Kestirimleri'}, ...
            'TextColor', 'black', 'Color', 'white', 'Location', 'southwest');
    end
end
title(tl, 'Korelasyon Sonrası MUSIC 2D Uzaysal Spektrum - Spoofer Senaryoları');

% 2. GRAFİK: Son durum çoklu hüzme MVDR diyagramları
PAT1 = pattern(dizi, fL1, azScan, elScan, PropagationSpeed=c, Weights=wMVDR1, ...
    Type="powerdb", Normalize=true);
PAT2 = pattern(dizi, fL1, azScan, elScan, PropagationSpeed=c, Weights=wMVDR2, ...
    Type="powerdb", Normalize=true);
if size(PAT1, 1) ~= numel(elScan), PAT1 = PAT1.'; PAT2 = PAT2.'; end
PAT1 = max(PAT1, -50);   % Derin null'ları (-Inf) görselleştirme için kırp
PAT2 = max(PAT2, -50);

figure('Name', 'MVDR Beamforming Pattern (Son Durum)', 'Position', [120, 120, 1300, 450]);

subplot(1, 2, 1);
imagesc(azScan, elScan, PAT1); axis xy; colorbar; hold on;
plot(sat1Az, sat1El, "wo", 'MarkerSize', 12, 'LineWidth', 2);
plot(senaryo(end, 1), spfEl, "ms", 'MarkerSize', 12, 'LineWidth', 2);
legend({'PRN 1', 'Spoofer'}, 'TextColor', 'black', 'Color', 'white', 'Location', 'southwest');
xlabel("Azimuth [°]"); ylabel("Elevation [°]");
title(sprintf("Hüzme 1 (PRN 1) | Spoofer (%d°,%d°, %d dB) sönümü: %.2f dB", ...
    senaryo(end, 1), spfEl, senaryo(end, 2), sonumSpf(end)));

subplot(1, 2, 2);
imagesc(azScan, elScan, PAT2); axis xy; colorbar; hold on;
plot(sat2Az, sat2El, "co", 'MarkerSize', 12, 'LineWidth', 2);
plot(senaryo(end, 1), spfEl, "ms", 'MarkerSize', 12, 'LineWidth', 2);
legend({'PRN 2', 'Spoofer'}, 'TextColor', 'black', 'Color', 'white', 'Location', 'southwest');
xlabel("Azimuth [°]"); ylabel("Elevation [°]");
title("Hüzme 2 (PRN 2'ye Yönelik MVDR)");

% 3. GRAFİK: Son durum zaman domeni sinyalleri
N_len = length(mvdrSinyalReal);
startIdx = D + 2000;
gosterimPenceresi = 1500;
endIdx = min(startIdx + gosterimPenceresi, N_len);
plotRange = startIdx:endIdx;

figure('Name', 'Sinyal Karşılaştırması (Son Durum)', 'Position', [100, 400, 1000, 400]);
plot(t(plotRange)*1e6, hamAnten1(plotRange), 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', 'Ham Anten Sinyali'); hold on;
plot(t(plotRange)*1e6, mvdrSinyalReal(plotRange), 'Color', [0 0.447 0.741], 'LineWidth', 1.5, 'DisplayName', 'MVDR Çıkışı');
plot(t(plotRange)*1e6, temizSinyal(plotRange), 'r', 'LineWidth', 2, 'DisplayName', 'NLMS Temizlenmiş Sinyal');
grid on;
legend('Location', 'best');
xlabel('Zaman (\mus)'); ylabel('Genlik');
title(sprintf('Ham / MVDR / NLMS - Son Durum (Spoofer Az=%d°, El=%d°, \\DeltaP=%d dB)', ...
    senaryo(end, 1), spfEl, senaryo(end, 2)));

%% ================================================================
%  KORELASYONLU (KORELASYON SONRASI) MVDR
%
%  Klasik MVDR kovaryansı ham X'ten hesaplar. Ham sinyalde spoofer
%  gürültünün ~17-22 dB altında olduğu için kovaryansta görünmez.
%  Burada kovaryans, spoofer replikasıyla korele edilmiş çıktıdan
%  (Cspf) kestirilir: despreading sonrası spoofer gürültünün ~20-25 dB
%  üstüne çıkar, kovaryansta baskın olur ve MVDR spoofer yönüne null koyar.
%  Hüzme yönü: MUSIC ile kestirilen PRN 1 yönü. Ağırlıklar ham X'e uygulanır.
% ================================================================

Nel = getNumElements(dizi);

korSonumSpf  = zeros(nSen, 1);   % Korelasyonlu MVDR - spoofer yönünde sönüm (hüzme diyagramından)
korPrn1Kayip = zeros(nSen, 1);   % Korelasyonlu MVDR - PRN 1 SNR kazanç kaybı (grafik başlığı için)

for k = 1:nSen
    spfAz      = senaryo(k, 1);
    spfGucFark = senaryo(k, 2);

    % Aynı senaryo sinyali (aynı gürültü ile) yeniden oluşturulur
    Xspf = collectPlaneWave(dizi, genlik(sat1CN0 + spfGucFark)*spf_if_c, [spfAz; spfEl], fL1);
    Xk   = Xsat1 + Xsat2 + Xspf + gurultu;

    Cspf_k = korele(Xk, spf_if_c, N_1ms, num_epochs);   % eğitim verisi

    % --- Korelasyonlu MVDR: kovaryans Cspf'den ---
    mvdrKor = phased.MVDRBeamformer(SensorArray=dizi, OperatingFrequency=fL1, ...
        Direction=estSat1(:, k), WeightsOutputPort=true, TrainingInputPort=true);
    [~, wKor] = mvdrKor(Xk, Cspf_k);

    aD = sv(fL1, estSat1(:, k));        % PRN 1 (hüzme) yönü
    aS = sv(fL1, [spfAz; spfEl]);       % Gerçek spoofer yönü

    % Hüzme diyagramından: hüzme yönüne göre spoofer yönündeki sönüm
    korSonumSpf(k) = 10*log10(abs(wKor'*aD)^2 / abs(wKor'*aS)^2);

    % PRN 1 için SNR kazancı kaybı (ideal 4 elemanlı kazanç = 10log10(Nel) dB'e göre)
    korPrn1Kayip(k) = 10*log10(Nel) - 10*log10(abs(wKor'*aD)^2 / real(wKor'*wKor));

    fprintf('Korelasyonlu MVDR - Durum %d: Sönüm=%.2f dB | PRN1 kaybı=%.2f dB\n', ...
        k, korSonumSpf(k), korPrn1Kayip(k));
end
% Döngü sonunda wKor -> SON DURUM korelasyonlu MVDR ağırlıkları

%% ---------------- KORELASYONLU MVDR TABLOSU ---------------- %%

% Sönüm değerleri hüzme diyagramından hesaplanır:
% 10log10( |w^H a_PRN1|^2 / |w^H a_spoofer|^2 )
T2 = table((1:nSen).', senaryo(:, 1), senaryo(:, 2), ...
    round(sonumSpf, 2), round(korSonumSpf, 2), ...
    'VariableNames', {'Durum', 'Spoofer_Azimut_deg', 'Spoofer_GucFark_dB', ...
                      'Klasik_MVDR_Sonum_dB', 'Korelasyonlu_MVDR_Sonum_dB'});

fprintf('\n=========== KLASİK vs KORELASYONLU MVDR - SPOOFER SÖNÜMÜ ===========\n');
disp(T2);

tabloFig2 = uifigure('Name', 'Korelasyonlu MVDR Tablosu', 'Position', [180 180 820 240]);
uitable(tabloFig2, 'Data', T2, 'Position', [10 10 800 220]);

%% ---------------- KORELASYONLU MVDR GRAFİĞİ ---------------- %%

PATk = pattern(dizi, fL1, azScan, elScan, PropagationSpeed=c, Weights=wKor, ...
    Type="powerdb", Normalize=true);
if size(PATk, 1) ~= numel(elScan), PATk = PATk.'; end
PATk = max(PATk, -50);

figure('Name', 'Klasik vs Korelasyonlu MVDR', 'Position', [60, 80, 1600, 480]);
tl2 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) Klasik MVDR - son durum
nexttile;
imagesc(azScan, elScan, PAT1); axis xy; colorbar; caxis([-50 0]); hold on;
plot(sat1Az, sat1El, "wo", 'MarkerSize', 12, 'LineWidth', 2);
plot(senaryo(end, 1), spfEl, "ms", 'MarkerSize', 12, 'LineWidth', 2);
legend({'PRN 1', 'Spoofer'}, 'TextColor', 'black', 'Color', 'white', 'Location', 'southwest');
xlabel("Azimuth [°]"); ylabel("Elevation [°]");
title(sprintf("Klasik MVDR | Sönüm: %.2f dB", sonumSpf(end)));

% (b) Korelasyonlu MVDR - son durum
nexttile;
imagesc(azScan, elScan, PATk); axis xy; colorbar; caxis([-50 0]); hold on;
plot(sat1Az, sat1El, "wo", 'MarkerSize', 12, 'LineWidth', 2);
plot(senaryo(end, 1), spfEl, "ms", 'MarkerSize', 12, 'LineWidth', 2);
legend({'PRN 1', 'Spoofer'}, 'TextColor', 'black', 'Color', 'white', 'Location', 'southwest');
xlabel("Azimuth [°]"); ylabel("Elevation [°]");
title(sprintf("Korelasyonlu MVDR | Sönüm: %.2f dB | PRN1 kaybı: %.2f dB", ...
    korSonumSpf(end), korPrn1Kayip(end)));

% (c) 6 senaryo için sönüm karşılaştırması
nexttile;
b = bar(1:nSen, [sonumSpf korSonumSpf]);
b(1).FaceColor = [0 0.447 0.741];
b(2).FaceColor = [0.85 0.325 0.098];
grid on;
xticks(1:nSen);
xticklabels(compose('%d°/%ddB', senaryo(:, 1), senaryo(:, 2)));
xtickangle(30);
xlabel('Spoofer Az / Güç Farkı');
ylabel('Spoofer Yönünde Sönüm [dB]');
legend({'Klasik MVDR', 'Korelasyonlu MVDR'}, 'Location', 'northwest');
title('6 Senaryo - Spoofer Sönümü');

title(tl2, sprintf('Klasik vs Korelasyonlu MVDR (Son Durum: Spoofer Az=%d°, El=%d°, \\DeltaP=%d dB)', ...
    senaryo(end, 1), spfEl, senaryo(end, 2)));

%% ================================================================
%  TEMİZ UYDU SİNYALİ vs NLMS ÇIKIŞI (SON DURUM)
%
%  Referans: gürültüsüz, spoofer'sız PRN 1 sinyali. Hüzme 1 PRN 1 yönünde
%  birim kazançlı olduğu için MVDR/NLMS çıkışında ideal olarak görmek
%  istediğimiz sinyal budur.
%  GPS sinyali gürültünün ~20 dB altında olduğundan zaman domeninde
%  benzerlik zayıf görünür; asıl karşılaştırma kod korelasyonu üzerinden
%  yapılır (PRN 1 tepesi 0 chip'te, spoofer tepesi +3 chip'te).
% ================================================================

temizPRN1 = real(genlik(sat1CN0) * sat1_if_c);   % Temiz PRN 1 (gerçel, NLMS ile aynı formatta)
[yKorSon, ~] = mvdrKor(Xk, Cspf_k);               % Son durum korelasyonlu MVDR çıkışı
korMvdrReal  = real(yKorSon);

% --- Korelasyonlu MVDR çıkışına NLMS (klasik MVDR ile aynı ayarlar: D = 1 ms, M = 64, mu = mumax/5) ---
xinKor = [zeros(D, 1); korMvdrReal(1:end-D)];
nlmsKor = dsp.LMSFilter(M, Method="Normalized LMS");
nlmsKor.StepSize = maxstep(nlmsKor, xinKor) / 5;
[temizSinyalKor, ~, ~] = nlmsKor(xinKor, korMvdrReal);

% --- Korelasyonlu MVDR + NLMS zaman domeni grafiği (klasik MVDR grafiğiyle aynı format) ---
figure('Name', 'Sinyal Karşılaştırması - Korelasyonlu MVDR (Son Durum)', 'Position', [140, 440, 1000, 400]);
plot(t(plotRange)*1e6, hamAnten1(plotRange), 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', 'Ham Anten Sinyali'); hold on;
plot(t(plotRange)*1e6, korMvdrReal(plotRange), 'Color', [0.2 0.6 0.2], 'LineWidth', 1.5, 'DisplayName', 'Korelasyonlu MVDR Çıkışı');
plot(t(plotRange)*1e6, temizSinyalKor(plotRange), 'Color', [0.5 0.2 0.7], 'LineWidth', 2, 'DisplayName', 'Korelasyonlu MVDR + NLMS');
grid on;
legend('Location', 'best');
xlabel('Zaman (\mus)'); ylabel('Genlik');
title(sprintf('Ham / Korelasyonlu MVDR / NLMS - Son Durum (Spoofer Az=%d°, El=%d°, \\DeltaP=%d dB)', ...
    senaryo(end, 1), spfEl, senaryo(end, 2)));

% --- (a) Zaman domeni benzerliği: normalize korelasyon katsayısı ---
gecerli  = (D + 1):N_len;                          % NLMS geçici rejimi hariç
rhoFn    = @(a, b) (a.' * b) / (norm(a) * norm(b));
rhoHam   = rhoFn(temizPRN1(gecerli), hamAnten1(gecerli));
rhoMVDR  = rhoFn(temizPRN1(gecerli), mvdrSinyalReal(gecerli));
rhoNLMS  = rhoFn(temizPRN1(gecerli), temizSinyal(gecerli));
rhoKor   = rhoFn(temizPRN1(gecerli), korMvdrReal(gecerli));
rhoKorNLMS = rhoFn(temizPRN1(gecerli), temizSinyalKor(gecerli));
birimRMS = @(s) s / sqrt(mean(s.^2));

% --- (b) Kod korelasyon fonksiyonu (PRN 1 replikası ile) ---
ornPerChip = fs / fChip;
lagChip    = -5:1/ornPerChip:8;
lagOrn     = round(lagChip * ornPerChip);
epochlar   = 3:num_epochs;                          % NLMS oturduktan sonraki epoch'lar

Pkor = [kodKorelasyonu(temizPRN1,      sat1_if_c, N_1ms, epochlar, lagOrn), ...
        kodKorelasyonu(hamAnten1,      sat1_if_c, N_1ms, epochlar, lagOrn), ...
        kodKorelasyonu(mvdrSinyalReal, sat1_if_c, N_1ms, epochlar, lagOrn), ...
        kodKorelasyonu(temizSinyal,    sat1_if_c, N_1ms, epochlar, lagOrn), ...
        kodKorelasyonu(korMvdrReal,    sat1_if_c, N_1ms, epochlar, lagOrn), ...
        kodKorelasyonu(temizSinyalKor, sat1_if_c, N_1ms, epochlar, lagOrn)];
PkorNorm = Pkor ./ max(Pkor, [], 1);                % Her eğri kendi tepesine normalize

i0   = find(lagOrn == 0, 1);                        % PRN 1 kod gecikmesi
iSpf = find(lagOrn == kaymaOrnek, 1);               % Spoofer kod gecikmesi (+3 chip)
tepeOrani_dB = 10*log10(Pkor(iSpf, :) ./ Pkor(i0, :));   % Spoofer tepe / PRN 1 tepe

sinyalAdlari = {'Temiz PRN 1', 'Ham Anten', 'Klasik MVDR', 'Klasik MVDR + NLMS', 'Korelasyonlu MVDR', 'Kor. MVDR + NLMS'};
renkler = [0 0 0; 0.6 0.6 0.6; 0 0.447 0.741; 0.85 0.1 0.1; 0.2 0.6 0.2; 0.5 0.2 0.7];

fprintf('\n=========== TEMİZ PRN 1 İLE KARŞILAŞTIRMA (SON DURUM) ===========\n');
fprintf('Zaman domeni korelasyon katsayısı (rho):\n');
fprintf('  Ham anten         : %.3f\n', rhoHam);
fprintf('  Klasik MVDR       : %.3f\n', rhoMVDR);
fprintf('  Klasik MVDR + NLMS: %.3f\n', rhoNLMS);
fprintf('  Korelasyonlu MVDR : %.3f\n', rhoKor);
fprintf('  Kor. MVDR + NLMS  : %.3f\n', rhoKorNLMS);
fprintf('Kod korelasyonunda spoofer tepesi / PRN 1 tepesi [dB]:\n');
for i = 1:numel(sinyalAdlari)
    fprintf('  %-18s: %6.2f dB\n', sinyalAdlari{i}, tepeOrani_dB(i));
end

% --- FİGÜR ---
figure('Name', 'Temiz Uydu Sinyali vs NLMS', 'Position', [60, 60, 1600, 500]);
tl3 = tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

% (a) Zaman domeni: temiz PRN 1 ve NLMS çıkışı (birim RMS)
nexttile;
plot(t(plotRange)*1e6, birimRMS(temizSinyal(plotRange)), 'Color', [0.95 0.55 0.55], ...
    'LineWidth', 1, 'DisplayName', 'Klasik MVDR + NLMS'); hold on;
plot(t(plotRange)*1e6, birimRMS(temizSinyalKor(plotRange)), 'Color', [0.75 0.6 0.9], ...
    'LineWidth', 1, 'DisplayName', 'Kor. MVDR + NLMS');
plot(t(plotRange)*1e6, birimRMS(temizPRN1(plotRange)), 'k', 'LineWidth', 1.5, ...
    'DisplayName', 'Temiz PRN 1');
grid on; legend('Location', 'best');
xlabel('Zaman (\mus)'); ylabel('Genlik (birim RMS)');
title(sprintf('Zaman Domeni | \\rho: MVDR=%.3f, NLMS=%.3f, KorMVDR=%.3f, KorMVDR+NLMS=%.3f', ...
    rhoMVDR, rhoNLMS, rhoKor, rhoKorNLMS));

% (b) Kod korelasyon fonksiyonları
nexttile;
for i = 1:numel(sinyalAdlari)
    plot(lagChip, PkorNorm(:, i), 'Color', renkler(i, :), 'LineWidth', 1.8, ...
        'DisplayName', sinyalAdlari{i}); hold on;
end
xline(0, '--k', 'PRN 1', 'HandleVisibility', 'off');
xline(spfKodKayma, '--m', 'Spoofer', 'HandleVisibility', 'off');
grid on; legend('Location', 'northeast');
xlabel('Kod Gecikmesi [chip]'); ylabel('Normalize Korelasyon Gücü');
title('PRN 1 Replikası ile Kod Korelasyonu');

% (c) Spoofer tepesi / PRN 1 tepesi
nexttile;
bb = bar(tepeOrani_dB, 'FaceColor', 'flat');
bb.CData = renkler;
yline(0, 'k', 'LineWidth', 1.2);
grid on;
xticks(1:numel(sinyalAdlari)); xticklabels(sinyalAdlari); xtickangle(30);
ylabel('Spoofer Tepe / PRN 1 Tepe [dB]');
title('> 0 dB: alıcı spoofer''a kilitlenir');
text(1:numel(sinyalAdlari), tepeOrani_dB, compose('%.1f', tepeOrani_dB), ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom');

title(tl3, sprintf('Temiz PRN 1 vs MVDR / NLMS Çıkışları (Son Durum: Spoofer Az=%d°, El=%d°, \\DeltaP=%d dB)', ...
    senaryo(end, 1), spfEl, senaryo(end, 2)));

%% ================================================================
%  DLL / PLL TAKİP DÖNGÜSÜ (PRN 1 KANALI)
%
%  Amaç: Alıcının hangi sinyale kilitlendiğini görmek. Her giriş için
%  (temiz PRN 1, ham anten, klasik/korelasyonlu MVDR, NLMS çıkışları):
%    1) Edinme   : Doppler x kod-fazı taraması, en güçlü tepe seçilir.
%                  (Ham antende spoofer 10 dB güçlü -> tepe spoofer'dadır.)
%    2) FLL      : İlk fllSure ms'de frekans çekme (edinme Doppler'i 1 ms
%                  entegrasyonla kaba olduğu için PLL'e girmeden önce).
%    3) PLL      : 2. derece Costas PLL (atan(Q/I)) -> nav bitine duyarsız.
%    4) DLL      : Non-coherent Early-Late güç ayırıcı, 1. derece süzgeç.
%
%  Takip sonunda kod gecikmesi ~0 chip ise PRN 1'e, ~+3 chip ise
%  spoofer'a kilitlenilmiştir.
%
%  Not: Sinyal üretiminde kod Doppler'i yok -> DLL'de taşıyıcı yardımı yok.
%  Tüm girişler (NLMS çıkışları gerçel olsa da) aynı kanaldan geçirilir.
% ================================================================

% --- Replika: PRN 1'in 1 ms'lik dönemi (kod + nav biti işareti, birim RMS) ---
kodPer = baseband1(1:N_1ms);
kodPer = kodPer / sqrt(mean(abs(kodPer).^2));

% --- Döngü parametreleri ---
prm.dEL      = 0.5;     % Early-Late aralığı [chip]
prm.dllBn    = 2;       % DLL gürültü bant genişliği [Hz]
prm.pllBn    = 40;      % PLL gürültü bant genişliği [Hz] (kısa 100 ms'lik kayıtta oturma süresini kısaltır)
prm.fllBn    = 40;      % FLL (çekme aşaması) bant genişliği [Hz] - yalnızca "ozel" yöntem
prm.fllBnTb  = 10;      % FLL bant genişliği [Hz] - gnssSignalTracker çekme aşaması (ilk ms'lerdeki sıçramayı azaltmak için düşük)
prm.fllSure  = 30;      % FLL çekme süresi [ms / epoch]
prm.zeta     = 0.707;   % PLL sönüm oranı

% --- Takip yöntemi ---
%   "toolbox": MATLAB gnssSignalTracker (Satellite Communications Toolbox, R2023b+)
%   "ozel"   : bu dosyadaki takipDLLPLL (FLL çekmeli Costas PLL + E-L DLL)
takipYontemi = "toolbox";

% --- Edinme parametreleri: Doppler penceresi PRN 1'in nominal değeri etrafında ---
fdAra = sat1Fd + (-100:25:100);   % [Hz]
nNC   = 10;                        % Non-coherent toplanan epoch sayısı

% --- Giriş sinyalleri (sinyalAdlari / renkler sırasıyla aynı) ---
girisler = {temizPRN1, X(:, 1), yMVDR1, temizSinyal, yKorSon, temizSinyalKor};
nG = numel(girisler);

takip  = cell(nG, 1);
acqTau = zeros(nG, 1);   % [örnek]
acqFd  = zeros(nG, 1);   % [Hz]  (tabloda gösterilen ilk edinme değeri)

% Yardımlı edinme: NLMS çıkışları, NLMS öncesi sinyalle zaman hizalıdır ve zayıf kaldığı için
% kendi başına edinilemiyor -> kod gecikmesi/Doppler NLMS öncesi sinyalden alınır.
acqKaynak = [1 2 3 3 5 5];

for i = 1:nG
    if acqKaynak(i) == i
        [acqTau(i), acqFd(i)] = kodEdinme(girisler{i}, kodPer, N_1ms, fs, fIF, fdAra, nNC);
    else
        acqTau(i) = acqTau(acqKaynak(i));
        acqFd(i)  = acqFd(acqKaynak(i));
    end

    % 1. geçiş
    takip{i} = takipSec(takipYontemi, girisler{i}, kodPer, N_1ms, num_epochs, fs, fChip, fIF, ...
                        prn1, acqTau(i), acqFd(i), ornPerChip, prm);

    % 2. geçiş (sıcak başlangıç): 1 ms'lik edinmenin Doppler hatası (5-25 Hz) 100 ms'lik kayıtta
    % PLL'in oturmasına yetmiyor. 1. geçişte kilitlenen kanallarda Doppler, 1. geçişin son
    % değeriyle yeniden başlatılır (kilitlenmeyen kanallarda yapılmaz).
    if mean(takip{i}.cos2(end-29:end)) > 0.3
        fdRafine = mean(takip{i}.fd(end-9:end));
        takip{i} = takipSec(takipYontemi, girisler{i}, kodPer, N_1ms, num_epochs, fs, fChip, fIF, ...
                            prn1, acqTau(i), fdRafine, ornPerChip, prm);
    end
end
fprintf('\nTakip yöntemi: %s\n', takipYontemi);

% --- Özet sonuçlar ---
son = num_epochs - 19 : num_epochs;          % son 20 ms
ort = num_epochs - 49 : num_epochs;          % son 50 ms (jitter için)
tauSon = zeros(nG, 1);  fdSon = zeros(nG, 1);
jitterM = zeros(nG, 1); kilitG = zeros(nG, 1);
hedef = strings(nG, 1);
hedefAdlari = ["PRN 1 (gerçek)", "SPOOFER"];
for i = 1:nG
    tauSon(i)  = mean(takip{i}.tauChip(son));
    fdSon(i)   = mean(takip{i}.fd(son));
    jitterM(i) = std(takip{i}.tauChip(ort)) * c / fChip;           % [m]
    kilitG(i)  = mean(takip{i}.cos2(prm.fllSure + 21 : end));      % PLL oturduktan sonra
    [~, j]     = min(abs(tauSon(i) - [0, spfKodKayma]));
    hedef(i)   = hedefAdlari(j);
    if kilitG(i) < 0.3                  % PLL kilitli değil -> kod fazına güvenilmez
        hedef(i) = "KİLİTSİZ / GÜVENİLMEZ";
    end
end

T3 = table(string(sinyalAdlari(:)), ...
    round(mod(acqTau/ornPerChip + 511.5, 1023) - 511.5, 2), acqFd, ...
    round(tauSon, 3), compose("%.1f", fdSon), round(jitterM, 2), round(kilitG, 3), hedef, ...
    'VariableNames', {'Giris', 'Edinme_Kod_chip', 'Edinme_Doppler_Hz', ...
                      'Takip_Kod_chip', 'Takip_Doppler_Hz', 'DLL_Jitter_m', ...
                      'Kilit_Gostergesi', 'Kilitlenen_Hedef'});

fprintf('\n=========== DLL / PLL TAKİP SONUÇLARI (SON DURUM) ===========\n');
fprintf('PRN 1: 0 chip / %d Hz | Spoofer: %d chip / %d Hz\n', sat1Fd, spfKodKayma, sat1Fd + spfFdFark);
disp(T3);

tabloFig3 = uifigure('Name', 'DLL/PLL Takip Tablosu', 'Position', [200 200 1100 260]);
uitable(tabloFig3, 'Data', T3, 'Position', [10 10 1080 240]);

% --- GRAFİK: takip büyüklüklerinin zamansal değişimi ---
tms = (1:num_epochs).';
figure('Name', 'DLL/PLL Takip', 'Position', [60, 60, 1500, 800]);
tl4 = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile;   % (a) DLL: kod gecikmesi
for i = 1:nG
    plot(tms, takip{i}.tauChip, 'Color', renkler(i, :), 'LineWidth', 1.6, ...
        'DisplayName', sinyalAdlari{i}); hold on;
end
yline(0, '--k', 'PRN 1', 'HandleVisibility', 'off');
yline(spfKodKayma, '--m', 'Spoofer', 'HandleVisibility', 'off');
grid on; legend('Location', 'best');
xlabel('Zaman [ms]'); ylabel('Kod gecikmesi [chip]');
title('DLL: Kod Gecikmesi Tahmini');

nexttile;   % (b) PLL/FLL: Doppler
for i = 1:nG
    plot(tms, takip{i}.fd, 'Color', renkler(i, :), 'LineWidth', 1.6, ...
        'DisplayName', sinyalAdlari{i}); hold on;
end
yline(sat1Fd, '--k', 'PRN 1', 'HandleVisibility', 'off');
yline(sat1Fd + spfFdFark, '--m', 'Spoofer', 'HandleVisibility', 'off');
xline(prm.fllSure, ':', 'FLL \rightarrow PLL', 'HandleVisibility', 'off');
grid on;
xlabel('Zaman [ms]'); ylabel('Doppler tahmini [Hz]');
title('FLL/PLL: Doppler Tahmini');

nexttile;   % (c) Kilit göstergesi
for i = 1:nG
    plot(tms, movmean(takip{i}.cos2, 10), 'Color', renkler(i, :), 'LineWidth', 1.6, ...
        'DisplayName', sinyalAdlari{i}); hold on;
end
xline(prm.fllSure, ':', 'HandleVisibility', 'off');
grid on; ylim([-1.05 1.05]);
xlabel('Zaman [ms]'); ylabel('cos(2\phi_e)  (10 ms hareketli ort.)');
title('PLL Kilit Göstergesi (1: kilitli)');

nexttile;   % (d) Prompt genliği
for i = 1:nG
    plot(tms, 20*log10(abs(takip{i}.P)), 'Color', renkler(i, :), 'LineWidth', 1.6, ...
        'DisplayName', sinyalAdlari{i}); hold on;
end
grid on;
xlabel('Zaman [ms]'); ylabel('|P| [dB]');
title('Prompt Korelatör Genliği');

title(tl4, sprintf('DLL / PLL Takip (Son Durum: Spoofer Az=%d°, El=%d°, \\DeltaP=%d dB)', ...
    senaryo(end, 1), spfEl, senaryo(end, 2)));

% --- GRAFİK: Prompt I/Q (Costas PLL kilitliyse noktalar I ekseni üzerinde) ---
figure('Name', 'Prompt I/Q', 'Position', [100, 100, 1400, 750]);
tl5 = tiledlayout(2, 3, 'TileSpacing', 'loose', 'Padding', 'loose');
for i = 1:nG
    nexttile;
    kk = prm.fllSure + 1 : num_epochs;
    scatter(real(takip{i}.P(kk)), imag(takip{i}.P(kk)), 18, renkler(i, :), 'filled'); hold on;
    xline(0, 'Color', [0.7 0.7 0.7]); yline(0, 'Color', [0.7 0.7 0.7]);
    axis equal; grid on;
    xlabel('I (Prompt)'); ylabel('Q (Prompt)');
    title(sprintf('%s | %s', sinyalAdlari{i}, hedef(i)));
end
title(tl5, 'Prompt Korelatör I/Q (PLL sonrası) - Veri bitleri ± I ekseninde');

%% ---------------- YEREL FONKSİYON ---------------- %%
function C = korele(Y, rep, N, E)
% 1 ms'lik bloklarla korelasyon: Y (L x K) sinyali, rep (L x 1) replika
% Çıkış: C (E x K), her satır bir epoch, her sütun bir anten elemanı
    L = N * E;
    Z = Y(1:L, :) .* conj(rep(1:L));
    C = squeeze(sum(reshape(Z, N, E, []), 1));
    if isvector(C), C = C(:); end
end

function P = kodKorelasyonu(s, rep, N, epochlar, lagOrn)
% Kod gecikmesine göre korelasyon gücü:
% her epoch'ta 1 ms koherent dairesel korelasyon (FFT ile),
% epoch'lar üzerinde koherent olmayan ortalama.
    P = zeros(numel(lagOrn), 1);
    idxLag = mod(lagOrn(:), N) + 1;
    for ep = epochlar
        idx = (ep-1)*N + (1:N);
        cc  = ifft(fft(s(idx)) .* conj(fft(rep(idx))));
        P   = P + abs(cc(idxLag)).^2;
    end
    P = P / numel(epochlar);
end

function [tauAcq, fdAcq] = kodEdinme(y, kodPer, N, fs, fIF, fdAra, nNC)
% Basit edinme: her Doppler hücresinde 1 ms'lik koherent dairesel korelasyon
% (FFT ile), nNC epoch üzerinde non-coherent toplama. Çıkış: en güçlü hücre.
% tauAcq: kod gecikmesi [örnek, 0..N-1], fdAcq: Doppler [Hz]
    Kf = conj(fft(kodPer));
    S  = zeros(numel(fdAra), N);
    for i = 1:numel(fdAra)
        for ep = 1:nNC
            idx = (ep-1)*N + (1:N).';
            yb  = y(idx) .* exp(-1j*2*pi*(fIF + fdAra(i)) * (idx-1)/fs);
            cc  = ifft(fft(yb) .* Kf);
            S(i, :) = S(i, :) + abs(cc).' .^ 2;
        end
    end
    [~, lin] = max(S(:));
    [iD, iT] = ind2sub(size(S), lin);
    fdAcq  = fdAra(iD);
    tauAcq = iT - 1;
end

function s = takipDLLPLL(y, kodPer, N, E, fs, fChip, fIF, fdBas, tauBas, prm)
% FLL-çekmeli 2. derece Costas PLL + non-coherent Early-Late DLL.
%  Zaman  : her epoch = 1 ms (T = N/fs)
%  Replika: kodPer'in dairesel kaydırılmış hali (örnek çözünürlüğü, 1/16 chip)
%  Kod modeli: y(n) ~ kod(n - tau)  ->  Early = tau - h, Late = tau + h gecikmeli
    T       = N / fs;
    ornChip = fs / fChip;                         % örnek/chip (=16)
    h       = round(prm.dEL/2 * ornChip);         % E/L yarı aralığı [örnek]
    n       = (0:N-1).';

    wn   = prm.pllBn / (0.5*(prm.zeta + 1/(4*prm.zeta)));   % PLL doğal frekansı [rad/s]
    kDLL = 4 * prm.dllBn * T;                                % 1. derece DLL kazancı
    kFLL = 4 * prm.fllBn * T;                                % 1. derece FLL kazancı

    fdEst  = fdBas;        % Doppler tahmini [Hz]
    tauEst = tauBas;       % kod gecikmesi tahmini [örnek]
    fazCyc = 0;            % taşıyıcı NCO fazı [cycle]
    wInt   = 0;            % PLL integratörü [rad/s]
    Pold   = 0;

    s.tauChip = zeros(E, 1);
    s.fd      = zeros(E, 1);
    s.P       = zeros(E, 1);
    s.cos2    = zeros(E, 1);
    s.phDeg   = zeros(E, 1);

    for ep = 1:E
        idx  = (ep-1)*N + (1:N).';
        fEst = fIF + fdEst;

        % --- Taşıyıcı sıyırma + E/P/L korelatörleri ---
        yb   = y(idx) .* exp(-1j*2*pi*(fazCyc + fEst*n/fs));
        tauO = round(tauEst);
        kP   = kodPer(mod(n - tauO,     N) + 1);
        kE   = kodPer(mod(n - tauO + h, N) + 1);
        kL   = kodPer(mod(n - tauO - h, N) + 1);
        P    = sum(yb .* conj(kP));
        Ee   = sum(yb .* conj(kE));
        Ll   = sum(yb .* conj(kL));

        % --- DLL: normalize non-coherent E-L ayırıcı (üçgen korelasyon için) ---
        dllErr = (abs(Ee) - abs(Ll)) / (abs(Ee) + abs(Ll));
        tauEst = tauEst - kDLL * (1 - prm.dEL/2) * dllErr * ornChip;

        % --- PLL ayırıcı (Costas): nav bitine duyarsız ---
        phiErr = atan(imag(P) / real(P));          % [rad]

        if ep <= prm.fllSure
            % FLL çekme: ardışık prompt'lar arası çapraz/nokta çarpım ayırıcısı
            if ep > 1
                dotP  = real(P * conj(Pold));
                crsP  = imag(conj(Pold) * P);
                fErr  = atan(crsP / dotP) / (2*pi*T);       % [Hz]
                fdEst = fdEst + kFLL * fErr;
            end
            wInt = 2*pi*(fdEst - fdBas);           % PLL'e geçişte süreklilik
        else
            % 2. derece PLL döngü süzgeci
            wInt  = wInt + wn^2 * T * phiErr;
            fdEst = fdBas + (wInt + 2*prm.zeta*wn*phiErr) / (2*pi);
        end
        Pold   = P;
        fazCyc = mod(fazCyc + fEst*N/fs, 1);

        s.tauChip(ep) = mod(tauEst/ornChip + (N/ornChip)/2, N/ornChip) - (N/ornChip)/2;  % (-511.5, 511.5]
        s.fd(ep)      = fdEst;
        s.P(ep)       = P;
        s.cos2(ep)    = cos(2*phiErr);
        s.phDeg(ep)   = phiErr * 180/pi;
    end
end

function s = takipToolbox(y, fs, fIF, prnID, N, E, acqChip, acqFd, prm)
% gnssSignalTracker ile takip (PLL + FLL + DLL birlikte çalışır).
%  - Nesne, gpsWaveformGenerator'ın temel bant formatını bekler (dokümandaki örnek gibi);
%    bu yüzden giriş analitik hale getirilip temel banda indirilir (IF = 0).
%  - InitialCodePhaseOffset C/A için TAMSAYI chip ister; edinme çıktısı yuvarlanır.
%  - FLL bant genişliği ilk fllSure ms'de geniş (çekme), sonra dar (tunable özellik).
%  Not: trInfo.DelayEstimate / FrequencyEstimate birimleri için aşağıdaki varsayımlar
%  geçerlidir; temiz PRN 1 satırıyla (gerçek: 0 chip, sat1Fd Hz) doğrulayın.
    if isreal(y), y = hilbert(y); end                 % gerçel IF -> analitik sinyal
    t   = (0:E*N-1).' / fs;
    yBB = y(1:E*N) .* exp(-1j*2*pi*fIF*t);            % temel banda indir (Doppler kalır)

    chip0 = mod(round(acqChip), 1023);                % [0, 1022]
    gst = gnssSignalTracker(GNSSSignalType="GPS C/A", SampleRate=fs, PRNID=prnID, ...
        InitialCodePhaseOffset=chip0, InitialFrequencyOffset=acqFd, ...
        PLLNoiseBandwidth=prm.pllBn, FLLNoiseBandwidth=prm.fllBnTb, ...
        DLLNoiseBandwidth=prm.dllBn);

    s.tauChip = zeros(E, 1);
    s.fd      = zeros(E, 1);
    s.P       = zeros(E, 1);
    s.cos2    = zeros(E, 1);
    s.phDeg   = zeros(E, 1);

    for ep = 1:E
        if ep == prm.fllSure + 1
            gst.FLLNoiseBandwidth = 4;                % çekme bitti -> FLL'i daralt
        end
        [P, info] = gst(yBB((ep-1)*N + (1:N)));
        P = -1j * P;                                  % C/A kuadratür (Q) eksenindedir -> I eksenine döndür (-90°)
        phiErr = atan(imag(P) / real(P));             % Costas: nav bitine duyarsız

        s.P(ep)       = P;
        s.cos2(ep)    = cos(2*phiErr);
        s.phDeg(ep)   = phiErr * 180/pi;
        % FrequencyEstimate, edinme Doppler'ine göre GÖRELİ ve işareti ters çıkıyor (çıktılarda 4 satırda
        % acqFd - tahmin değeri gerçek Doppler'e 3-6 Hz içinde uyuyor). Doğrulama: Doppler grafiğinde
        % eğriler kesikli PRN 1 / Spoofer çizgilerine oturmalı.
        s.fd(ep)      = acqFd - info.FrequencyEstimate;
        s.tauChip(ep) = mod(chip0 + info.DelayEstimate + 511.5, 1023) - 511.5;  % varsayım: chip, edinme ofsetine göre bağıl
    end
end

function s = takipSec(yontem, y, kodPer, N, E, fs, fChip, fIF, prnID, tauOrn, fdBas, ornChip, prm)
% Seçilen yönteme göre takip kanalını çalıştırır (tauOrn: edinme kod gecikmesi [örnek]).
    switch yontem
        case "toolbox"
            s = takipToolbox(y, fs, fIF, prnID, N, E, tauOrn/ornChip, fdBas, prm);
        otherwise
            s = takipDLLPLL(y, kodPer, N, E, fs, fChip, fIF, fdBas, tauOrn, prm);
    end
end
