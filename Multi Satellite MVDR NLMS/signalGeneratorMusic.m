clear;
clc;
close all;

%% GPS L1 C/A - Aynı PRN Kodunu Kullanan Uydu ve Spoofer - KORELASYONLU MUSIC & MVDR

numbits = 5;            
navdata = randi([0 1], numbits, 1);

fs    = 16.368e6;       
fIF   = 4.092e6;        
fL1   = 1575.42e6;      
fChip = 1.023e6;        
sigma2 = 1;             

% --- 1. Gerçek Uydu (PRN 1) ---
prn1     = 1;
sat1Az   = 30;          % [derece]
sat1El   = 50;          % [derece]
sat1Fd   = 1500;        % Doppler [Hz]
sat1CN0  = 45;          % [dB-Hz]

% --- 2. Gerçek Uydu (PRN 2) ---
prn2     = 2;
sat2Az   = 120;         % [derece]
sat2El   = 70;          % [derece]
sat2Fd   = -500;        % Doppler [Hz]
sat2CN0  = 44;          % [dB-Hz]

% --- SPOOFER: 1. Uydu ile AYNI PRN KODUNU (PRN 1) Kullanır ---
spfPRN       = 1;       % Spoofer'a açıkça Uydu 1'in PRN kodunu atıyoruz
spfAz        = 75;      % PRN 1'e uzaysal olarak yakın [derece]
spfEl        = 60;      % [derece] 
spfFdFark    = 30;      % [Hz] (Gerçek uyduya göre Doppler farkı)
spfKodKayma  = 3;       % [chip] (Spoofer her zaman biraz daha gecikmeli gelir)
spfGucFark   = 1;      % Gerçek uydudan ne kadar güçlü [dB]
spfFaz       = 0;       

%% Temel Bant ve IF Sinyal Üretimi

% Uydu 1 ve Uydu 2 için PRN üretimi
gpswaveobj1 = gpsWaveformGenerator(SignalType="legacy", PRNID=prn1, EnablePCode=false, SampleRate=fs);
baseband1 = gpswaveobj1(navdata);

gpswaveobj2 = gpsWaveformGenerator(SignalType="legacy", PRNID=prn2, EnablePCode=false, SampleRate=fs);
baseband2 = gpswaveobj2(navdata);

% Spoofer için aynı PRN Kodu (spfPRN = 1) ile açık üretim
gpswaveobjSpf = gpsWaveformGenerator(SignalType="legacy", PRNID=spfPRN, EnablePCode=false, SampleRate=fs);
basebandSpf = gpswaveobjSpf(navdata);

t = (0:length(baseband1)-1).' / fs;

% IF Frekansına Taşıma
sat1_if_c = baseband1 .* exp(1j*2*pi*(fIF + sat1Fd)*t);     
sat2_if_c = baseband2 .* exp(1j*2*pi*(fIF + sat2Fd)*t);     

% Spoofer'ın kodu, gerçek uydu 1 ile aynı PRN'den üretilmiş olmasına rağmen
% fiziksel olarak alıcıya 3 chip gecikmeli (kaymış) ulaşır.
kaymaOrnek = round(spfKodKayma * fs / fChip);
spf_bb     = circshift(basebandSpf, kaymaOrnek); % Spoofer baseband'i kaydırılıyor
spf_if_c   = spf_bb .* exp(1j*(2*pi*(fIF + sat1Fd + spfFdFark)*t + spfFaz));

%% 4 Elemanlı CRPA Dizisi Sinyal Karışımı

lambda = physconst("LightSpeed") / fL1;          
dizi = phased.URA(Size=[2 2], ElementSpacing=lambda/2, ArrayNormal="z");

Ps = mean(abs(baseband1).^2);
genlik = @(CN0dB) sqrt(10^(CN0dB/10) * sigma2 / (fs * Ps));

Asat1 = genlik(sat1CN0);
Asat2 = genlik(sat2CN0);
Aspf  = genlik(sat1CN0 + spfGucFark);

Xsat1 = collectPlaneWave(dizi, Asat1*sat1_if_c, [sat1Az; sat1El], fL1);
Xsat2 = collectPlaneWave(dizi, Asat2*sat2_if_c, [sat2Az; sat2El], fL1);
Xspf  = collectPlaneWave(dizi, Aspf*spf_if_c,   [spfAz; spfEl], fL1);

X = Xsat1 + Xsat2 + Xspf + sqrt(sigma2/2) * (randn(size(Xsat1)) + 1j*randn(size(Xsat1)));

%% SİNYAL KORELASYONU (DESPREADING)
N_1ms = round(fs * 1e-3); 
num_epochs = floor(size(X, 1) / N_1ms); 

X_corr_prn1 = zeros(num_epochs, 4);
X_corr_prn2 = zeros(num_epochs, 4);
X_corr_spf  = zeros(num_epochs, 4);

for ep = 1:num_epochs
    idx = (ep-1)*N_1ms + 1 : ep*N_1ms;
    
    % Lokal Kopyalar (Alıcının içeride ürettiği referans sinyaller)
    loc1 = baseband1(idx) .* exp(1j*2*pi*(fIF + sat1Fd)*t(idx));
    loc2 = baseband2(idx) .* exp(1j*2*pi*(fIF + sat2Fd)*t(idx));
    
    % Alıcı, PRN 1 için arama yaparken gecikmeli gelen spoofer'ın kopyasına kilitlendiğini varsayalım
    locSpf = spf_bb(idx) .* exp(1j*(2*pi*(fIF + sat1Fd + spfFdFark)*t(idx) + spfFaz));
    
    X_corr_prn1(ep, :) = sum(X(idx, :) .* conj(loc1), 1);
    X_corr_prn2(ep, :) = sum(X(idx, :) .* conj(loc2), 1);
    X_corr_spf(ep, :)  = sum(X(idx, :) .* conj(locSpf), 1);
end

%% KORELASYON SONRASI BAĞIMSIZ MUSIC KESTİRİMLERİ

azScan = -180:1:180;
elScan = 0:1:90;

musicEst = phased.MUSICEstimator2D(SensorArray=dizi, OperatingFrequency=fL1, ...
    AzimuthScanAngles=azScan, ElevationScanAngles=elScan, ...
    DOAOutputPort=true, NumSignalsSource="Property", NumSignals=1);

[spek1, aci1]     = musicEst(X_corr_prn1);
[spek2, aci2]     = musicEst(X_corr_prn2);
[spekSpf, aciSpf] = musicEst(X_corr_spf);

estAz1 = aci1(1,1); estEl1 = aci1(2,1);
estAz2 = aci2(1,1); estEl2 = aci2(2,1);
estAzSpf = aciSpf(1,1); estElSpf = aciSpf(2,1);

fprintf("\n--- KORELASYONLU MUSIC KESTİRİMLERİ ---\n");
fprintf("1. Uydu (PRN 1) (Gerçek: %d°, %d°) -> Kestirim: %.1f°, %.1f°\n", sat1Az, sat1El, estAz1, estEl1);
fprintf("2. Uydu (PRN 2) (Gerçek: %d°, %d°) -> Kestirim: %.1f°, %.1f°\n", sat2Az, sat2El, estAz2, estEl2);
fprintf("Spoofer (PRN 1) (Gerçek: %d°, %d°) -> Kestirim: %.1f°, %.1f°\n", spfAz, spfEl, estAzSpf, estElSpf);

%% BİRLEŞTİRİLMİŞ MUSIC HEATMAP ÇİZDİRİLMESİ

spektrumToplam = spek1 + spek2 + spekSpf;
spektrumdB = 10*log10(spektrumToplam / max(spektrumToplam(:)));

if size(spektrumdB,1) ~= numel(elScan)
    spektrumdB = spektrumdB.';
end

figure('Position', [150, 150, 900, 500]);
imagesc(azScan, elScan, spektrumdB); 
axis xy; 
c = colorbar; c.Label.String = 'Normalize Güç [dB]';
hold on;

p1 = plot(sat1Az, sat1El, "wo", 'MarkerSize', 12, 'LineWidth', 2);
p2 = plot(sat2Az, sat2El, "co", 'MarkerSize', 12, 'LineWidth', 2);
p3 = plot(spfAz,  spfEl,  "ms", 'MarkerSize', 12, 'LineWidth', 2);
p4 = plot(estAz1, estEl1, "rx", 'MarkerSize', 14, 'LineWidth', 2);
plot(estAz2, estEl2, "rx", 'MarkerSize', 14, 'LineWidth', 2);
plot(estAzSpf, estElSpf, "rx", 'MarkerSize', 14, 'LineWidth', 2);

legend([p1, p2, p3, p4], ...
       {'PRN 1 (Gerçek Uydu)', 'PRN 2 (Gerçek Uydu)', 'PRN 1 (Spoofer)', 'MUSIC Kestirimleri'}, ...
       'TextColor', 'black', 'Color', 'white', 'Location', 'southwest');
       
xlabel("Azimuth [°]");
ylabel("Elevation [°]");
title("Korelasyon Sonrası MUSIC 2D Uzaysal Spektrum Heatmap");

%% ÇOKLU HÜZME İÇİN LCMV (Zorlanmış Null) UYGULAMASI
% MVDR'ın korelasyon zafiyetini aşmak için Spoofer AoA'sını zorla sıfırlıyoruz.

sv = phased.SteeringVector(SensorArray=dizi, PropagationSpeed=physconst("LightSpeed"));

% Kısıtlama (Constraint) vektörlerinin MUSIC kestirimlerinden oluşturulması
cSat1 = sv(fL1, [estAz1; estEl1]);
cSat2 = sv(fL1, [estAz2; estEl2]);
cSpf  = sv(fL1, [estAzSpf; estElSpf]);

% --- Hüzme 1 (PRN 1) ---
% Constraint: [Uydu 1, Spoofer]
% DesiredResponse: [1; 0] -> Yani Uydu 1'e 1 kazanç ver, Spoofer'a zorla 0 kazanç (Null) ver!
lcmv1 = phased.LCMVBeamformer(Constraint=[cSat1, cSpf], DesiredResponse=[1; 0], WeightsOutputPort=true);
[yLCMV1, wLCMV1] = lcmv1(X);

% --- Hüzme 2 (PRN 2) ---
% Constraint: [Uydu 2, Spoofer]
% DesiredResponse: [1; 0] -> Uydu 2'ye 1 kazanç ver, Spoofer'a zorla 0 kazanç (Null) ver!
lcmv2 = phased.LCMVBeamformer(Constraint=[cSat2, cSpf], DesiredResponse=[1; 0], WeightsOutputPort=true);
[yLCMV2, wLCMV2] = lcmv2(X);

%% Anten Diyagramlarının Çizdirilmesi (Zorlanmış Null'lar)

figure('Position', [100, 100, 1000, 400]);

subplot(1,2,1);
pattern(dizi, fL1, -180:1:180, 0:1:90, PropagationSpeed=physconst("LightSpeed"), ...
    Weights=wLCMV1, CoordinateSystem="rectangular", Type="powerdb", Normalize=true);
title("Hüzme 1 (PRN 1 - LCMV)");
xlabel("Azimuth [°]"); ylabel("Elevation [°]");

subplot(1,2,2);
pattern(dizi, fL1, -180:1:180, 0:1:90, PropagationSpeed=physconst("LightSpeed"), ...
    Weights=wLCMV2, CoordinateSystem="rectangular", Type="powerdb", Normalize=true);
title("Hüzme 2 (PRN 2 - LCMV)");
xlabel("Azimuth [°]"); ylabel("Elevation [°]");
% % ÇOKLU HÜZME (MULTI-BEAM) MVDR UYGULAMASI
% 
% Çoklu MVDR algoritmasıyla, her uydu için izole çıkış sinyalleri elde ediliyor[cite: 1, 2]
% MVDR hüzme oluşturucular tüm yönlü gürültü filtrelemesi için ham anten sinyalini (X) kullanır.
% 
% --- Hüzme 1: 1. Uydu (PRN 1) ---
% mvdr1 = phased.MVDRBeamformer(SensorArray=dizi, OperatingFrequency=fL1, ...
%     Direction=[estAz1; estEl1], WeightsOutputPort=true);
% [yMVDR1, wMVDR1] = mvdr1(X_corr_spf);
% 
% --- Hüzme 2: 2. Uydu (PRN 2) ---
% mvdr2 = phased.MVDRBeamformer(SensorArray=dizi, OperatingFrequency=fL1, ...
%     Direction=[estAz2; estEl2], WeightsOutputPort=true);
% [yMVDR2, wMVDR2] = mvdr2(X_corr_spf);
% 
% % Anten Diyagramlarının Çizdirilmesi (İki Bağımsız Hüzme)
% 
% figure('Position', [100, 100, 1000, 400]);
% 
% subplot(1,2,1);
% pattern(dizi, fL1, -180:1:180, 0:1:90, PropagationSpeed=physconst("LightSpeed"), ...
%     Weights=wMVDR1, CoordinateSystem="rectangular", Type="powerdb", Normalize=true);
% title("Hüzme 1 (PRN 1'e Yönelik)");
% xlabel("Azimuth [°]"); ylabel("Elevation [°]");
% 
% subplot(1,2,2);
% pattern(dizi, fL1, -180:1:180, 0:1:90, PropagationSpeed=physconst("LightSpeed"), ...
%     Weights=wMVDR2, CoordinateSystem="rectangular", Type="powerdb", Normalize=true);
% title("Hüzme 2 (PRN 2'ye Yönelik)");
% xlabel("Azimuth [°]"); ylabel("Elevation [°]");