clear;
clc;
close all;

%% GPS L1 C/A - 1 Uydu + 1 Spoofer (Aynı PRN) - KORELASYONLU MUSIC & STANDART MVDR
% Bu kodda LCMV kısıtlaması yoktur. MVDR sadece uydu yönünü bilir ve 
% yüksek enerjiyi kendi kendine (otomatik) ezmeye çalışır[cite: 1, 2].

numbits = 5;            
navdata = randi([0 1], numbits, 1);

fs    = 16.368e6;       
fIF   = 4.092e6;        
fL1   = 1575.42e6;      
fChip = 1.023e6;        
sigma2 = 1;             

% --- 1. Gerçek Uydu (PRN 1) ---
prn1     = 1;
satAz    = 30;          % [derece]
satEl    = 50;          % [derece]
satFd    = 1500;        % Doppler [Hz]
satCN0   = 45;          % [dB-Hz]

% --- SPOOFER: Gerçek uydu ile AYNI PRN KODUNU (PRN 1) Kullanır ---
spfAz        = -10;     % [derece] 
spfEl        = 55;      % [derece] 
spfFdFark    = 30;      % [Hz]
spfKodKayma  = 3;      % Sinyal gecikmesi [chip]
spfGucFark   = 10;      % Gerçek uydudan ne kadar güçlü [dB]
spfFaz       = 0;       

%% Temel Bant ve IF Sinyal Üretimi

gpswaveobj = gpsWaveformGenerator(SignalType="legacy", PRNID=prn1, EnablePCode=false, SampleRate=fs);
baseband = gpswaveobj(navdata);

t = (0:length(baseband)-1).' / fs;

sat_if_c = baseband .* exp(1j*2*pi*(fIF + satFd)*t);     

% Spoofer kopyası
kaymaOrnek = round(spfKodKayma * fs / fChip);
spf_bb     = circshift(baseband, kaymaOrnek);
spf_if_c   = spf_bb .* exp(1j*(2*pi*(fIF + satFd + spfFdFark)*t + spfFaz));

%% 4 Elemanlı CRPA Dizisi Sinyal Karışımı

lambda = physconst("LightSpeed") / fL1;          
dizi = phased.URA(Size=[2 2], ElementSpacing=lambda/2, ArrayNormal="z");

Ps = mean(abs(baseband).^2);
genlik = @(CN0dB) sqrt(10^(CN0dB/10) * sigma2 / (fs * Ps));

Asat = genlik(satCN0);
Aspf = genlik(satCN0 + spfGucFark);

Xsat = collectPlaneWave(dizi, Asat*sat_if_c, [satAz; satEl], fL1);
Xspf = collectPlaneWave(dizi, Aspf*spf_if_c, [spfAz; spfEl], fL1);

X = Xsat + Xspf + sqrt(sigma2/2) * (randn(size(Xsat)) + 1j*randn(size(Xsat)));

%% SİNYAL KORELASYONU (DESPREADING)
N_1ms = round(fs * 1e-3); 
num_epochs = floor(size(X, 1) / N_1ms); 

X_corr_sat = zeros(num_epochs, 4);
X_corr_spf = zeros(num_epochs, 4);

for ep = 1:num_epochs
    idx = (ep-1)*N_1ms + 1 : ep*N_1ms;
    
    locSat = baseband(idx) .* exp(1j*2*pi*(fIF + satFd)*t(idx));
    locSpf = spf_bb(idx) .* exp(1j*(2*pi*(fIF + satFd + spfFdFark)*t(idx) + spfFaz));
    
    X_corr_sat(ep, :) = sum(X(idx, :) .* conj(locSat), 1);
    X_corr_spf(ep, :) = sum(X(idx, :) .* conj(locSpf), 1);
end

%% KORELASYON SONRASI BAĞIMSIZ MUSIC KESTİRİMLERİ

azScan = -180:1:180;
elScan = 0:1:90;

musicEst = phased.MUSICEstimator2D(SensorArray=dizi, OperatingFrequency=fL1, ...
    AzimuthScanAngles=azScan, ElevationScanAngles=elScan, ...
    DOAOutputPort=true, NumSignalsSource="Property", NumSignals=1);

[spekSat, aciSat] = musicEst(X_corr_sat);
[spekSpf, aciSpf] = musicEst(X_corr_spf);

estAzSat = aciSat(1,1); estElSat = aciSat(2,1);
estAzSpf = aciSpf(1,1); estElSpf = aciSpf(2,1);

fprintf("\n--- KESTİRİMLER ---\n");
fprintf("Gerçek Uydu -> Kestirim: %.1f°, %.1f°\n", estAzSat, estElSat);
fprintf("Spoofer     -> Kestirim: %.1f°, %.1f°\n", estAzSpf, estElSpf);

%% BİRLEŞTİRİLMİŞ MUSIC HEATMAP

spektrumToplam = spekSat + spekSpf;
spektrumdB = 10*log10(spektrumToplam / max(spektrumToplam(:)));

if size(spektrumdB,1) ~= numel(elScan)
    spektrumdB = spektrumdB.';
end

figure('Position', [150, 150, 800, 500]);
imagesc(azScan, elScan, spektrumdB); 
axis xy; colorbar; hold on;
plot(satAz, satEl, "wo", 'MarkerSize', 12, 'LineWidth', 2);
plot(spfAz, spfEl, "ms", 'MarkerSize', 12, 'LineWidth', 2);
plot(estAzSat, estElSat, "rx", 'MarkerSize', 14, 'LineWidth', 2);
plot(estAzSpf, estElSpf, "rx", 'MarkerSize', 14, 'LineWidth', 2);
legend('Gerçek Uydu', 'Spoofer', 'MUSIC Kestirimleri', 'TextColor', 'black', 'Color', 'white', 'Location', 'southwest');
xlabel("Azimuth [°]"); ylabel("Elevation [°]");
title("Korelasyon Sonrası MUSIC 2D Uzaysal Spektrum Heatmap");

%% STANDART MVDR UYGULAMASI (LCMV YOK)

% Algoritma sadece uydunun yönünü biliyor, Spoofer'ı genel kovaryans matrisi (X) 
% üzerinden otomatik ezmeye çalışacak[cite: 1, 2].
mvdr = phased.MVDRBeamformer( ...
    SensorArray=dizi, ...
    OperatingFrequency=fL1, ...
    Direction=[estAzSat; estElSat], ...
    WeightsOutputPort=true);

[yMVDR, wMVDR] = mvdr(X);

%% MVDR Anten Diyagramının Çizdirilmesi

figure('Position', [300, 200, 600, 500]);
pattern(dizi, fL1, -180:1:180, 0:1:90, PropagationSpeed=physconst("LightSpeed"), ...
    Weights=wMVDR, CoordinateSystem="rectangular", Type="powerdb", Normalize=true);
title("Standart MVDR Anten Diyagramı (LCMV Yok)");
xlabel("Azimuth [°]"); ylabel("Elevation [°]"); zlabel("Normalize Güç [dB]");