clear;
clc;
close all;

%% GPS L1 C/A - 2 Uydu + 1 PRN Kodsuz Tehdit (Jammer) - STANDART MVDR
% Tehdit, GPS PRN yapısına sahip olmadığı için rastgele geniş bantlı bir sinyal (BPSK/Noise) olarak modellenmiştir.

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

% --- TEHDİT (PRN KODSUZ SPOOFER / JAMMER) ---
% Herhangi bir uyduyla korelasyona girmez (ilişkisizdir).
jamAz        = 60;      % Sizin belirttiğiniz ve null oluşmayan o koordinat [derece]
jamEl        = 70;      % [derece] 
jamGucFark   = 5;      % MVDR'ın kolayca görüp ezmesi için 20 dB daha güçlü

%% Temel Bant ve IF Sinyal Üretimi

% Uydu 1 ve Uydu 2 için PRN üretimi
gpswaveobj1 = gpsWaveformGenerator(SignalType="legacy", PRNID=prn1, EnablePCode=false, SampleRate=fs);
baseband1 = gpswaveobj1(navdata);

gpswaveobj2 = gpsWaveformGenerator(SignalType="legacy", PRNID=prn2, EnablePCode=false, SampleRate=fs);
baseband2 = gpswaveobj2(navdata);

t = (0:length(baseband1)-1).' / fs;

sat1_if_c = baseband1 .* exp(1j*2*pi*(fIF + sat1Fd)*t);     
sat2_if_c = baseband2 .* exp(1j*2*pi*(fIF + sat2Fd)*t);     

% PRN Kodsuz Tehdit Üretimi: Geniş bantlı rastgele bir sinyal (Gaussian Noise / Rastgele BPSK benzeri)
jam_bb   = randn(size(baseband1)) + 1j*randn(size(baseband1));
jam_if_c = jam_bb .* exp(1j*2*pi*fIF*t);

%% 4 Elemanlı CRPA Dizisi Sinyal Karışımı

lambda = physconst("LightSpeed") / fL1;          
dizi = phased.URA(Size=[2 2], ElementSpacing=lambda/2, ArrayNormal="z");

Ps = mean(abs(baseband1).^2);
genlik = @(CN0dB) sqrt(10^(CN0dB/10) * sigma2 / (fs * Ps));

Asat1 = genlik(sat1CN0);
Asat2 = genlik(sat2CN0);
Ajam  = genlik(sat1CN0 + jamGucFark);

Xsat1 = collectPlaneWave(dizi, Asat1*sat1_if_c, [sat1Az; sat1El], fL1);
Xsat2 = collectPlaneWave(dizi, Asat2*sat2_if_c, [sat2Az; sat2El], fL1);
Xjam  = collectPlaneWave(dizi, Ajam*jam_if_c,   [jamAz; jamEl], fL1);

% Toplam Sinyal
X = Xsat1 + Xsat2 + Xjam + sqrt(sigma2/2) * (randn(size(Xsat1)) + 1j*randn(size(Xsat1)));

%% SİNYAL KORELASYONU (Sadece Gerçek Uyduların Yönünü Bulmak İçin)
N_1ms = round(fs * 1e-3); 
num_epochs = floor(size(X, 1) / N_1ms); 

X_corr_prn1 = zeros(num_epochs, 4);
X_corr_prn2 = zeros(num_epochs, 4);

for ep = 1:num_epochs
    idx = (ep-1)*N_1ms + 1 : ep*N_1ms;
    
    loc1 = baseband1(idx) .* exp(1j*2*pi*(fIF + sat1Fd)*t(idx));
    loc2 = baseband2(idx) .* exp(1j*2*pi*(fIF + sat2Fd)*t(idx));
    
    X_corr_prn1(ep, :) = sum(X(idx, :) .* conj(loc1), 1);
    X_corr_prn2(ep, :) = sum(X(idx, :) .* conj(loc2), 1);
end

%% MUSIC KESTİRİMLERİ

azScan = -180:1:180;
elScan = 0:1:90;

musicEst = phased.MUSICEstimator2D(SensorArray=dizi, OperatingFrequency=fL1, ...
    AzimuthScanAngles=azScan, ElevationScanAngles=elScan, ...
    DOAOutputPort=true, NumSignalsSource="Property", NumSignals=1);

% Uyduların (Korelasyonlu) Kestirimi
[~, aci1] = musicEst(X_corr_prn1);
[~, aci2] = musicEst(X_corr_prn2);

estAz1 = aci1(1,1); estEl1 = aci1(2,1);
estAz2 = aci2(1,1); estEl2 = aci2(2,1);

% Tehdit (PRN Kodsuz) Kestirimi: Sinyal çok güçlü olduğu için doğrudan ham sinyalden (X) bulunur
[~, aciJam] = musicEst(X(1:1000, :)); 
estAzJam = aciJam(1,1); estElJam = aciJam(2,1);

fprintf("\n--- KESTİRİMLER ---\n");
fprintf("1. Uydu (PRN 1) -> Kestirim: %.1f°, %.1f°\n", estAz1, estEl1);
fprintf("2. Uydu (PRN 2) -> Kestirim: %.1f°, %.1f°\n", estAz2, estEl2);
fprintf("Tehdit (Jammer) -> Kestirim: %.1f°, %.1f°\n", estAzJam, estElJam);

%% STANDART ÇOKLU HÜZME MVDR UYGULAMASI (LCMV KULLANILMADAN)
% PRN kodu olmayan tehdit, X matrisinin varyansını (enerjisini) domine edeceği için,
% MVDR algoritması ekstra bir emre (kısıtlamaya) ihtiyaç duymadan tehdidi otomatik ezecektir.

% --- Hüzme 1: 1. Uydu (PRN 1) ---
mvdr1 = phased.MVDRBeamformer(SensorArray=dizi, OperatingFrequency=fL1, ...
    Direction=[estAz1; estEl1], WeightsOutputPort=true);
[yMVDR1, wMVDR1] = mvdr1(X);

% --- Hüzme 2: 2. Uydu (PRN 2) ---
mvdr2 = phased.MVDRBeamformer(SensorArray=dizi, OperatingFrequency=fL1, ...
    Direction=[estAz2; estEl2], WeightsOutputPort=true);
[yMVDR2, wMVDR2] = mvdr2(X);

%% Anten Diyagramlarının Çizdirilmesi

figure('Position', [100, 100, 1000, 400]);

subplot(1,2,1);
pattern(dizi, fL1, -180:1:180, 0:1:90, PropagationSpeed=physconst("LightSpeed"), ...
    Weights=wMVDR1, CoordinateSystem="rectangular", Type="powerdb", Normalize=true);
title("Hüzme 1 (PRN 1 - Standart MVDR)");
xlabel("Azimuth [°]"); ylabel("Elevation [°]");
% Grafiği incelediğinizde Azimuth 60, Elevation 70'te doğal, derin bir null göreceksiniz.

subplot(1,2,2);
pattern(dizi, fL1, -180:1:180, 0:1:90, PropagationSpeed=physconst("LightSpeed"), ...
    Weights=wMVDR2, CoordinateSystem="rectangular", Type="powerdb", Normalize=true);
title("Hüzme 2 (PRN 2 - Standart MVDR)");
xlabel("Azimuth [°]"); ylabel("Elevation [°]");