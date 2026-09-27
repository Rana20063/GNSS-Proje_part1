clear;
clc;
close all;

%% GPS L1 C/A - Aynı PRN Kodunu Kullanan Uydu ve Spoofer 
% KORELASYONLU MUSIC + LCMV (Zorlanmış Null) + NLMS (Zamansal Filtre)

numbits = 5;            
navdata = randi([0 1], numbits, 1);

fs    = 16.368e6;       
fIF   = 4.092e6;        
fL1   = 1575.42e6;      
fChip = 1.023e6;        
sigma2 = 1;             

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
spfAz        = -120;    
spfEl        = 60;      
spfFdFark    = 30;      
spfKodKayma  = 3;       
spfGucFark   = 5;       
spfFaz       = 0;       

%% Temel Bant ve IF Sinyal Üretimi

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

%% 4 Elemanlı CRPA Dizisi Sinyal Karışımı

lambda = physconst("LightSpeed") / fL1;          
dizi = phased.URA(Size=[2 2], ElementSpacing=lambda/2, ArrayNormal="z");

Ps = mean(abs(baseband1).^2);
genlik = @(CN0dB) sqrt(10^(CN0dB/10) * sigma2 / (fs * Ps));

Xsat1 = collectPlaneWave(dizi, genlik(sat1CN0)*sat1_if_c, [sat1Az; sat1El], fL1);
Xsat2 = collectPlaneWave(dizi, genlik(sat2CN0)*sat2_if_c, [sat2Az; sat2El], fL1);
Xspf  = collectPlaneWave(dizi, genlik(sat1CN0 + spfGucFark)*spf_if_c, [spfAz; spfEl], fL1);

X = Xsat1 + Xsat2 + Xspf + sqrt(sigma2/2) * (randn(size(Xsat1)) + 1j*randn(size(Xsat1)));

%% SİNYAL KORELASYONU (DESPREADING)
N_1ms = round(fs * 1e-3); 
num_epochs = floor(size(X, 1) / N_1ms); 

X_corr_prn1 = zeros(num_epochs, 4);
X_corr_prn2 = zeros(num_epochs, 4);
X_corr_spf  = zeros(num_epochs, 4);

for ep = 1:num_epochs
    idx = (ep-1)*N_1ms + 1 : ep*N_1ms;
    
    loc1 = baseband1(idx) .* exp(1j*2*pi*(fIF + sat1Fd)*t(idx));
    loc2 = baseband2(idx) .* exp(1j*2*pi*(fIF + sat2Fd)*t(idx));
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

%% ÇOKLU HÜZME (MULTI-BEAM) LCMV UYGULAMASI (UZAYSAL FİLTRE)
% MUSIC istihbaratı ile Spoofer yönüne zorla Null atılır[cite: 1, 2].

sv = phased.SteeringVector(SensorArray=dizi, PropagationSpeed=physconst("LightSpeed"));

cSat1 = sv(fL1, [estAz1; estEl1]);
cSat2 = sv(fL1, [estAz2; estEl2]);
cSpf  = sv(fL1, [estAzSpf; estElSpf]);

% --- Hüzme 1: 1. Uydu (PRN 1) ---
% Constraint: Uydu 1 ve Spoofer; DesiredResponse: Uydu 1'e 1, Spoofer'a 0 kazanç
lcmv1 = phased.LCMVBeamformer(Constraint=[cSat1, cSpf], DesiredResponse=[1; 0], WeightsOutputPort=true);
[yLCMV1, wLCMV1] = lcmv1(X);

% --- Hüzme 2: 2. Uydu (PRN 2) ---
% Constraint: Uydu 2 ve Spoofer; DesiredResponse: Uydu 2'ye 1, Spoofer'a 0 kazanç
lcmv2 = phased.LCMVBeamformer(Constraint=[cSat2, cSpf], DesiredResponse=[1; 0], WeightsOutputPort=true);
[yLCMV2, wLCMV2] = lcmv2(X);

%% ZAMANSAL FİLTRELEME: LCMV ÇIKIŞINA NLMS UYGULAMASI

lcmvSinyalReal = real(yLCMV1);    
sat1_if_real   = real(sat1_if_c); 
hamAnten1      = real(X(:,1));    

D = round(1e-3 * fs);       
M = 64;                     

xin = [zeros(D,1); lcmvSinyalReal(1:end-D)];  
d_sig = lcmvSinyalReal;                                 

nlms = dsp.LMSFilter(M, Method="Normalized LMS");
mumax = maxstep(nlms, xin);
nlms.StepSize = mumax / 5; 

[temizSinyal, ~, ~] = nlms(xin, d_sig);   

%% ---------------- GRAFİKLER ---------------- %%

N_len = length(lcmvSinyalReal);
startIdx = D + 2000; 
gosterimPenceresi = 1500;
endIdx = min(startIdx + gosterimPenceresi, N_len);
plotRange = startIdx:endIdx;

% 1. GRAFİK: Zaman Domeni Sinyalleri (Üst Üste)
figure('Name', 'Sinyal Karşılaştırması', 'Position', [100, 400, 1000, 400]);
plot(t(plotRange)*1e6, hamAnten1(plotRange), 'Color', [0.7 0.7 0.7], 'LineWidth', 1, 'DisplayName', 'Ham Anten Sinyali'); hold on;
plot(t(plotRange)*1e6, lcmvSinyalReal(plotRange), 'Color', [0 0.447 0.741], 'LineWidth', 1.5, 'DisplayName', 'LCMV Çıkışı');
plot(t(plotRange)*1e6, temizSinyal(plotRange), 'r', 'LineWidth', 2, 'DisplayName', 'NLMS Temizlenmiş Sinyal');
grid on;
legend('Location', 'best');
xlabel('Zaman (\mus)'); ylabel('Genlik');
title('Ham Sinyal, LCMV ve NLMS Çıktılarının Üst Üste Karşılaştırması');

% 2. GRAFİK: MUSIC Heatmap
spektrumToplam = spek1 + spek2 + spekSpf;
spektrumdB = 10*log10(spektrumToplam / max(spektrumToplam(:)));
if size(spektrumdB,1) ~= numel(elScan), spektrumdB = spektrumdB.'; end

figure('Name', 'MUSIC Heatmap', 'Position', [150, 100, 800, 500]);
imagesc(azScan, elScan, spektrumdB); 
axis xy; colorbar; hold on;
plot(sat1Az, sat1El, "wo", 'MarkerSize', 12, 'LineWidth', 2);
plot(sat2Az, sat2El, "co", 'MarkerSize', 12, 'LineWidth', 2);
plot(spfAz,  spfEl,  "ms", 'MarkerSize', 12, 'LineWidth', 2);
plot(estAz1, estEl1, "rx", 'MarkerSize', 14, 'LineWidth', 2);
plot(estAz2, estEl2, "rx", 'MarkerSize', 14, 'LineWidth', 2);
plot(estAzSpf, estElSpf, "rx", 'MarkerSize', 14, 'LineWidth', 2);
legend({'PRN 1 (Uydu)', 'PRN 2 (Uydu)', 'PRN 1 (Spoofer)', 'MUSIC Kestirimleri'}, 'TextColor', 'black', 'Color', 'white', 'Location', 'southwest');
xlabel("Azimuth [°]"); ylabel("Elevation [°]");
title("Korelasyon Sonrası MUSIC 2D Uzaysal Spektrum Heatmap");

% 3. GRAFİK: Çoklu Hüzme LCMV Anten Diyagramları
figure('Name', 'LCMV Beamforming Pattern', 'Position', [200, 150, 1000, 400]);

subplot(1,2,1);
pattern(dizi, fL1, -180:1:180, 0:1:90, PropagationSpeed=physconst("LightSpeed"), ...
    Weights=wLCMV1, CoordinateSystem="rectangular", Type="powerdb", Normalize=true);
title("Hüzme 1 (PRN 1'e Yönelik LCMV)");
xlabel("Azimuth [°]"); ylabel("Elevation [°]");

subplot(1,2,2);
pattern(dizi, fL1, -180:1:180, 0:1:90, PropagationSpeed=physconst("LightSpeed"), ...
    Weights=wLCMV2, CoordinateSystem="rectangular", Type="powerdb", Normalize=true);
title("Hüzme 2 (PRN 2'ye Yönelik LCMV)");
xlabel("Azimuth [°]"); ylabel("Elevation [°]");