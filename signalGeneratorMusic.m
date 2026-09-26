clear;
clc;
close all;

%% GPS L1 C/A - Gerçek Uydu + Spoofer - 4 Elemanlı CRPA + MUSIC

prn = 1;                % GPS PRN 1
numbits = 5;            % 5 navigasyon biti = 100 ms
navdata = randi([0 1], numbits, 1);

fs    = 16.368e6;       % örnekleme hızı
fIF   = 4.092e6;        % ara frekans
fL1   = 1575.42e6;      % L1 taşıyıcı (anten fazları RF dalga boyuna göre oluşur)
fChip = 1.023e6;        % C/A chip hızı
sigma2 = 1;             % eleman başına gürültü gücü

% --- Gerçek uydu (navigasyon kuralı: azimut kuzeyden saat yönünde) ---
satAz  = 30;            % [derece]
satEl  = 50;            % [derece]
satFd  = 1500;          % Doppler [Hz]
satCN0 = 45;            % [dB-Hz]

% --- Spoofer: aynı PRN'i taklit eder, yerden/alçaktan ve daha güçlü gelir ---
spfAz        = 200;     % [derece]
spfEl        = 10;      % [derece] (yer tabanlı spoofer genelde ufka yakındır)
spfFdFark    = 30;      % gerçek uydunun Doppler'ine göre fark [Hz]
spfKodKayma  = 3;       % gerçek uydunun koduna göre gecikme [chip]
spfGucFark   = 10;      % gerçek uydudan ne kadar güçlü [dB]

%% GPS waveform generator (C/A + nav verisi, temel bant)

gpswaveobj = gpsWaveformGenerator(SignalType="legacy", PRNID=prn, EnablePCode=false, SampleRate=fs);

baseband = gpswaveobj(navdata);

disp("Temel bant boyutu:");
disp(size(baseband));

t = (0:length(baseband)-1).' / fs;

%% Gerçek uydu sinyali - IF'e taşı (Doppler dahil)

sat_if_c = baseband .* exp(1j*2*pi*(fIF + satFd)*t);     % kompleks IF
sat_if   = real(sat_if_c);                                % gerçel IF (ADC çıkışı gibi)

%% Spoofer sinyali - aynı kod ve aynı nav verisi, kaydırılmış
% Spoofer gerçek sinyalin kopyasını üretir: kod fazı spfKodKayma chip geride,
% Doppler'i biraz farklı, taşıyıcı fazı bağımsız.

kaymaOrnek = round(spfKodKayma * fs / fChip);
spf_bb     = circshift(baseband, kaymaOrnek);
spf_if_c   = spf_bb .* exp(1j*(2*pi*(fIF + satFd + spfFdFark)*t + 2*pi*rand));

%% Zaman domeninde gerçek uydu sinyali (kısa bir kesit)

figure;
plot(t(1:2000)*1e6, sat_if(1:2000));
grid on;
xlabel("Zaman (\mus)");
ylabel("Genlik");
title("GPS L1 C/A Sinyali - PRN 1 - IF (4.092 MHz)");

%% Spektrum (gerçek uydu)

txscope = spectrumAnalyzer( ...
    SampleRate=fs, ...
    SpectrumType="power-density", ...
    SpectrumUnits="dBW/Hz");

txscope(sat_if);

%% 4 elemanlı CRPA dizisi (2x2 kare, gökyüzüne bakıyor)
% ArrayNormal="z": elemanlar yatay düzlemde (x = Doğu, y = Kuzey, z = Yukarı)

lambda = physconst("LightSpeed") / fL1;          % ~19 cm
dizi = phased.URA(Size=[2 2], ElementSpacing=lambda/2, ArrayNormal="z");

% Açı kuralı dönüşümü: navigasyon azimutu <-> Phased Array azimutu
nav2mat = @(az) mod(90 - az + 180, 360) - 180;
mat2nav = @(az) mod(90 - az, 360);

%% Sinyallerin 4 anten elemanına ulaşması

Ps = mean(abs(baseband).^2);
genlik = @(CN0dB) sqrt(10^(CN0dB/10) * sigma2 / (fs * Ps));

Asat = genlik(satCN0);
Aspf = genlik(satCN0 + spfGucFark);

Xsat = collectPlaneWave(dizi, Asat*sat_if_c, [nav2mat(satAz); satEl], fL1);
Xspf = collectPlaneWave(dizi, Aspf*spf_if_c, [nav2mat(spfAz); spfEl], fL1);

X = Xsat + Xspf + sqrt(sigma2/2) * (randn(size(Xsat)) + 1j*randn(size(Xsat)));

%% Kovaryans özdeğerleri (2 tanesi diğerlerinden büyük olmalı)

R = (X' * X) / size(X,1);
disp("Kovaryans özdeğerleri:");
disp(sort(real(eig(R))).');

%% MUSIC 2D ile geliş açısı kestirimi (2 kaynak)

azScan = -180:1:180;
elScan = 0:1:90;

musicEst = phased.MUSICEstimator2D( ...
    SensorArray=dizi, ...
    OperatingFrequency=fL1, ...
    AzimuthScanAngles=azScan, ...
    ElevationScanAngles=elScan, ...
    DOAOutputPort=true, ...
    NumSignalsSource="Property", ...
    NumSignals=2);

[spektrum, aci] = musicEst(X);

%% Kestirimleri gerçek yönlerle eşleştir ve karşılaştır

u = @(az,el) [cosd(el)*sind(az); cosd(el)*cosd(az); sind(el)];   % Doğu-Kuzey-Yukarı
acisalFark = @(az1,el1,az2,el2) acosd(min(1, dot(u(az1,el1), u(az2,el2))));

estAz = mat2nav(aci(1,:));
estEl = aci(2,:);

% Her kestirimin gerçek uyduya mı spoofer'a mı ait olduğunu en yakın yönle belirle
d11 = acisalFark(satAz,satEl, estAz(1),estEl(1)) + acisalFark(spfAz,spfEl, estAz(2),estEl(2));
d12 = acisalFark(satAz,satEl, estAz(2),estEl(2)) + acisalFark(spfAz,spfEl, estAz(1),estEl(1));
if d11 <= d12, iSat = 1; iSpf = 2; else, iSat = 2; iSpf = 1; end

fprintf("\n                      Azimut   Elevasyon  |  Kestirim Az   El    |  Hata\n");
fprintf("Gerçek uydu (PRN %d): %6.1f°   %6.1f°    |  %8.1f°  %5.1f°  |  %5.2f°\n", ...
    prn, satAz, satEl, estAz(iSat), estEl(iSat), acisalFark(satAz,satEl,estAz(iSat),estEl(iSat)));
fprintf("Spoofer            : %6.1f°   %6.1f°    |  %8.1f°  %5.1f°  |  %5.2f°\n", ...
    spfAz, spfEl, estAz(iSpf), estEl(iSpf), acisalFark(spfAz,spfEl,estAz(iSpf),estEl(iSpf)));

%% MUSIC spektrumu + gerçek ve kestirilen yönler

if size(spektrum,1) ~= numel(elScan), spektrum = spektrum.'; end
spektrumdB = 10*log10(spektrum / max(spektrum(:)));

figure;
imagesc(azScan, elScan, spektrumdB); axis xy; colorbar;
hold on;
plot(nav2mat(satAz), satEl, "wo", MarkerSize=12, LineWidth=2);
plot(nav2mat(spfAz), spfEl, "ms", MarkerSize=12, LineWidth=2);
plot(aci(1,:), aci(2,:), "rx", MarkerSize=12, LineWidth=2);
legend("Gerçek uydu", "Spoofer", "MUSIC kestirimleri", TextColor="w", Color="none");
xlabel("Azimut [°] (Phased Array kuralı: Doğu'dan saat yönü tersi)");
ylabel("Elevasyon [°]");
title("MUSIC 2D uzamsal spektrum [dB] - gerçek uydu + spoofer");
