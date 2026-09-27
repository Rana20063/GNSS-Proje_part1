clear;
clc;
close all;

%% GPS L1 C/A - 2 Gerçek Uydu + 1 Spoofer - 4 Elemanlı CRPA + MUSIC + Multi-Beam MVDR

numbits = 5;            % 5 navigasyon biti = 100 ms
navdata = randi([0 1], numbits, 1);

fs    = 16.368e6;       % örnekleme hızı
fIF   = 4.092e6;        % ara frekans
fL1   = 1575.42e6;      % L1 taşıyıcı
fChip = 1.023e6;        % C/A chip hızı
sigma2 = 1;             % eleman başına gürültü gücü

% --- 1. Gerçek Uydu (PRN 1) ---
prn1     = 1;
sat1Az   = 60;          % [derece]
sat1El   = 50;          % [derece]
sat1Fd   = 1500;        % Doppler [Hz]
sat1CN0  = 50;          % [dB-Hz]

% --- 2. Gerçek Uydu (PRN 2) ---
prn2     = 2;
sat2Az   = 120;         % [derece]
sat2El   = 70;          % [derece]
sat2Fd   = -500;        % Doppler [Hz]
sat2CN0  = 44;          % [dB-Hz]

% --- Spoofer: PRN 1'i taklit eder, daha güçlüdür ---
spfAz        = 65;      % PRN 1'e açısal olarak yakın [derece]
spfEl        = 55;      % [derece] 
spfFdFark    = 30;      % [Hz]
spfKodKayma  = 3;       % [chip]
spfGucFark   = 10;      % 1. uydudan ne kadar güçlü [dB]

%% GPS waveform generator (Temel Bant Üretimi)

% PRN 1 için temel bant
gpswaveobj1 = gpsWaveformGenerator(SignalType="legacy", PRNID=prn1, EnablePCode=false, SampleRate=fs);
baseband1 = gpswaveobj1(navdata);

% PRN 2 için temel bant
gpswaveobj2 = gpsWaveformGenerator(SignalType="legacy", PRNID=prn2, EnablePCode=false, SampleRate=fs);
baseband2 = gpswaveobj2(navdata);

t = (0:length(baseband1)-1).' / fs;

%% Sinyalleri IF'e taşı (Doppler dahil)

sat1_if_c = baseband1 .* exp(1j*2*pi*(fIF + sat1Fd)*t);     
sat2_if_c = baseband2 .* exp(1j*2*pi*(fIF + sat2Fd)*t);     

% Spoofer (PRN 1 kopyası)
kaymaOrnek = round(spfKodKayma * fs / fChip);
spf_bb     = circshift(baseband1, kaymaOrnek);
spf_if_c   = spf_bb .* exp(1j*(2*pi*(fIF + sat1Fd + spfFdFark)*t + 2*pi*rand));

%% 4 elemanlı CRPA dizisi 

lambda = physconst("LightSpeed") / fL1;          
dizi = phased.URA(Size=[2 2], ElementSpacing=lambda/2, ArrayNormal="z");

% Açı kuralı dönüşümü (Navigasyon <-> Phased Array)
nav2mat = @(az) mod(90 - az + 180, 360) - 180;
mat2nav = @(az) mod(90 - az, 360);

%% Sinyallerin Anten Elemanlarında Toplanması

Ps = mean(abs(baseband1).^2);
genlik = @(CN0dB) sqrt(10^(CN0dB/10) * sigma2 / (fs * Ps));

Asat1 = genlik(sat1CN0);
Asat2 = genlik(sat2CN0);
Aspf  = genlik(sat1CN0 + spfGucFark);

Xsat1 = collectPlaneWave(dizi, Asat1*sat1_if_c, [nav2mat(sat1Az); sat1El], fL1);
Xsat2 = collectPlaneWave(dizi, Asat2*sat2_if_c, [nav2mat(sat2Az); sat2El], fL1);
Xspf  = collectPlaneWave(dizi, Aspf*spf_if_c,   [nav2mat(spfAz); spfEl], fL1);

% Toplam sinyal matrisi (N anten x K örnek)
X = Xsat1 + Xsat2 + Xspf + sqrt(sigma2/2) * (randn(size(Xsat1)) + 1j*randn(size(Xsat1)));

%% MUSIC 2D Kestirimi (3 Kaynak)

azScan = -180:1:180;
elScan = 0:1:90;

musicEst = phased.MUSICEstimator2D( ...
    SensorArray=dizi, ...
    OperatingFrequency=fL1, ...
    AzimuthScanAngles=azScan, ...
    ElevationScanAngles=elScan, ...
    DOAOutputPort=true, ...
    NumSignalsSource="Property", ...
    NumSignals=3); % Artık 3 kaynak aranıyor

[spektrum, aci] = musicEst(X);

%% Kestirimleri Gerçek Yönlerle Eşleştirme (Uzaklık Matrisi ile)

u = @(az,el) [cosd(el)*sind(az); cosd(el)*cosd(az); sind(el)];   
acisalFark = @(az1,el1,az2,el2) acosd(min(1, dot(u(az1,el1), u(az2,el2))));

estAz = mat2nav(aci(1,:));
estEl = aci(2,:);

gercekler = [sat1Az, sat1El; sat2Az, sat2El; spfAz, spfEl];
eslesme = zeros(1,3); 
atanan = [];

% Her bir gerçek hedef için en yakın MUSIC kestirimini bul
for i = 1:3
    minMesafe = inf;
    enIyiIndis = -1;
    for j = 1:3
        if ~ismember(j, atanan)
            d = acisalFark(gercekler(i,1), gercekler(i,2), estAz(j), estEl(j));
            if d < minMesafe
                minMesafe = d;
                enIyiIndis = j;
            end
        end
    end
    eslesme(i) = enIyiIndis;
    atanan = [atanan, enIyiIndis];
end

iSat1 = eslesme(1);
iSat2 = eslesme(2);
iSpf  = eslesme(3);

fprintf("\nMUSIC KESTİRİMLERİ:\n");
fprintf("1. Uydu (PRN 1): Az: %.1f°, El: %.1f°\n", estAz(iSat1), estEl(iSat1));
fprintf("2. Uydu (PRN 2): Az: %.1f°, El: %.1f°\n", estAz(iSat2), estEl(iSat2));
fprintf("Spoofer        : Az: %.1f°, El: %.1f°\n\n", estAz(iSpf), estEl(iSpf));

%% ÇOKLU HÜZME (MULTI-BEAM) MVDR UYGULAMASI

% --- HÜZME 1: 1. Uyduya (PRN 1) Yönelik MVDR ---
hedefYon1 = [nav2mat(estAz(iSat1)); estEl(iSat1)];
mvdr1 = phased.MVDRBeamformer( ...
    SensorArray=dizi, ...
    OperatingFrequency=fL1, ...
    Direction=hedefYon1, ...
    WeightsOutputPort=true);

[yMVDR1, wMVDR1] = mvdr1(X);

% --- HÜZME 2: 2. Uyduya (PRN 2) Yönelik MVDR ---
hedefYon2 = [nav2mat(estAz(iSat2)); estEl(iSat2)];
mvdr2 = phased.MVDRBeamformer( ...
    SensorArray=dizi, ...
    OperatingFrequency=fL1, ...
    Direction=hedefYon2, ...
    WeightsOutputPort=true);

[yMVDR2, wMVDR2] = mvdr2(X);

fprintf("Çoklu MVDR Çıkışları Hazır.\n");
fprintf("1. Hüzme (PRN 1) çıkış boyutu: %d x %d\n", size(yMVDR1,1), size(yMVDR1,2));
fprintf("2. Hüzme (PRN 2) çıkış boyutu: %d x %d\n", size(yMVDR2,1), size(yMVDR2,2));

%% MUSIC Spektrumu 2D Heatmap ve Kestirimlerin Çizdirilmesi

% Spektrum boyutlarını kontrol et ve imagesc eksenlerine uydurmak için transpoze al
if size(spektrum,1) ~= numel(elScan)
    spektrum = spektrum.'; 
end

% Spektrumu dB formatına çevir ve en yüksek tepe noktasına göre normalize et
spektrumdB = 10*log10(spektrum / max(spektrum(:)));

figure('Position', [150, 150, 800, 500]);
imagesc(azScan, elScan, spektrumdB); 
axis xy; 
c = colorbar;
c.Label.String = 'Normalize Güç [dB]';
hold on;

% Gerçek konumların çizdirilmesi (Navigasyon açısından Phased Array açısına dönüştürülerek)
p1 = plot(nav2mat(sat1Az), sat1El, "wo", 'MarkerSize', 12, 'LineWidth', 2);
p2 = plot(nav2mat(sat2Az), sat2El, "co", 'MarkerSize', 12, 'LineWidth', 2);
p3 = plot(nav2mat(spfAz),  spfEl,  "ms", 'MarkerSize', 12, 'LineWidth', 2);

% MUSIC kestirimlerinin çizdirilmesi (aci matrisinin 1. satırı azimut, 2. satırı elevasyon)
p4 = plot(aci(1,:), aci(2,:), "rx", 'MarkerSize', 12, 'LineWidth', 2);

% Açıklamalar ve etiketler
legend([p1, p2, p3, p4], ...
       {'Gerçek Uydu 1', 'Gerçek Uydu 2', 'Gerçek Spoofer', 'MUSIC Kestirimleri'}, ...
       'TextColor', 'black', 'Color', 'white', 'Location', 'southwest');
       
xlabel("Azimut [°] (Phased Array Kuralı: Doğu'dan saat yönünün tersine)");
ylabel("Elevasyon [°]");
title("MUSIC 2D Uzaysal Spektrum Heatmap (2 Uydu + 1 Spoofer)");
%% Anten Diyagramlarının Çizdirilmesi (İki Bağımsız Hüzme)

figure('Position', [100, 100, 1000, 400]);

% 1. Uydu için Pattern
subplot(1,2,1);
pattern(dizi, fL1, -180:1:180, 0:1:90, ...
    PropagationSpeed=physconst("LightSpeed"), ...
    Weights=wMVDR1, ...
    CoordinateSystem="rectangular", ...
    Type="powerdb", ...
    Normalize=true);
title("Hüzme 1 (PRN 1'e Yönelik)");
xlabel("Azimut (Phased Array)"); ylabel("Elevasyon"); zlabel("dB");

% 2. Uydu için Pattern
subplot(1,2,2);
pattern(dizi, fL1, -180:1:180, 0:1:90, ...
    PropagationSpeed=physconst("LightSpeed"), ...
    Weights=wMVDR2, ...
    CoordinateSystem="rectangular", ...
    Type="powerdb", ...
    Normalize=true);
title("Hüzme 2 (PRN 2'ye Yönelik)");
xlabel("Azimut (Phased Array)"); ylabel("Elevasyon"); zlabel("dB");