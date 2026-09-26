%% ============================================================
%  nlmsTemizleme.m 
%  MUSIC (DOA) çıktısını uzamsal olarak demetleyip (beamforming)
%  ardından NLMS uyarlamalı filtresi ile temizleyen betik.
%% ============================================================

clc;

requiredVars = {"dizi","X","fL1","nav2mat","estAz","estEl","iSat","t","fs"};
for k = 1:numel(requiredVars)
    if ~exist(requiredVars{k}, "var")
        error("Değişken bulunamadı: '%s'. Önce signalGeneratorMusic.m çalıştırılmalı.", requiredVars{k});
    end
end


%% 1) MUSIC'in kestirdiği yönde MVDR Demetleme (Uzamsal Filtre)
% MUSIC'ten elde edilen gerçek uydu yönü Phased Array formatına çevrilir
hedefYon = [nav2mat(estAz(iSat)); estEl(iSat)];

% Hazır MVDR Beamformer paketinin (System Object) tanımlanması
mvdrBeamformer = phased.MVDRBeamformer( ...
    SensorArray=dizi, ...
    OperatingFrequency=fL1, ...
    PropagationSpeed=physconst("LightSpeed"), ...
    Direction=hedefYon);

% Ham anten matrisi (X) doğrudan MVDR filtresinden geçirilir
musicSinyal = mvdrBeamformer(X);   
musicSinyalReal = real(musicSinyal);         

%% 2) NLMS Uyarlamalı Hat Güçlendirici (Adaptive Line Enhancer - ALE)
D  = round(1e-3 * fs);      % 1 ms gecikme ≈ 1 C/A kod periyodu
M  = 64;                    % Filtre uzunluğu (daha iyi öğrenme için artırıldı)

xin = [zeros(D,1); musicSinyalReal(1:end-D)];  
d   = musicSinyalReal;                          

nlms = dsp.LMSFilter(M, Method="Normalized LMS");

mumax = maxstep(nlms, xin);
nlms.StepSize = mumax / 5; % Yakınsama hızını artırmak için adım boyu bir miktar büyütüldü

[temizSinyal, hataSinyal, agirliklar] = nlms(xin, d);   

%% 3) Ham (MUSIC demetlenmiş) sinyal ile NLMS çıktısının karşılaştırması
% Düzeltme: Çizim aralığı gecikme (D) sonrasına, filtrenin oturduğu bölgeye kaydırıldı.
N = length(musicSinyalReal);
startIdx = D + 2000; 
gosterimPenceresi = 2000;
endIdx = min(startIdx + gosterimPenceresi, N);
plotRange = startIdx:endIdx;

figure;
plot(t(plotRange)*1e6, musicSinyalReal(plotRange)); hold on;
plot(t(plotRange)*1e6, temizSinyal(plotRange), "r", LineWidth=1.5);
grid on;
legend("Ham sinyal (MUSIC demetleme)", "NLMS ile temizlenmiş sinyal");
xlabel("Zaman (\mus)");
ylabel("Genlik");
title("MUSIC Çıkışı ile NLMS Temizlenmiş Sinyalin Karşılaştırması (Yakınsama Sonrası)");


%% 4) (Opsiyonel) Simülasyon doğrulaması
if exist("sat_if", "var")
    mseHam   = mean(abs(musicSinyalReal(plotRange) - sat_if(plotRange)).^2);
    mseTemiz = mean(abs(temizSinyal(plotRange)    - sat_if(plotRange)).^2);
    fprintf("\n--- NLMS Performans Özeti ---\n");
    fprintf("Ham sinyal MSE   : %.4e\n", mseHam);
    fprintf("Temiz sinyal MSE : %.4e\n", mseTemiz);
    fprintf("İyileşme         : %.2f dB\n", 10*log10(mseHam / mseTemiz));
end