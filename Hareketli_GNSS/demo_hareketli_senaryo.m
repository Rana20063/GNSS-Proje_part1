% İlk aşama: hareketli verici/kanal modeli. Henüz alıcı/filtre zinciri değil.
% Bütün X kaydını belleğe almak yerine 1 ms bloklar halinde üretir.
clear; clc;
[durum, plan] = gps_hareketli_senaryo();
disp(plan);
Nepoch = durum.toplamEpoch;
zaman = (0:Nepoch-1).' + 0.5;
azimut = interp1(plan.Zaman_ms,plan.Spoofer_Az_deg,zaman);
gucFarki = interp1(plan.Zaman_ms,plan.GucFarki_dB,zaman);
figure('Name','Hareketli Spoofer - Tek Zaman Akışı');
tiledlayout(2,1);
nexttile; plot(zaman,azimut,'LineWidth',1.5); grid on;
ylabel('Spoofer azimut [derece]'); title('1 derece aralıklı hareket düğümleri ve duruşlar');
nexttile; plot(zaman,gucFarki,'LineWidth',1.5); grid on;
xlabel('Zaman [ms]'); ylabel('Uyduya göre güç farkı [dB]');

% Bu döngünün içine edinme/takip, korelasyon, MUSIC ve ağırlık güncellemesi
% eklenecek. gercek değişkeni yalnızca sonuç doğrulamasına ayrılmıştır.
antenGucu = zeros(Nepoch,4);
for ep = 1:Nepoch
    [X,durum,gercek] = gps_hareketli_blok(durum);
    antenGucu(ep,:) = mean(abs(X).^2,1);
    if mod(ep,200) == 0
        fprintf('%d / %d ms üretildi.\n',ep,Nepoch);
    end
end
fprintf('Tek akış tamamlandı: %.2f saniye. Alıcı henüz bağlanmadı.\n',Nepoch/1000);
