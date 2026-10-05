function gps_sonuc_goster(s)
    figure('Name','Hareketli GNSS - Geri Beslemeli Alıcı','Position',[80 80 1300 750]);
    tiledlayout(3,2,'TileSpacing','compact');
    nexttile; plot(s.zamanMs,s.gercekSpf(:,1),'k',s.zamanMs,s.yonSpf(:,1),'r');
    grid on; ylabel('Azimut [derece]'); legend('Gerçek spoofer','MUSIC güçlü aday');
    nexttile; plot(s.zamanMs,s.gercekSpf(:,3),'k',s.zamanMs,s.gucFarkiOlculen,'r');
    grid on; ylabel('Güç farkı [dB]'); legend('Gerçek','Alıcı kestirimi');
    nexttile; plot(s.zamanMs,[s.yonUydu1(:,1),s.yonUydu2(:,1)]);
    grid on; ylabel('Seçilen yön [derece]'); legend('PRN1','PRN2');
    nexttile; plot(s.zamanMs,s.sonum); grid on;
    ylabel('MVDR hüzme sönümü [dB]'); title('Yönsel kazanç oranı; çıkış SSR ölçümü değil');
    nexttile; plot(s.zamanMs,s.takipTau); grid on;
    xlabel('Zaman [ms]'); ylabel('Takip gecikmesi [chip]'); legend('PRN1','PRN2');
    nexttile; plot(s.zamanMs,s.kilit); grid on;
    xlabel('Zaman [ms]'); ylabel('PLL faz uyumu'); legend('PRN1','PRN2'); ylim([-1 1]);
    figure('Name','İki Ayrı MVDR / NLMS Çıkışı - Son Kayıt','Position',[100 100 1300 600]);
    tiledlayout(2,1,'TileSpacing','compact');
    goster = max(1,size(s.mvdr,1)-1500):size(s.mvdr,1);
    for g = 1:2
        nexttile;
        plot(s.kayitZamanS(goster)*1e3,real(s.mvdr(goster,g))); hold on;
        plot(s.kayitZamanS(goster)*1e3,real(s.nlms(goster,g)));
        grid on; ylabel('Gerçel IF genliği');
        title(sprintf('PRN %d: MVDR%d -> NLMS%d',g,g,g));
        legend('MVDR','NLMS');
    end
    xlabel('Zaman [ms]');
    figure('Name','Başarı Ölçümü - Kod Tepeleri','Position',[120 120 1300 600]);
    tiledlayout(1,2,'TileSpacing','compact');
    renkler=[0.6 0.6 0.6; 0 0.447 0.741; 0.85 0.325 0.098];
    nexttile; h=plot(s.zamanMs,s.ssr); grid on;
    for k=1:3, h(k).Color=renkler(k,:); end
    yline(0,'--k'); xlabel('Zaman [ms]'); ylabel('Spoofer / uydu tepe gücü [dB]');
    legend('Ham anten','MVDR1','NLMS1','Location','best');
    title('Daha düşük: spoofer tepesi daha zayıf');
    nexttile;
    [chip,P] = son_kod_korelasyonu(s);
    h=plot(chip,10*log10(max(P,eps))); grid on;
    for k=1:3, h(k).Color=renkler(k,:); end
    xline(0,'--k','Uydu','HandleVisibility','off');
    xline(3,'--m','Spoofer','HandleVisibility','off');
    xlabel('Kod gecikmesi [chip]'); ylabel('Korelasyon gücü [dB]');
    legend('Ham anten','MVDR1','NLMS1','Location','best');
    title('Son kayıtta PRN1 kod korelasyonu');
    gps_son_uzamsal_grafikler(s);
end

function [chip,P] = son_kod_korelasyonu(s)
% Yerel kod ve ALINAN veriden kestirilen taşıyıcı; verici replikası yok.
    N = round(s.fs*1e-3); ornChip=s.fs/1.023e6;
    chip=(-5:1/ornChip:8).'; lag=mod(round(chip*ornChip),N)+1;
    kod=gps_ca_replika(1,s.fs); K=conj(fft(kod));
    sinyaller=[s.ham,s.mvdr(:,1),s.nlms(:,1)];
    fd=mean(s.takipFd(end-min(99,size(s.takipFd,1)-1):end,1),'omitnan');
    P=zeros(numel(chip),3);
    for ep=1:floor(size(sinyaller,1)/N)
        idx=(ep-1)*N+(1:N);
        bb=sinyaller(idx,:).*exp(-1j*2*pi*(s.fIF+fd)*s.kayitZamanS(idx));
        c=ifft(fft(bb,[],1).*K,[],1)/N;
        P=P+abs(c(lag,:)).^2;
    end
    P=P/floor(size(sinyaller,1)/N);
end
