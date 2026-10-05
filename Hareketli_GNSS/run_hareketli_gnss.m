function sonuc = run_hareketli_gnss(secenek)
%RUN_HAREKETLI_GNSS İki uydu + hareketli spoofer, iki ayrı MVDR/NLMS çıkışı.
% sonuc = run_hareketli_gnss();
% sonuc = run_hareketli_gnss(struct('maxEpoch',200,'grafik',false));
% Senaryo ayarları: secenek.senaryo = struct('hareketAdimiMs',100,...).
% Alıcı için yerel PRN kodu dışında vericinin referans sinyali kullanılmaz.
    %% 1. Deney seçenekleri
    if nargin == 0, secenek = struct(); end
    if ~isfield(secenek,'maxEpoch'), secenek.maxEpoch = inf; end
    if ~isfield(secenek,'grafik'), secenek.grafik = true; end
    if ~isfield(secenek,'kayitMs'), secenek.kayitMs = 20; end
    if ~isfield(secenek,'senaryo'), secenek.senaryo = struct(); end
    validateattributes(secenek.kayitMs,{'numeric'},{'scalar','integer','positive'});
    validateattributes(secenek.maxEpoch,{'numeric'},{'scalar','positive'});
    %% 2. Toolbox vericileri ve toolbox alıcısı
    [verici,plan] = gps_hareketli_senaryo(secenek.senaryo);
    alici = gps_alici_baslat(verici.dizi,verici.ayar.fs,verici.ayar.fIF);
    E = min(verici.toplamEpoch,floor(secenek.maxEpoch));
    assert(E >= alici.edinmeMs+alici.pencereMs+1, ...
        'En az 61 ms çalıştırılmalıdır.');
    kayitMs = min(secenek.kayitMs,E);
    N = verici.N;
    sonuc.plan = plan;
    sonuc.surum = 'toolbox-v1';
    sonuc.fs = verici.ayar.fs;
    sonuc.fIF = verici.ayar.fIF;
    sonuc.zamanMs = (0:E-1).'+0.5;
    sonuc.gercekSpf = nan(E,3);
    sonuc.yonUydu1 = nan(E,2); sonuc.yonUydu2 = nan(E,2);
    sonuc.yonSpf = nan(E,2);
    sonuc.gucFarkiOlculen = nan(E,1);
    sonuc.sonum = nan(E,1); % Yalnızca MVDR'nin sabit uzamsal ağırlık cevabı.
    sonuc.ssr = nan(E,3);   % Ham / MVDR1 / NLMS1 ölçülen tepe güç oranları.
    sonuc.takipTau = nan(E,2); sonuc.takipFd = nan(E,2);
    sonuc.kilit = nan(E,2);
    sonuc.promptGucu = nan(E,2,2);
    sonuc.mvdr = complex(nan(N*kayitMs,2));
    sonuc.nlms = complex(nan(N*kayitMs,2));
    sonuc.ham = complex(nan(N*kayitMs,1));
    sonuc.temizUydu = complex(nan(N*kayitMs,2));
    sonuc.kayitZamanS = ((E-kayitMs)*N+(0:N*kayitMs-1)).'/sonuc.fs;
    %% 3. Tek zaman akışı: sinyal üret -> alıcıyı bir ms ilerlet
    for ep = 1:E
        [X,verici,gercek] = gps_hareketli_blok(verici);
        % Kritik sınır: alıcı gercek/verici durumuna erişmez.
        [alici,cikti] = gps_alici_blok(alici,X);
        sonuc.gercekSpf(ep,:) = [gercek.spfAz,gercek.spfEl,gercek.spfGucFarki_dB];
        if cikti.hazir
            for g = 1:2
                i = cikti.secim(g);
                if g == 1, sonuc.yonUydu1(ep,:) = cikti.yonler(:,i).';
                else, sonuc.yonUydu2(ep,:) = cikti.yonler(:,i).'; end
                sonuc.takipTau(ep,g) = cikti.takip(i).tauChip;
                sonuc.takipFd(ep,g) = cikti.takip(i).fd;
                sonuc.kilit(ep,g) = cikti.takip(i).lock;
                sonuc.promptGucu(ep,g,:) = reshape(cikti.promptGucu(g,:),1,1,2);
            end
            ids = alici.gruplar{1};
            if numel(ids) == 2
                diger = setdiff(ids,cikti.secim(1));
                sonuc.yonSpf(ep,:) = cikti.yonler(:,diger).';
                sonuc.gucFarkiOlculen(ep) = 10*log10(cikti.guc(diger)/cikti.guc(cikti.secim(1)));
            end
            % Gerçek açı SADECE başarı ölçümünde kullanılır; ağırlık hesabında yok.
            aD = alici.sv(alici.fL1,gercek.uyduAzEl(:,1));
            aS = alici.sv(alici.fL1,[gercek.spfAz;gercek.spfEl]);
            w = cikti.wUygulanan(:,1);
            sonuc.sonum(ep) = 10*log10(abs(w'*aD)^2/max(abs(w'*aS)^2,eps));
            sonuc.ssr(ep,:) = cikti.ssr;
        end
        if ep > E-kayitMs
            idx = (ep-(E-kayitMs)-1)*N+(1:N);
            sonuc.mvdr(idx,:) = cikti.mvdr;
            sonuc.nlms(idx,:) = cikti.nlms;
            sonuc.ham(idx) = X(:,1);
            sonuc.temizUydu(idx,:) = gercek.temizUydu;
        end
        if mod(ep,200) == 0
            fprintf('%d/%d ms | MUSIC güncelleme: %d | PRN1 kilit: %.2f\n', ...
                ep,E,cikti.guncellendi,sonuc.kilit(ep,1));
        end
    end
    %% 4. Sonuçlar ve grafikler (alıcı kararları için kullanılmaz)
    sonuc.alici = alici;
    sonuc.gercekUydu = gercek.uyduAzEl;
    sonuc.varsayim = 'Aynı PRN içinde daha güçlü aday spoofer kabul edilir.';
    if secenek.grafik, gps_sonuc_goster(sonuc); end
end
