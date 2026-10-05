function [alici, cikti] = gps_alici_blok(alici, X)
%GPS_ALICI_BLOK Nedensel edinme -> takip -> korelasyon -> MUSIC -> MVDR/NLMS.
% Yeni MUSIC/MVDR kestirimleri SONRAKİ bloğa uygulanır. Gerçek senaryo verisi yok.
    N = alici.N;
    assert(isequal(size(X),[N 4]),'X bir ms ve dört anten içermelidir.');
    alici.epoch = alici.epoch+1;
    ep = alici.epoch;
    cikti = struct('mvdr',complex(nan(N,2)), 'nlms',complex(nan(N,2)), ...
        'hazir',alici.hazir,'guncellendi',false,'yonler',alici.yonler, ...
        'guc',alici.guc,'secim',alici.secim,'takip',[], 'promptGucu',nan(2,2), ...
        'wUygulanan',complex(nan(4,2)), 'wNlmsUygulanan',complex(nan(4,2)), ...
        'ssr',nan(1,3));
    if ep <= alici.edinmeMs
        alici.hamBuffer((ep-1)*N+(1:N),:) = X;
        if ep == alici.edinmeMs
            alici.gurultu = median(mean(abs(alici.hamBuffer).^2,1));
            for g = 1:2
                adaylar = gps_prn_edinme(alici.hamBuffer,alici.prnler(g), ...
                    alici.fs,alici.fIF,alici.fdGrid,2);
                for j = 1:numel(adaylar)
                    idx = numel(alici.kanallar)+1;
                    alici.kanallar{idx} = struct('grup',g,'takip',adaylar(j));
                    alici.gruplar{g}(end+1) = idx;
                    alici.C{idx} = complex(zeros(alici.pencereMs,4));
                    alici.mvdr{idx} = phased.MVDRBeamformer('SensorArray',alici.dizi, ...
                        'OperatingFrequency',alici.fL1,'TrainingInputPort',true, ...
                        'DirectionSource','Input port','WeightsOutputPort',true, ...
                        'DiagonalLoadingFactor',alici.gurultu/N);
                    alici.training{idx} = [];
                end
            end
            assert(all(~cellfun(@isempty,alici.gruplar)), ...
                'Bir PRN edinilemedi; veri penceresi/eşik gözden geçirilmeli.');
            nK = numel(alici.kanallar);
            alici.yonler = nan(2,nK);
            alici.guc = nan(1,nK);
            alici.wAday = complex(zeros(4,nK));
            alici.hamBuffer = [];
            alici.korGuc = nan(alici.pencereMs,nK,3);
            fprintf('Edinme: PRN1 %d aday, PRN2 %d aday.\n', ...
                numel(alici.gruplar{1}),numel(alici.gruplar{2}));
        end
        return;
    end

    % Önceki pencereden gelen ağırlıklar: güncel bloğun ana hüzme çıkışları.
    if alici.hazir
        for g = 1:2
            i = alici.secim(g);
            a = alici.sv(alici.fL1,alici.yonler(:,i));
            korunan = a;
            diger = setdiff(alici.gruplar{g},i);
            if ~isempty(diger)
                korunan = [korunan,alici.sv(alici.fL1,alici.yonler(:,diger))];
            elseif numel(alici.gruplar{1}) == 2
                % PRN2 hüzmesinde de saptanan güçlü PRN1 adayının bastırmasını koru.
                [~,j] = max(alici.guc(alici.gruplar{1}));
                spfIdx = alici.gruplar{1}(j);
                korunan = [korunan,alici.sv(alici.fL1,alici.yonler(:,spfIdx))];
            end
            % Ağırlık hesabı ve beamforming doğrudan toolbox nesnesindedir.
            [cikti.mvdr(:,g),cikti.wUygulanan(:,g)] = alici.mvdr{i}( ...
                X,alici.training{i},alici.yonler(:,i));
            [cikti.nlms(:,g),alici.nlms{g},cikti.wNlmsUygulanan(:,g)] = gps_uzamsal_nlms( ...
                X,cikti.wUygulanan(:,g),korunan,alici.nlms{g},cikti.mvdr(:,g));
        end
    end

    satir = mod(ep-alici.edinmeMs-1,alici.pencereMs)+1;
    bilgi = repmat(struct('P',0,'lock',0,'tauChip',0,'fd',0),1,numel(alici.kanallar));
    for i = 1:numel(alici.kanallar)
        g = alici.kanallar{i}.grup;
        yTakip = X(:,1);
        if alici.hazir
            a = alici.sv(alici.fL1,alici.yonler(:,i));
            % NCO başlangıcında anten1 fazı vardı. Birim kazançlı hüzmenin
            % geri beslemesine aynı anten fazını ekleyerek faz sıçramasını önle.
            if alici.secim(g) == i
                yTakip = cikti.nlms(:,g)*exp(1j*angle(a(1)));
            else
                yTakip = (X*conj(alici.wAday(:,i)))*exp(1j*angle(a(1)));
            end
        end
        [C,alici.kanallar{i}.takip,bilgi(i),rep] = gps_dll_pll_blok( ...
            X,yTakip,alici.kodlar{g},alici.kanallar{i}.takip,alici.fs,alici.fIF,ep);
        alici.C{i}(satir,:) = C;
        if alici.hazir
            deger = abs(sum([X(:,1),cikti.mvdr(:,g),cikti.nlms(:,g)].*conj(rep),1)/N).^2;
            alici.korGuc(satir,i,:) = reshape(deger,1,1,3);
        end
        if alici.hazir && alici.secim(g) == i
            cikti.promptGucu(g,:) = abs([sum(cikti.mvdr(:,g).*conj(rep)), ...
                sum(cikti.nlms(:,g).*conj(rep))]/N).^2;
        end
    end
    cikti.takip = bilgi;

    if ep >= alici.edinmeMs+alici.pencereMs && ...
            mod(ep-alici.edinmeMs,alici.guncellemeMs) == 0
        for g = 1:2
            ids = alici.gruplar{g};
            for i = ids
                Ci = alici.C{i};
                enerji = mean(abs(Ci(:)).^2);
                alici.guc(i) = max(enerji-alici.gurultu/N,eps);
                % Aynı PRN'nin iki kod tepesini ayrı kapılarda tut. Ortak
                % MUSIC2 taramasında güçlü tepenin yanındaki iki maksimumun
                % iki kaynak sanılması, zayıf adaya yanlış yön atayabiliyor.
                % Burada her edinilmiş adayın kendi baskın yönü kestirilir.
                % Bu yöntem kod fazıyla ayrılan kaynaklar içindir; çakışık
                % coherent iki kaynak için genel çözüm iddiası yoktur.
                est = alici.music{g,1};
                [spektrum,doa] = est(Ci);
                alici.musicSpektrum{i} = spektrum;
                alici.yonler(:,i) = doa(:,1);
            end
            % AÇIK DENEY VARSAYIMI: aynı PRN'de daha zayıf aday otantiktir.
            [~,j] = min(alici.guc(ids));
            yeniSecim = ids(j);
            eski = alici.secim(g);
            if eski ~= 0 && eski ~= yeniSecim && ...
                    10*log10(alici.guc(eski)/alici.guc(yeniSecim)) < 2
                yeniSecim = eski; % 2 dB histerezis: gürültüyle seçim salınmasın.
            end
            if eski ~= yeniSecim, alici.nlms{g} = []; end
            alici.secim(g) = yeniSecim;
        end
        for i = 1:numel(alici.kanallar)
            % MVDR eğitim kovaryansı: diğer adayların PRN korelasyonları.
            % Zayıf RF kaynaklar için ham X kovaryansı yerine despreading kullanılır.
            egitim = complex(zeros(0,4));
            for j = setdiff(1:numel(alici.kanallar),i)
                egitim = [egitim;alici.C{j}]; %#ok<AGROW>
            end
            alici.training{i} = sqrt(numel(alici.kanallar)-1)*egitim;
            [~,alici.wAday(:,i)] = alici.mvdr{i}( ...
                complex(zeros(1,4)),alici.training{i},alici.yonler(:,i));
        end
        for g = 1:2, alici.w(:,g) = alici.wAday(:,alici.secim(g)); end
        alici.hazir = true;
        alici.sonUzamsalMs = ep;
        cikti.guncellendi = true;
    end
    cikti.yonler = alici.yonler;
    cikti.guc = alici.guc;
    cikti.secim = alici.secim;
    ids = alici.gruplar{1};
    if cikti.hazir && numel(ids)==2
        zayif = cikti.secim(1); guclu = setdiff(ids,zayif);
        Pz = reshape(mean(alici.korGuc(:,zayif,:),1,'omitnan'),1,3);
        Pg = reshape(mean(alici.korGuc(:,guclu,:),1,'omitnan'),1,3);
        cikti.ssr = 10*log10(max(Pg,eps)./max(Pz,eps));
    end
end
