function [durum, plan] = gps_hareketli_senaryo(ayar)
%GPS_HAREKETLI_SENARYO Sürekli deneyin verici/kanal durumunu oluşturur.
% [durum, plan] = gps_hareketli_senaryo()
% [X, durum, gercek] = gps_hareketli_blok(durum)
% Alıcıya sadece X, örnekleme/IF bilgisi ve aranacak PRN kimlikleri verilir.
% durum ve gercek, alıcı kestirimine değil yalnızca simülasyona aittir.
% Açı değişimi bir kinematik yol modelidir; menzil/Doppler eşleşmesi yoktur.
    if nargin == 0, ayar = struct(); end
    varsayilan = struct('fs',16.368e6, 'fIF',4.092e6, ...
        'sigma2',1, 'seed',1, 'hareketAdimiMs',100, ...
        'baslangicBeklemeMs',200, 'araBeklemeMs',300, ...
        'sonBeklemeMs',200, 'gucGecisMs',200, 'bitisBeklemeMs',300);
    adlar = fieldnames(ayar);
    for i = 1:numel(adlar)
        assert(isfield(varsayilan, adlar{i}), 'Bilinmeyen ayar: %s', adlar{i});
        varsayilan.(adlar{i}) = ayar.(adlar{i});
    end
    ayar = varsayilan;
    sureAdlari = {'hareketAdimiMs','baslangicBeklemeMs','araBeklemeMs', ...
        'sonBeklemeMs','gucGecisMs','bitisBeklemeMs'};
    for i = 1:numel(sureAdlari)
        validateattributes(ayar.(sureAdlari{i}), {'numeric'}, ...
            {'scalar','integer','positive','finite'});
    end
    validateattributes(ayar.sigma2, {'numeric'}, {'scalar','positive','finite'});
    validateattributes(ayar.fIF, {'numeric'}, {'scalar','nonnegative','finite'});
    validateattributes(ayar.seed, {'numeric'}, {'scalar','integer','nonnegative','finite'});
    assert(ayar.fIF + 3000 + 1.023e6 < ayar.fs/2, ...
        'IF ve C/A ana lobu Nyquist aralığına sığmalıdır.');

    % Her hareket düğümü 1 derece ilerler. Düğümler arasında doğrusal hareket,
    % aynı açılı düğümler arasında duruş vardır; ayrı deneyler oluşturulmaz.
    tMs = [0; ayar.baslangicBeklemeMs];
    az = [30; 30];
    guc = [5; 5];
    for aci = 31:45
        tMs(end+1,1) = tMs(end) + ayar.hareketAdimiMs;
        az(end+1,1) = aci;
        guc(end+1,1) = 5;
        if ismember(aci, [35 40])
            tMs(end+1,1) = tMs(end) + ayar.araBeklemeMs;
            az(end+1,1) = aci;
            guc(end+1,1) = 5;
        end
    end
    tMs(end+1,1) = tMs(end) + ayar.sonBeklemeMs;
    az(end+1,1) = 45; guc(end+1,1) = 5;
    for p = [7 10]
        tMs(end+1,1) = tMs(end) + ayar.gucGecisMs;
        az(end+1,1) = 45; guc(end+1,1) = p;
    end
    tMs(end+1,1) = tMs(end) + ayar.bitisBeklemeMs;
    az(end+1,1) = 45; guc(end+1,1) = 10;
    plan = table(tMs, az, repmat(20,size(az)), guc, ...
        'VariableNames', {'Zaman_ms','Spoofer_Az_deg','Spoofer_El_deg','GucFarki_dB'});

    c = physconst('LightSpeed');
    fL1 = 1575.42e6;
    durum.ayar = ayar;
    durum.plan = plan;
    durum.dizi = phased.URA('Size',[2 2], 'ElementSpacing',c/fL1/2, 'ArrayNormal','z');
    durum.fL1 = fL1;
    durum.N = round(ayar.fs*1e-3);
    durum.epoch = 0;
    durum.toplamEpoch = tMs(end);
    durum.rastgele = RandStream('mt19937ar','Seed',ayar.seed);
    % GPS C/A dalga biçimini elle kurmuyoruz: üç toolbox verici nesnesi.
    durum.gps = {gpsWaveformGenerator('SignalType','legacy','PRNID',1, ...
        'EnablePCode',false,'SampleRate',ayar.fs), ...
        gpsWaveformGenerator('SignalType','legacy','PRNID',2, ...
        'EnablePCode',false,'SampleRate',ayar.fs), ...
        gpsWaveformGenerator('SignalType','legacy','PRNID',1, ...
        'EnablePCode',false,'SampleRate',ayar.fs)};
    durum.navdata = randi(durum.rastgele,[0 1],ceil(tMs(end)/20),2);
    durum.bb = [];
    durum.uyduAzEl = [30 120; 50 70];
    durum.uyduFd = [1500 -500];
    durum.uyduCN0 = [45 44];
    durum.spfKodChip = 3;
    durum.spfFd = 1530;
    durum.spfFaz = 0;
    durum.spfGecikme = dsp.Delay(round(durum.spfKodChip*ayar.fs/1.023e6));
    durum.uyduFazDizisi = [collectPlaneWave(durum.dizi,1, ...
        durum.uyduAzEl(:,1),fL1); collectPlaneWave(durum.dizi,1, ...
        durum.uyduAzEl(:,2),fL1)];
end
