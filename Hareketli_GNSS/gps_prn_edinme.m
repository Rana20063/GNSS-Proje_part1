function adaylar = gps_prn_edinme(X, prn, fs, fIF, fdGrid, maxAday)
%GPS_PRN_EDINME Toolbox haritasından aynı PRN'nin en fazla iki tepesini seç.
% Tablodaki tek tepe yerine ikinci çıktı corrmat kullanılır; aynı PRN'li
% uydu ve spoofer böylece ayrı aday olabilir. Kod fazı ve Doppler kestirilir.
    N = round(fs*1e-3); E = floor(size(X,1)/N);
    acquirer = gnssSignalAcquirer('GNSSSignalType','GPS C/A','SampleRate',fs, ...
        'IntermediateFrequency',fIF,'FrequencyRange',[fdGrid(1) fdGrid(end)], ...
        'FrequencyResolution',fdGrid(2)-fdGrid(1));
    S = zeros(N,numel(fdGrid));
    for ep = 1:E
        idx = (ep-1)*N+(1:N);
        for anten = 1:4
            [~,harita] = acquirer(X(idx,anten),prn);
            S = S+harita.^2; % Epoch ve antenlerde noncoherent güç toplama.
        end
    end
    taban = median(S(:)); calisma = S;
    kod = gps_ca_replika(prn,fs);
    t = (0:E*N-1).'/fs;
    ornChip = fs/1.023e6;
    adaylar = struct('tauChip',{},'fd',{},'phase',{},'metric',{}, ...
        'fdBase',{},'carrierFd',{},'phaseNCO',{},'age',{},'lock',{},'tracker',{});
    for k = 1:min(maxAday,2)
        [tepe,lin] = max(calisma(:));
        if tepe/max(taban,eps)<8, break; end
        [itau,ifd] = ind2sub(size(calisma),lin);
        tau = itau-1; rep = circshift(kod,tau);
        C = complex(zeros(E,4));
        for ep = 1:E
            idx = (ep-1)*N+(1:N);
            z = X(idx,:).*exp(-1j*2*pi*(fIF+fdGrid(ifd))*t(idx));
            C(ep,:) = sum(z.*rep,1)/N;
        end
        % Kaba Doppler düzeltmesi alınan korelasyon fazından ölçülür.
        fark = C(2:end,:).*conj(C(1:end-1,:));
        fd = fdGrid(ifd)+angle(sum(fark(:).^2))/(4*pi*1e-3);
        chip0 = mod(round(tau/ornChip),1023);
        tracker = gnssSignalTracker('GNSSSignalType','GPS C/A', ...
            'SampleRate',fs,'PRNID',prn,'InitialCodePhaseOffset',chip0, ...
            'InitialFrequencyOffset',fd,'PLLNoiseBandwidth',18, ...
            'DisableFLL',true,'DLLNoiseBandwidth',2);
        adaylar(end+1) = struct('tauChip',chip0,'fd',fd, ...
            'phase',mod(2*pi*fd*E*1e-3,2*pi),'metric',tepe/max(taban,eps), ...
            'fdBase',fd,'carrierFd',fd,'phaseNCO',0,'age',0,'lock',0,'tracker',tracker); %#ok<AGROW>
        disla = mod(tau+(-round(ornChip):round(ornChip)),N)+1;
        calisma(disla,:) = 0;
    end
end
