function [fMusic,fMvdr] = gps_son_uzamsal_grafikler(s)
%GPS_SON_UZAMSAL_GRAFIKLER Son MUSIC güncellemesi ve iki MVDR ağırlık seti.
% Gerçek yönler yalnızca çizimde işaretlenir; alıcıya geri verilmez.
    a = s.alici;
    assert(a.hazir,'Uzamsal kestirim henüz hazır değil.');
    az = a.music{1}.AzimuthScanAngles;
    el = a.music{1}.ElevationScanAngles;
    if isfield(a,'sonUzamsalMs') && isfinite(a.sonUzamsalMs)
        t = a.sonUzamsalMs;
    else
        t = a.edinmeMs + floor((a.epoch-a.edinmeMs)/a.guncellemeMs)*a.guncellemeMs;
    end
    if isfield(s,'gercekUydu')
        uydu = s.gercekUydu;
    else
        % toolbox-v1 eski kayıtlarında saklanmayan sabit senaryo açıları.
        uydu = [30 120;50 70];
    end
    spf = s.gercekSpf(t,1:2);

    %% Aday kapılarının MUSIC haritaları: üç kaynağı tek MUSIC'e vermiyoruz.
    harita = zeros(numel(el),numel(az));
    for i = 1:numel(a.kanallar)
        if isfield(a,'musicSpektrum') && numel(a.musicSpektrum)>=i && ...
                ~isempty(a.musicSpektrum{i})
            p = a.musicSpektrum{i};
        else
            % Eski .mat dosyaları da yeniden deney koşmadan çizilebilir.
            g = a.kanallar{i}.grup;
            p = a.music{g}(a.C{i});
        end
        p = el_satirlari(p,az,el);
        harita = harita + p/max(max(p(:)),eps);
    end
    db = 10*log10(max(harita/max(max(harita(:)),eps),eps));
    fMusic = figure('Name','Son Durum - MUSIC 2D','Position',[100 100 1200 650]);
    ax = axes(fMusic);
    imagesc(ax,az,el,db); axis(ax,'xy'); hold(ax,'on');
    colormap(ax,parula); clim(ax,[-15 0]);
    cb = colorbar(ax); cb.Label.String = 'Normalize MUSIC spektrumu [dB]';
    h1 = plot(ax,uydu(1,1),uydu(2,1),'ko','MarkerFaceColor','w', ...
        'MarkerSize',13,'LineWidth',2.5);
    h2 = plot(ax,uydu(1,2),uydu(2,2),'co','MarkerSize',13,'LineWidth',2.5);
    hs = plot(ax,spf(1),spf(2),'ms','MarkerSize',13,'LineWidth',2.5);
    hx = plot(ax,a.yonler(1,:),a.yonler(2,:),'rx','MarkerSize',13,'LineWidth',2.5);
    xlabel(ax,'Azimut [derece]'); ylabel(ax,'Elevasyon [derece]');
    title(ax,{sprintf('Son MUSIC kestirimleri | %d ms',t), ...
        'Her aday spektrumu kendi tepesine normalize edilerek birleştirildi'});
    legend(ax,[h1 h2 hs hx],{'PRN1 gerçek uydu','PRN2 gerçek uydu', ...
        'Gerçek spoofer','MUSIC kestirimleri'},'Location','southwest');
    xlim(ax,[-180 180]); ylim(ax,[0 90]);

    %% Toolbox array pattern: son ağırlıklarla |w'' a(az,el)|^2 yüzeyleri.
    fMvdr = figure('Name','Son Durum - MVDR 3D','Position',[80 80 1450 700]);
    tiles = tiledlayout(fMvdr,1,2,'TileSpacing','compact','Padding','compact');
    for g = 1:2
        pat = pattern(a.dizi,a.fL1,az,el,'Weights',a.w(:,g), ...
            'PropagationSpeed',physconst('LightSpeed'),'Type','powerdb','Normalize',true);
        pat = el_satirlari(pat,az,el);
        pat = max(pat,-50); % Sadece çizim tabanı; gerçek sönüm ölçümünü değiştirmez.
        ax = nexttile(tiles);
        surf(ax,az,el,pat,'EdgeColor','none'); hold(ax,'on');
        zUydu = interp2(az,el,pat,uydu(1,g),uydu(2,g));
        zSpf = interp2(az,el,pat,spf(1),spf(2));
        hd = plot3(ax,uydu(1,g),uydu(2,g),zUydu,'ko', ...
            'MarkerFaceColor','w','MarkerSize',10,'LineWidth',2.5);
        hs = plot3(ax,spf(1),spf(2),zSpf,'ms','MarkerSize',10,'LineWidth',2.5);
        colormap(ax,parula); clim(ax,[-50 0]);
        cb = colorbar(ax); cb.Label.String = 'Normalize güç [dB]';
        xlabel(ax,'Azimut [derece]'); ylabel(ax,'Elevasyon [derece]');
        zlabel(ax,'Normalize güç [dB]');
        title(ax,sprintf('MVDR%d - PRN%d uydu hüzmesi',g,g));
        legend(ax,[hd hs],{'Gerçek hedef uydu','Gerçek spoofer'},'Location','southoutside');
        xlim(ax,[-180 180]); ylim(ax,[0 90]); zlim(ax,[-50 0]);
        xticks(ax,-180:90:180); yticks(ax,0:30:90);
        grid(ax,'on'); view(ax,[-35 35]);
    end
    title(tiles,sprintf('Son MVDR ağırlıkları | %d ms | Spoofer: %.1f / %.1f derece', ...
        t,spf(1),spf(2)));
end

function p = el_satirlari(p,az,el)
% MUSIC/pattern çıktısını imagesc ve surf için elevasyon x azimut düzenine getir.
    if isequal(size(p),[numel(az),numel(el)]), p = p.'; end
    assert(isequal(size(p),[numel(el),numel(az)]),'Beklenmeyen uzamsal spektrum boyutu.');
end
