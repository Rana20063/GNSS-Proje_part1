function [y, durum, wSon] = gps_uzamsal_nlms(X,wMVDR,korunanVektorler,durum,d)
%GPS_UZAMSAL_NLMS Toolbox NLMS ile hedefi koruyan uzamsal girişim giderme.
% B'nin her sütunu bir yardımcı anten hüzmesidir; hedef yönü içermez.
% Her yardımcı hüzme için dsp.LMSFilter(1,Method="Normalized LMS") kullanılır.
% Bir katsayı zamansal eşitleme değil, o yardımcı hüzmenin uzamsal katsayısıdır.
    if nargin<5, d=X*conj(wMVDR); end
    B = null(korunanVektorler');
    if isempty(durum) || ~isequal(size(durum.B),size(B)) || norm(durum.B-B,'fro')>1e-8
        durum = struct('B',B,'filtre',{{}});
        for k=1:size(B,2)
            durum.filtre{k}=dsp.LMSFilter(1,'Method','Normalized LMS', ...
                'StepSize',0.002,'LeakageFactor',0.99999);
        end
    end
    y=d; wSon=wMVDR;
    for k=1:size(B,2)
        referans=X*conj(B(:,k));
        % yHat girişim kestirimi; hata=ana hüzme-yHat temizleme çıkışıdır.
        [~,y,katsayi]=durum.filtre{k}(referans,y);
        wSon=wSon-B(:,k)*conj(katsayi);
    end
    % wSon son örneğin katsayılarıdır; bütün bloğun sabit ağırlığı değildir.
    % NLMS başarısı hüzme diyagramıyla değil çıkış korelasyonlarıyla ölçülür.
end
