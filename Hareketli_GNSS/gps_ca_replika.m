function kod = gps_ca_replika(prn, fs)
%GPS_CA_REPLIKA Alıcının bildiği, veri bitlerinden bağımsız 1 ms C/A kodu.
% Bu fonksiyon uydu/spoofer dalga biçimi, Doppler, gecikme veya navdata almaz.
% Edinme ve takip, bu kodun gecikmesini/taşıyıcısını X üzerinden kestirmelidir.
    validateattributes(prn, {'numeric'}, {'scalar','integer','>=',1,'<=',210});
    validateattributes(fs, {'numeric'}, {'scalar','positive','finite'});
    ornekChip = fs / 1.023e6;
    assert(abs(ornekChip - round(ornekChip)) < 1e-10, ...
        'Bu sürümde örnekleme hızı 1.023 MHz''in tam katı olmalıdır.');
    chipler = 1 - 2*double(gnssCACode(prn, 'GPS'));
    kod = repelem(chipler(:), round(ornekChip));
end
