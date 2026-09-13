# GNSS Anti-Spoofing Projesi — Aşama 1 Teknik Dokümanı

## 1. Proje Özeti

Proje, **GPS L1 C/A sinyallerine yönelik akıllı aldatma (spoofing) saldırılarını tespit etmek** amacıyla geliştirilmiş MATLAB/Simulink tabanlı arayüze sahip bir araç oluşturmayı hedefliyor şeklinde anladım.

**Ana Hedef:** Tek anten, tek otantik uydu (Yani gerçek GPS uydusu) sinyali ve tek spoofer sinyali kombinasyonunda iki ayrı sinyalin varlığını korelasyon analizi ile tespit etmek; bu tespitin üzerine sonraki aşamalarda (CRPA, S-curve, ARAIM/FDE) katmanlar eklemek.



---

## 2. Projenin Genel Adımları 

Projeyi bir bütün olarak dört aşamada düşünebiliriz. Her aşama bir öncekinin üzerine bir "katman" ekliyor — yani önce sahte bir GPS ortamı kuruyoruz, sonra bu ortamda saldırıyı fiziksel olarak bastırmaya çalışıyoruz, sonra sinyal seviyesinde saldırıyı tespit ediyoruz, en son da konum hesabında hâlâ sızan hataları temizliyoruz. Biz şu an sadece **1. aşamadayız**.

**Aşama 1 — Sahte Bir Uydu Ortamı Kurmak**
Önce gerçek bir GPS sinyalinin nasıl göründüğünü bilgisayarda üretiyoruz: uydunun kod imzası (PRN), gerçek bir alıcının göreceği frekans kayması, gecikme gibi bozulmalar. Sonra bu temiz ortama bir de sahte cihaz ekliyoruz: zamanla gücünü artırıp alıcıyı kandırmaya çalışan bir sahte sinyal (spoofer) ve doğal bir yansıma (multipath). Yani bu aşamanın amacı, sonraki aşamaları test edebileceğimiz gerçekçi bir modelleme ortamı hazırlamak.

**Aşama 2 — Anten ile Saldırıyı Fiziksel Olarak Bastırmak (CRPA)**
Tek bir anten yerine 4 antenden oluşan bir dizi (CRPA) kullanıyoruz. Bu antenler birlikte çalışarak, sahte sinyalin geldiği yöne doğru bir tür sağırlık ? oluşturabiliyor — yani spoofer'ın sesini kısıp gerçek uydunun sesini olduğu gibi bırakıyor. Buna *null-steering* deniyor. Bu, saldırıyı sinyal daha işlenmeden, donanım seviyesinde zayıflatmanın yolu.

**Aşama 3 — Sinyal Seviyesinde Doğrulama Yapmak**
Anten katmanı saldırıyı tamamen susturamayabilir; hâlâ sinyalde iz kalabilir. Bu yüzden alıcının izleme devrelerinden (Early-Prompt-Late gibi) gelen verilerle bir korelasyon haritası çıkarıyoruz. Bu harita bize "burada tek bir sinyal mi var, yoksa üst üste binmiş birden fazla sinyal mi var?" sorusunun cevabını veriyor.

**Aşama 4 — Konum Hesabını Güvenceye Almak (ARAIM/FDE)**
En küçük ve fark edilmesi zor sapmalar bile buraya kadar sızabilir. Son katmanda, konum hesaplamasına giren tüm uydu ölçümlerini birbiriyle karşılaştırıp, tutarsız/bozulmuş olanı tespit edip devre dışı bırakan algoritmalar (ARAIM, FDE) çalışıyor. Yani en sonunda "hangi ölçüme güvenebiliriz?" sorusuna cevap bulunuyor ve nihai konum bu temiz ölçümlerle hesaplanıyor.

**Özetle:** RF ortamı kur → anten ile bastır → sinyalde doğrula → konumda son kontrolü yap. Bu doküman, sadece ilk aşamada geliştirilen beş bileşeni anlatıyor.

Bu Repo'da 1. Aşama için bir simülasyon oluşturdum:

---

## 3. Sistemin Beş Temel Bileşeni ve Dosyalar
## Proje Dosyaları


- **`build_simulink_model.m`**
  `adim1_basit_antispoofing.slx` modelini arayüze hiç dokunmadan sıfırdan ve otonom olarak inşa eden "wrapper" (sarmalayıcı) scripttir. Taşınabilirlik sorunlarını çözer ve sistem parametrelerini (IF, Doppler, drift hızı) kurgular.

- **`correlation_detector.m`**
  Simulink'ten bağımsız, saf MATLAB tabanlı test simülasyonunda korelasyon sinyalindeki tepeleri (peaks) tespit eden ve spoofing kararını veren temel algoritma fonksiyonudur.

- **`generate_ca_code.m`**
  GPS ICD-200 standartlarına uygun olarak PRN uyduları (örn. PRN 7) için 1023 çiplik Gold kodunu üreten ve sinyalin en temel yapısını (baz bant) sentezleyen jeneratördür.

- **`generate_if_signal.m`**
  Üretilen baz bant PRN çiplerini alıp, ara frekansa (IF) ve Doppler kaymasına sahip fiziksel taşıyıcı dalgalara (RF ortamına) dönüştüren sinyal modülasyon dosyasıdır.

- **`main_basic_simulation.m`**
  Tüm sistemi Simulink'e taşımadan önce saf MATLAB koduyla hızlıca test etmek, parametreleri doğrulamak ve korelasyon ile tepe tespit yeteneklerini görmek için kullanılan ana başlatıcı test betiğidir. İlk önce kontrol amaçlı çalıştırılması önerilir.

- **`run_simulink_correlation.m`**
  Simulink çalışırken her 10.000 veride bir çağrılan köprü fonksiyondur. Gelen IF sinyaline donanım seviyesinde I/Q demodülasyonu uygulayarak taşıyıcıyı yok eder (Envelope Detector) ve ortamdaki sahte/gerçek sinyal sayısını hesaplayıp Simulink'e geri gönderir.

### 1) GPS L1 C/A PRN Kod Üretici

**Dosya:** `generate_ca_code.m`

**Ne yapar?**
Her GPS uydusuna özgü, 1023 chip uzunluğunda bir Gold kodu üretir. Bu kod, uydunun "imzası" gibidir — alıcı bu kodu bilerek doğru uyduyu diğerlerinden ayırt edebilir.

**Neden 1023 chip, neden Gold kodu?**
- GPS standardı (ICD-200), C/A kodunu 1.023 MHz chip hızında, 1 milisaniyelik bir periyotta yayınlar (1023 = 1.023 MHz × 1 ms).
- Gold kodları, iki ayrı shifted register (G1 ve G2) çıkışının XOR'lanmasıyla üretilir ve çok iyi bir **oto-korelasyon** özelliğine sahiptir — yani kod kendisiyle çakıştığında keskin bir tepe verir, kaymış hâldeyken ise neredeyse sıfıra yakın bir değer üretir. Bu özellik, spoofing tespitinin matematiksel temelini oluşturur.

**Nasıl çalışır?**
```matlab
% G1: Polinom x^10 + x^3 + 1
g1_fb = mod(g1(3) + g1(10), 2);  % Tap 3 ve 10

% G2: Polinom x^10 + x^9 + x^8 + x^6 + x^3 + x^2 + 1
g2_fb = mod(sum(g2([2 3 6 8 9 10])), 2);  % 6 tap
```
G1 ve G2 kaydırıcıları GPS ICD-200 Tablo 3-I'deki spesifikasyona birebir uyacak şekilde tasarlanmıştır. PRN 7 için üretilen 1023-chiplik vektör, gerçek GPS uydu 7 sinyaliyle **bitwise özdeştir** — yani doğrulanabilir, standarda uygun bir referans koddur.

---

### 2)  IF Modülasyon Motoru

**Dosya:** `generate_if_signal.m`, `build_simulink_model.m` içindeki sinyal üretici scriptler

**Ne yapar?**
Üretilen baseband (taban bant) sinyalini, gerçek bir RF alıcısının göreceği fiziksel ortama taşır: 1.25 MHz ara frekans (IF), 1500 Hz Doppler kayması ve 0.5 chip atmosferik (iyonosferik) gecikme ekleyerek.

**Neden gerekli?**
Baseband benzetimi matematiksel olarak temizdir ama gerçekçi değildir. Gerçek bir GPS alıcısı sinyali doğrudan taban bantta almaz; önce bir taşıyıcı dalga (IF) üzerine binmiş hâlde, Doppler kaymasıyla bozulmuş ve atmosferik gecikmeyle kaymış olarak alır. Bu modül, sistemi bu gerçekçi koşullara maruz bırakarak sonraki bileşenlerin (özellikle I/Q demodülasyonun) neden gerekli olduğunu da ortaya koyar.

**Nasıl çalışır?**
```matlab
% Otantik: y_auth(t) = ca(t) × cos(2π×(f_IF + f_doppler)×t)
% Spoofer: y_spoof(t) = 1.4×ca(t) × cos(2π×f_IF×t)
% Multipath: y_mp(t) = 0.4×ca(t) × cos(2π×(f_IF+f_doppler)×t)
```
Sentezlenen baz bant sinyalleri, 1.25 MHz'lik kosinüs taşıyıcı dalgasına bindirilerek fiziksel RF donanım gerçekliğine taşınır. Bu adım olmadan, sistem sadece "kağıt üzerinde" çalışan idealize bir model olarak kalır.

**Parametreler:**
| Parametre | Değer | Anlamı |
|---|---|---|
| f_IF | 1.25 MHz | Tipik GPS alıcı ara frekansı |
| f_doppler | 1500 Hz | Platform hareketini simüle eder |
| atmos_delay | 0.5 chip | İyonosferik gecikme |

---

### 3) Dinamik Spoofer Modeli

**Dosya:** `build_simulink_model.m` içindeki `spoofScript`

**Ne yapar?**
Gerçekçi, modern bir GPS aldatma saldırısını modeller. Saldırı, otantik sinyalle **tam aynı fazda ve aynı güçte** başlar — böylece alıcı tarafından fark edilmez. Zamanla saldırı sinyali hem faz olarak kayar (sürüklenir) hem de gücünü artırır.

**Neden bu şekilde tasarlanmış?**
Gerçek spoofing saldırılarının klasik stratejisi budur: alıcının DLL (Delay Locked Loop) döngüsü otantik sinyale kilitliyken, spoofer önce o kilitle eşleşir, sonra alıcıyı yavaş yavaş kendi (sahte) fazına doğru çeker. Ani bir sıçrama alıcının alarm vermesine yol açar; kademeli sürüklenme ise tespit edilmeden konum bilgisini manipüle etmeyi mümkün kılar.

**Nasıl çalışır?**
```matlab
persistent drift genlik

% Her adımda:
drift = drift + 0.002;                 % Kod fazını kademeli kaydır
genlik = min(genlik + 0.0001, 2.0);    % Gücü kademeli artır (max 2×)

idx = mod(floor(phase + atmos_delay + drift), 1023) + 1;
baseband = ca(idx);
y = genlik × baseband × cos(...);
```
- **Drift:** Spoofer, alıcının izleme döngüsünü takip edip onu otantik sinyalden uzaklaştırmaya çalışır (kod fazı sürüklemesi).
- **Genlik artışı:** Başlangıçta gizli kalmak için otantik güçte başlar; alıcı sahte sinyali takip etmeye başladıkça gücünü artırarak baskın hâle gelir.

Bu davranış, `drift = 0` yapılarak devre dışı bırakılabilir ve statik (sürüklenmeyen) bir spoofer senaryosu test edilebilir.

---

### 4 ) Multipath (Yankı) Simülatörü

**Dosya:** `build_simulink_model.m` içindeki `mpScript`

**Ne yapar?**
Otantik sinyalin, çevredeki yüzeylerden (bina, yer, su vb.) yansıyarak alıcıya gecikmeli ve zayıflamış şekilde ikinci kez ulaşmasını modeller. Bu, gerçek dünyada her GPS alıcısının karşılaştığı doğal bir bozulma kaynağıdır.

**Neden gerekli?**
Multipath, spoofing ile karıştırılabilecek bir "sahte tepe" kaynağıdır. Bir anti-spoofing sisteminin güvenilir olması için, gerçek bir saldırıyı sıradan bir yansımadan ayırt edebilmesi gerekir. Bu blok, sistemin bu ayrımı yapabildiğini test etmek için kasıtlı olarak eklenmiştir.

**Nasıl çalışır?**
Otantik sinyalin **0.4 genlikli** (yani zayıflatılmış), **302. chip fazına** kaymış (yaklaşık 2 chip gecikmeli) bir kopyası ayrı bir blokta üretilip anten toplama noktasına eklenir. Bu genlik, spoofer sinyaline göre (genellikle 1.0 ve üzeri) belirgin şekilde daha düşük tutulur; böylece doğru kalibre edilmiş bir sistem, multipath tepesini spoofing tepesinden ayırt edebilir.

---

### 5) I/Q Demodülasyon ve Zarf Dedektörü

**Dosya:** `run_simulink_correlation.m`

**Ne yapar?**
IF taşıyıcı dalgasının korelasyon sonucuna kattığı yapay bozulmaları (yüksek frekanslı salınımları) matematiksel olarak yok eder ve geriye yalnızca sinyalin gerçek genlik profilini ("zarf") bırakır.

**Neden gerekli? (Çözdüğü problem)**
Sinyal IF taşıyıcısına bindirildikten sonra, doğrudan korelasyon alınırsa (naif yöntem) sonuç şu şekilde bozulur:
```
corr(t) = rxSignal(t) × ca(t)
        = 1×cos(2π×f_IF×t)×ca²(t)      [otantik tepe]
        + 1.4×cos(2π×f_IF×t)×ca²(t)    [spoofer tepe]
        + [yüksek frekans kalıntıları]  ← SAHTE TEPE!
```
1.25 MHz'lik taşıyıcı, korelasyon profilinde saniyede milyonlarca kez tekrarlanan bir "testere dişi" desen oluşturur. Bu durumda `findpeaks()` fonksiyonu gerçek sinyal tepelerini değil, bu yapay salınımları da tepe olarak sayar (örneğin 2 yerine 7-10 tepe tespit edilir).

**Çözüm nasıl çalışır?**
```matlab
% In-phase (I) kanal:
local_sig_I = baseband × cos(2π×(f_IF + doppler)×t)
corr_I = xcorr(rxSig, local_sig_I)

% Quadrature (Q) kanal:
local_sig_Q = baseband × sin(2π×(f_IF + doppler)×t)
corr_Q = xcorr(rxSig, local_sig_Q)

% Zarf (Envelope):
corrResult = sqrt(corr_I² + corr_Q²)
```
Matematiksel olarak:
```
I² + Q² = [baseband×cos(θ)]² + [baseband×sin(θ)]²
        = baseband² × [cos²(θ) + sin²(θ)]
        = baseband²   ← taşıyıcı terimi (θ) tamamen kayboldu
```
`cos²(θ) + sin²(θ) = 1` trigonometrik özdeşliği sayesinde taşıyıcı fazı denklemden çıkar ve geriye yalnızca sinyalin gerçek genlik bilgisi kalır. Sonuç: pürüzsüz, kalıntısız bir korelasyon profili — ve `findpeaks()` artık doğru sayıda tepeyi (gerçek sinyal kaynağı sayısını) tespit eder.

---

## 4. Bileşenlerin Bir Arada Çalışması

Bu beş bileşen, işlem hattında şu sırayla birleşir:

1. **PRN Kod Üretici**, referans Gold kodunu üretir.
2. **IF Modülasyon Motoru**, bu kodu otantik sinyal, spoofer sinyali ve multipath yankısı olarak üç ayrı fiziksel bileşene dönüştürüp taşıyıcı üzerine bindirir.
3. **Dinamik Spoofer Modeli**, spoofer bileşenini zamanla sürükleyip güçlendirerek gerçekçi bir saldırı senaryosu oluşturur.
4. **Multipath Simülatörü**, ortama doğal bir yankı bileşeni ekler.
5. Anten toplama noktasında bu üç sinyal (otantik + spoofer + multipath) toplanır ve gürültüyle karışır.
6. **I/Q Demodülasyon ve Zarf Dedektörü**, bu karışık sinyali taşıyıcıdan arındırıp temiz bir korelasyon profiline dönüştürür.
7. Korelasyon tabanlı tepe dedektörü (`findpeaks`), bu temiz profildeki tepe sayısını sayar: **1 tepe = normal, 2+ tepe = spoofing/multipath şüphesi**.

## Simulink  Scope Output yorumlanması
Simülasyon çalışırken tespit edilen korelasyon tepe sayısının zaman içindeki değişimi, senaryonun dinamik yapısını matematiksel olarak doğrular. Grafiğin 1'e hiç uğramadan basamaklar halinde 0'dan 2'ye, ardından 3'e çıkması şu şekilde açıklanır:

* **0 Seviyesi:** Simülasyonun ilk anlarıdır. Henüz veri tamponu dolmadığı veya sistem kilitlenmediği için tespit yapılamaz.
* **2 Seviyesi (Neden 1 Değil?):** Ortamda en başından beri kasıtlı bir **Multipath** (yankı) sinyali bulunmaktadır. Başlangıçta **Spoofer** sinyali, **Otantik** sinyal ile tamamen aynı fazda olduğu için bu ikisi üst üste binerek tek bir tepe gibi algılanır. Yankı (Multipath) sinyali ise ayrı bir tepe oluşturur. Bu nedenle sistem ilk ölçümde 1 değil, 2 tepe (Otantik+Spoofer ve Multipath) tespit eder.
* **3 Seviyesi (Saldırının Belirginleşmesi):** İlerleyen zamanlarda Dinamik Spoofer modeli, alıcıyı kendi üzerine çekmek için fazını kaydırmaya (drift) başlar. Fazlar ayrıştığı anda, üst üste binen tepe ikiye bölünür. Sistem artık **Otantik**, **kaymış Spoofer** ve **Multipath** olmak üzere birbirinden bağımsız 3 farklı sinyali tespit ederek saldırıyı tamamen açığa çıkarır .

---

## 5. Teknik Terimler Açıklama

| Terim | Anlam |
|---|---|
| **Spoofing** | GPS sinyalinin yapay taklit saldırısı |
| **C/A Kod** | Coarse/Acquisition — 1023 chiplik, 1.023 MHz'de, 1 ms periyotlu açık GPS kodu |
| **PRN** | Pseudo-Random Noise — her uyduya özgü kod üretim anahtarı (1-32) |
| **IF** | Intermediate Frequency — ara frekans (GPS alıcılarında tipik 1-2 MHz) |
| **Doppler** | Hareketin neden olduğu frekans kayması |
| **Multipath** | Sinyalin yansıması sonucu gecikmeli ikinci kez alınması |
| **I/Q Demodülasyon** | Taşıyıcı dalgasını matematiksel olarak kaldırma yöntemi |
| **Zarf (Envelope)** | Taşıyıcısız, sinyalin gerçek genlik profili — sqrt(I² + Q²) |
