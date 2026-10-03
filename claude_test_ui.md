Testasin CellShieldin selaimessa: kaikki kolme presetiä, ääriarvot (w 0.5–1.0 mm, R_in 0.2–2.0 mm, i 2.5 A/cm², V 0.5 V), kaikki välilehdet ja validointiosion. Allsolve-nappia en painanut, koska se lähettää työn ulkoiseen palveluun. Sovellus on nyt jäänyt viimeisiin testiarvoihin, ja sivun päivitys palauttaa sen.

Mikä toimii
Kaikki toimii virheittä: laskenta kestää 7–15 ms, eikä konsolissa ole virheitä.
Tarkistin luvut käsin, ja ne täsmäävät:
Re ≈ 120 ja Dean ≈ 66.
Koko kanavan painehäviö 118 Pa × 40 mutkaa = 4,73 kPa.
Kompressoriteho 83 W = 4,73 kPa × 0,0105 m³/s / 0,6.
f·Re = 62,23 (tarkka arvo 62,19).
Arviointiteksti päivittyy tilanteen mukaan, ja Pareto-suositus ("same pressure budget → 6.2 %") on hyödyllinen.
Löydökset
Esiasetus "High Load (Critical Flooding & Heat)" näyttää kaiken vihreänä: 77,6 °C ja tulvimisriski 7,9 % "LOW". Low Load antaa tulvimisriskiksi 24 % "ELEVATED", eli tulokset ovat päinvastaiset kuin esiasetusten nimet lupaavat.
Tulvimismittari ei riipu kuormasta. Kun nostin virrantiheyden 2,0 → 2,5 A/cm², tulvimisriski pysyi 6,9 %:ssa. Malli katsoo vain, missä nopeus on alle 1 m/s, eikä huomioi syntyvää vettä, kosteutta tai lämpötilaa.
Stoikiometria puuttuu. Inlet-nopeus ja virrantiheys ovat toisistaan riippumattomia. Karkea arvioni (oletuksina noin 1,5 bar, 80 °C ja pinta-ala = legit × jakoväli) on, että High Load -presetin λ on noin 0,2–0,4, eli kenno kärsisi hapen puutteesta. Arvio kannattaa tarkistaa koodissa.
Geometria vaikuttaa tuskin lainkaan:
Kun R_in kasvaa 0,6 → 2,0 mm, tulvimisriski muuttuu vain 23,7 → 23,6 %, mutta painehäviö kasvaa 40 → 48 Pa.
Pareto-kuvassa pisteet ryhmittyvät lähes pelkästään inlet-nopeuden mukaan.
Mutkahäviö on vain 0,4–0,7 % painehäviöstä. Turbulenttiin virtaukseen tarkoitettu Idelchikin kerroin todennäköisesti aliarvioi häviön, kun Re on 40–120.
Ribin leveys ei vaikuta lämpötilaan. Ribin lämpötilanousu on 7,15 K, kun ribin leveys on 0,4 mm, ja täsmälleen sama, kun se on 4 mm.
Painehäviökortti on epäjohdonmukainen. Kortti näyttää arvon per U-mutka, mutta sen väri tulee koko kanavan arvosta (5,58 kPa > 5 kPa). Kortista puuttuu myös status-merkki, joka muissa korteissa on.
Tulvimisprosentti voi riippua laskentahilasta. Low Loadilla kenttä näyttää lähes tasaiselta (1,4–1,5 m/s), mutta silti 24 % merkitään riskialueeksi. Todennäköisesti seinän vierussolut lasketaan mukaan.
Kuvaajien ongelmat:
Otsikko lupaa punaisen 1 m/s -isoviivan, mutta viiva piirtyy oranssina ja sulautuu seinäviivoihin.
Legin ja mutkan liitoskohdassa näkyy sauma.
Plotlyn työkalupalkki peittää otsikon.
GDL-poikkileikkauksen otsikko katkeaa.
Membraanikartta on pieni: akseli on ±20 mm, vaikka geometria on noin ±5 mm.
Lämpötilakartoilla on kaksi eri väriskaalaa.
Kortin "rib 72,2 °C" ja poikkileikkauksen noin 68 °C rib-pinnassa ovat ristiriidassa.
Jännite ja virrantiheys ovat toisistaan irrallaan. Esimerkiksi i = 2,5 A/cm² yhdessä V = 0,9 V:n kanssa on sallittu yhdistelmä. Esiasetus ei myöskään palauta geometriaa, eikä muokattua operointipistettä merkitä mitenkään ("custom").
Prompti parannuksia varten
Paranna CellShield-sovellusta (app.py + physics.py, Streamlit). Säilytä nykyinen arkkitehtuuri (Design-dataclass → evaluate() → Plotly) ja alle 50 ms:n laskenta-aika. UI-tekstit pysyvät englanniksi. Toteuta prioriteettijärjestyksessä ja lisää jokaiselle kohdalle pytest-testi.

P1 – FYSIIKAN JOHDONMUKAISUUS
1. Stoikiometria: laske katodin λ = (ilman O2-moolivirta) / (i·A_active/(4F)), jossa A_active = kanavan pituus × jakoväli (w + w_rib). Lisää paineelle ja lämpötilalle syötteet (oletus 1.5 bar, 80 °C). Näytä λ KPI-korttina: varoitus, kun λ < 1.5; kriittinen, kun λ < 1.0. Lisää vaihtoehtoinen tila, jossa käyttäjä antaa λ:n ja v_in lasketaan siitä.
2. Tulvimismalli: korvaa pelkkä |u| < 1 m/s -kynnys yhdistelmäindeksillä, joka huomioi (a) tuotetun veden i/(2F), (b) vesihöyryn kyllästymisen ulostulossa (RH_in-syöte, p_sat(T)), josta saadaan nestemäisen veden osuus, ja (c) pisaroiden irtoamiskriteerin (esim. leikkausjännitys vs. pintajännitys/kontaktikulma) matalan nopeuden alueilla. Tulvimisriskin täytyy kasvaa virrantiheyden kasvaessa, kun nopeus pysyy samana.
3. Polarisaatiokäyrä: laske V_cell(i) yksinkertaisella mallilla (OCV, Tafel, ohminen, konsentraatio) tai lisää vähintään varoitus epäfysikaalisille (i, V)-pareille. Lämpövuo q = i·(E_tn − V) käyttää tätä jännitettä.
4. Ribin leveyden vaikutus lämpöön: lisää ribin poikittainen johtuminen ja kontaktiresistanssi, jotta w_rib = 2·R_in vaikuttaa membraanin huippulämpötilaan. Nyt ribin lämpötilanousu on sama 7,15 K, kun w_rib on 0,4 mm tai 4 mm.
5. Laminaarin mutkan häviö: käytä Re-riippuvaa K:ta (esim. K = K_turb + C/Re tai Dean-luvun korjaus) Idelchikin turbulentin kertoimen sijaan, kun Re < 2000.
6. Hilariippumattomuus: tarkista, että tulvimisprosentti ei muutu yli 1 prosenttiyksikköä, kun hila tihennetään 2×. Jätä seinän vierussolut pois tai painota ne pinta-alalla. Lisää tarkistus validointitaulukkoon.

P2 – KALIBROINTI JA ESIASETUKSET
7. Viritä esiasetukset niin, että "High Load (Critical Flooding & Heat)" todella näyttää vähintään yhden kriittisen KPI:n, Nominal näyttää vihreää ja Low Load näyttää tulvimisriskin matalan λ:n tai nopeuden vuoksi. Merkitse "(custom)", kun käyttäjä muuttaa presetin arvoja, ja lisää "Reset to preset" -nappi.
8. Painehäviökortti: saman suureen pitää ohjata arvoa ja väriä. Näytä koko kanavan painehäviö pääarvona (raja 5/10 kPa) tai skaalaa raja per mutka. Lisää status-merkki kuten muihin kortteihin.
9. Engineering assessment: listaa kaikki rajan ylittävät KPI:t (nyt näkyy vain yksi). Suosittele vain toimenpiteitä, joiden vaikutus on mallissa todennettu: laske suositellun muutoksen vaikutus ennen sen näyttämistä (esim. "R_in 0.6 → 1.2 mm: flooding −X pp, ΔP +Y Pa").

P3 – PARETO
10. Erota suunnittelumuuttujat (w, R_in) operointimuuttujista (v_in/λ). Lisää 3-tavoitteinen tarkastelu (tulvimisriski, painehäviö, T_max) ja w sweepiin. Näytä varoitusraja 15 % sekä kriittinen raja ja korosta Pareto-rintama selvemmin.

P4 – VISUALISOINNIT
11. Piirrä 1 m/s -isoviiva oikeasti punaisena ja seinät eri värillä (esim. valkoinen ohut viiva). Poista sauma legin ja mutkan liitoksesta (sileä siirtymä). Siirrä Plotlyn modebar pois otsikon päältä.
12. Lämpötilakartat: yhteinen väriskaala molemmille kuville. Korjaa GDL-poikkileikkauksen otsikon katkeaminen ja käytä yhtenäistä akselisuhdetta (scaleanchor). Rajaa membraanikartan akselit geometriaan. Selvennä, mitä kortin "rib xx °C" tarkoittaa (rib-pinta vs. membraani ribin alla), ja varmista, että arvo vastaa poikkileikkauksen kenttää.
13. Velocity-kartan akselirajat skaalautumaan geometrian mukaan (nyt kiinteä ±4 mm).

P5 – ALLSOLVE
14. Näytä ennen lähetystä vahvistusikkuna (arvioitu ajoaika ja verkko). Näytä jobin tila ja palauta 3D-tulokset rinnakkain surrogaatin kanssa (T_max, painehäviö, tulvimisalue, poikkeama %). Tallenna verrokit, jotta surrogaattia voi myöhemmin kalibroida.

Hyväksymiskriteerit:
- High Load -preset näyttää vähintään yhden punaisen tai keltaisen KPI:n.
- Tulvimisriski kasvaa virrantiheyden kasvaessa, kun v_in pysyy samana.
- T_max muuttuu ribin leveyden mukana.
- Painehäviökortin väri ja arvo vastaavat samaa suuretta.
- λ on näkyvissä.
- Kaikki testit menevät läpi.
- evaluate() kestää alle 50 ms.