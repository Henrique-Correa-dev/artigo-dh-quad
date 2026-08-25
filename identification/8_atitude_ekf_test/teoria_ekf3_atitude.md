# Como o Pixhawk estima a atitude — Análise teórica (ArduPilot EKF3)

> Fase 1 (teoria) do estudo em `8_atitude_ekf_test/`.
> Objetivo: entender **exatamente** de onde vem o `Roll/Pitch/Yaw` do `ATT` no teu
> `logs_concat.mat`, quais medidas e algoritmos o firmware usa, e quais "pesos"
> ele aplica — para depois (fase 2) **testar essas hipóteses** contra os dados.

---

## 0. TL;DR — a configuração REAL do teu drone

Extraído do `PARM` e do `VER` do log (`5 ...bin-155264.mat`):

| Item | Valor no teu log | O que significa |
|---|---|---|
| Firmware | ArduPilot **4.6.3** (QuadPlane VTOL) | `VER`: Maj=4, Min=6, Pat=3 |
| `AHRS_EKF_TYPE` | **3** | A atitude (`ATT`) vem do **EKF3** |
| `EK3_ENABLE` | 1 (EK2 ausente) | Só o EKF3 roda |
| `EK3_IMU_MASK` | **3** (binário `11`) | **2 cores** rodando em paralelo: IMU0 e IMU1 → `XKF1_0` e `XKF1_1` |
| `EK3_SRC1_POSXY/VELXY` | 3 / 3 | Posição e velocidade horizontal → **GPS** |
| `EK3_SRC1_POSZ` | 1 | Altitude → **Barômetro** |
| `EK3_SRC1_VELZ` | 3 | Velocidade vertical → **GPS** |
| `EK3_SRC1_YAW` | 1 (Compass) | Yaw vem do **magnetômetro** |
| `COMPASS_USE` / `COMPASS_USE2` | **0 / 1** | Compass **1 desligado**, mas **compass 2 (externo) LIGADO e usado** ✔️ |

> ✔️ **Verificado nos dados (correção de uma hipótese inicial errada):** o drone tem
> **dois magnetômetros**. O compass 1 está desligado (`COMPASS_USE=0`), mas o
> **compass 2 (externo, `COMPASS_USE2=1`) está LIGADO e é fundido pelo EKF** —
> confirmado pelas inovações de mag ativas em `XKF3` (`IMX/IMY/IMZ`, std 40–66 mGauss)
> e pelo test ratio `SM` (0–1.2) em `XKF4`. Ou seja: **o yaw vem do magnetômetro
> (compass 2)**, como manda `EK3_SRC1_YAW=1`. O estimador **EKF-GSF** (`XKY0/XKY1`)
> roda **em paralelo** como *backup* e monitor de consistência, mas não é a fonte
> primária aqui. As inovações de mag são **grandes/ruidosas**, o que ajuda a explicar
> yaw ser o canal mais "mole" (ver `project_quad_sysid`): roll/pitch são amarrados com
> força pela velocidade do GPS, enquanto yaw depende de um mag ruidoso. **H4 quantifica.**

A "atitude" no teu log é, portanto, o **quatérnio do core primário do EKF3**
convertido para ângulos de Euler. **Não** é integração de giroscópio, **não** é
`atan2` do acelerômetro. É fusão sensorial.

---

## 1. O que é "atitude" dentro do log — e onde mais ela aparece

O `test_flight.m` só extrai `ATT` (Roll/Pitch/Yaw em graus). Mas o log bruto (`.bin`)
guarda **a atitude em várias formas**, que vamos usar para validação cruzada:

| Mensagem | Conteúdo | Papel |
|---|---|---|
| `ATT` | Roll, Pitch, Yaw (deg) + DesRoll/DesPitch/DesYaw | **Saída oficial** do AHRS (= EKF3 core primário) — o que o controlador usa |
| `XKQ` | Q1..Q4 | **Quatérnio** cru do EKF3 (por core) — a representação interna |
| `XKF1` | Roll, Pitch, Yaw, VN/VE/VD, PN/PE/PD, **GX/GY/GZ** (viés do giro), OH | Estado completo do EKF3, **por core** |
| `AHR2` | Roll, Pitch, Yaw | AHRS de **backup** (DCM — filtro complementar, não-EKF) |
| `DCM` | atitude + erros | Solução DCM pura (referência independente do EKF) |

Ter `XKQ` + `XKF1` + `AHR2` + `DCM` é ótimo: dá pra checar `ATT == euler(XKQ)`,
comparar core 0 vs core 1, e ver quanto o EKF diverge do DCM puro.

---

## 2. Por que não basta integrar o giroscópio (a motivação do Kalman)

Você já tem o `analyze_gyro_integration.m`, que mostra na prática: integrar `p,q,r`
puros **diverge** da atitude do EKF. Cada sensor isolado tem uma falha:

| Sensor | Dá o quê | Problema isolado |
|---|---|---|
| Giroscópio (`IMU.Gyr`) | velocidade angular `ω` | Integrar acumula **drift** = viés·t + random walk |
| Acelerômetro (`IMU.Acc`) | força específica `f = a − g` | Em repouso dá a vertical (gravidade) → roll/pitch; mas durante aceleração linear `a≠0` **mente** sobre a vertical; ruidoso/vibração |
| Magnetômetro (`MAG`) | campo magnético | Dá heading (yaw) absoluto; mas distorcido por correntes/motores → no teu drone o **compass 2 (externo) é o usado**, e é ruidoso |
| GPS (`GPS`) | posição/velocidade NED | Absoluto, sem drift; mas **lento** (~5–10 Hz), com atraso, e **não mede atitude diretamente** |
| Barômetro (`BARO`) | pressão → altitude | Boa altitude relativa; deriva lenta, ruído |

A ideia do EKF: **fundir** todos, dando a cada um um peso proporcional à sua
confiabilidade, e estimar de quebra os **viéses** dos sensores. O giro entra como
"motor" da predição (alta taxa, suave) e os demais "ancoram" essa predição (baixa
taxa, sem drift). É um filtro complementar generalizado e ótimo.

---

## 3. Arquitetura do EKF3

**EKF = Extended Kalman Filter.** "Extended" porque o modelo é **não-linear** (a
cinemática do quatérnio e as projeções de gravidade/campo não são lineares); o EKF
lineariza localmente (Jacobianos `F`, `H`) a cada passo.

Estrutura de duas etapas, repetidas em loop:

```
        ┌───────────────────────────────────────────────┐
  IMU → │  PREDIÇÃO (strapdown INS, ~alta taxa)          │
 (ω,a)  │   x⁻ = f(x, ω, a)     (propaga estado)         │
        │   P⁻ = F P Fᵀ + Q     (cresce a incerteza)     │
        └───────────────────────────────────────────────┘
                          │
        ┌─────────────────▼─────────────────────────────┐
GPS,    │  CORREÇÃO (a cada medida disponível)           │
Baro,   │   y = z − h(x⁻)            (inovação)          │
(Mag),  │   S = H P⁻ Hᵀ + R         (incerteza da inov.) │
GSF-yaw │   K = P⁻ Hᵀ S⁻¹           (ganho de Kalman)    │
        │   x⁺ = x⁻ + K y           (corrige estado)     │
        │   P⁺ = (I − K H) P⁻       (reduz a incerteza)  │
        └────────────────────────────────────────────────┘
```

Detalhes específicos do EKF3 do ArduPilot:

- **Fusion horizon / buffers de atraso.** Cada sensor tem latência diferente (GPS
  ~100–200 ms, baro/mag menos). O EKF3 mantém um buffer e funde cada medida no
  instante correto do passado, depois "re-avança" o estado. Por isso o `ATT` é
  causal e consistente apesar dos atrasos.
- **Predição roda a cada amostra de IMU** (após o filtro `INS_GYRO_FILTER`/
  `INS_ACCEL_FILTER` = 20 Hz no teu caso — então o que está no `IMU` do log **já
  vem filtrado a 20 Hz**, detalhe importante para a fase 2).
- **Correções rodam quando cada sensor entrega** (GPS lento, baro médio, etc.).

---

## 4. O vetor de estados (24 estados)

O EKF3 estima simultaneamente **24 estados**. Mapeando para o que está logado:

| # | Estado | Símbolo | Onde aparece no log |
|---|---|---|---|
| 1–4 | Quatérnio (NED→corpo) | `q0..q3` | `XKQ.Q1..Q4` → vira Roll/Pitch/Yaw |
| 5–7 | Velocidade NED | `vN,vE,vD` | `XKF1.VN/VE/VD` |
| 8–10 | Posição NED | `pN,pE,pD` | `XKF1.PN/PE/PD` |
| 11–13 | **Viés do giroscópio** | `b_g` | `XKF1.GX/GY/GZ` |
| 14–16 | **Viés do acelerômetro** | `b_a` | `XKF2.AX/AY/AZ` |
| 17–19 | Campo magnético terrestre NED | `m_E` | `XKF2.MN/ME/MD` |
| 20–22 | Campo magnético no corpo (hard-iron) | `m_B` | `XKF2.MX/MY/MZ` |
| 23–24 | Vento NE | `w_N,w_E` | `XKF2.VWN/VWE` |

**Insight direto para o teu projeto:** os estados 11–13 (`XKF1.GX/GY/GZ`) são
**o viés do giroscópio estimado online pelo próprio EKF**. É exatamente o que o teu
`estimate_bias.m` tenta calcular por janelas calmas. Na fase 2 vamos **comparar os
dois** (H2/H5) — o EKF te dá a resposta "oficial" de graça.

> Como o compass 2 está ativo, os estados 17–22 (campo magnético terrestre e do
> corpo) **são observados e atualizados** pelo mag — dá pra ver a evolução deles em
> `XKF2` (`MN/ME/MD`, `MX/MY/MZ`) e as inovações correspondentes em `XKF3`.

---

## 5. Passo de predição — giro e acelerômetro são ENTRADAS, não medidas

Ponto conceitual que confunde muita gente: no EKF de navegação, **o giroscópio e o
acelerômetro não são "medidas" no sentido do passo de correção** — eles são as
**entradas de controle `u`** que movem o modelo (mecanização strapdown INS):

**Atitude (quatérnio):**
```
Δθ = (ω_medido − b_g) · Δt          ← giro menos viés estimado
q⁻ = q ⊗ Δq(Δθ)                      ← propaga rotação (produto de quatérnios)
```

**Velocidade e posição:**
```
f_corpo = a_medido − b_a              ← força específica menos viés
v⁻ = v + ( R(q)·f_corpo + g_NED )·Δt  ← gira p/ NED, soma gravidade, integra
p⁻ = p + v·Δt
```

onde `R(q)` é a matriz de rotação corpo→NED e `g_NED = [0,0,+9.81]`.

**Covariância:**
```
P⁻ = F·P·Fᵀ + Q
```
`F = ∂f/∂x` é o Jacobiano (linearização). `Q` é o **ruído de processo** — quanto a
predição "perde confiança" por amostra. É aqui que entram teus parâmetros:

| Parâmetro (teu valor) | Estado afetado | Efeito |
|---|---|---|
| `EK3_GYRO_P_NSE = 0.015` rad/s | atitude | quanto o ruído do giro infla a incerteza de atitude |
| `EK3_ACC_P_NSE = 0.35` m/s² | velocidade | idem para o accel |
| `EK3_GBIAS_P_NSE = 0.001` rad/s | viés giro | random walk: quão rápido o viés do giro pode mudar |
| `EK3_ABIAS_P_NSE = 0.02` m/s² | viés accel | idem viés do accel |
| `EK3_WIND_P_NSE = 0.2` | vento | quão rápido o vento pode mudar |

Quanto **maior** o `P_NSE`, **menos** o EKF confia na predição → ele se apoia mais
nas medidas (mais "ágil", mais ruidoso). Quanto **menor**, mais suave e mais
dependente da qualidade do giro/accel.

---

## 6. Passo de correção — onde as medidas "puxam" o estado

Quando chega uma medida `z` (posição GPS, velocidade GPS, altitude baro, yaw do GSF),
o EKF calcula:

```
y = z − h(x⁻)              inovação (o quanto a medida discorda da predição)
S = H·P⁻·Hᵀ + R            incerteza da inovação  (H = ∂h/∂x)
K = P⁻·Hᵀ·S⁻¹              ganho de Kalman
x⁺ = x⁻ + K·y              corrige TODOS os estados correlacionados
P⁺ = (I − K·H)·P⁻          encolhe a incerteza
```

O `R` é o **ruído de medida** (os `*_M_NSE`) — teus valores reais:

| Parâmetro (teu valor) | Medida | Interpretação |
|---|---|---|
| `EK3_VELNE_M_NSE = 0.3` m/s | vel. horizontal GPS | desvio-padrão assumido |
| `EK3_VELD_M_NSE = 0.5` m/s | vel. vertical GPS | |
| `EK3_POSNE_M_NSE = 0.5` m | pos. horizontal GPS | |
| `EK3_ALT_M_NSE = 2.0` m | altitude baro | baro é "frouxo" (R alto) → corrige devagar |
| `EK3_MAG_M_NSE = 0.05` Gauss | magnetômetro | **irrelevante aqui** (compass off) |
| `EK3_YAW_M_NSE = 0.5` rad | yaw (GSF/externo) | ~28° de 1σ → yaw entra "fraco" |

**A intuição do ganho de Kalman** (caso escalar): `K = P/(P+R)`.
- `R` grande (medida ruim) → `K→0` → ignora a medida, confia no modelo.
- `R` pequeno (medida boa) → `K→1` → segue a medida.
- `P` grande (modelo incerto) → `K→1` → aceita a medida.

É literalmente um **filtro complementar com peso ótimo e adaptativo** (o peso varia
porque `P` evolui).

**Gating (rejeição de outliers):** antes de aceitar, o EKF testa a inovação
normalizada `yᵀ S⁻¹ y`. Se passar de um limiar, a medida é **rejeitada** (ex.: glitch
de GPS, `EK3_GLITCH_RAD = 25`). Esses "test ratios" estão logados em **`XKF4`**:
`SV` (velocidade), `SP` (posição), `SH` (altura), `SM` (mag). Valor `<1` = saudável,
`>1` = rejeitada. Ótimo para diagnosticar.

---

## 7. O ponto-chave: como cada ângulo fica observável

Esta é a parte mais importante para o teu sysid — e a mais contra-intuitiva.

### Roll e Pitch (tilt) — observáveis via **velocidade**, não "accel = gravidade"
A intuição comum ("o EKF usa o accel como referência de gravidade para roll/pitch")
está **certa só no repouso / alinhamento inicial**. Em voo:

- Um **erro de tilt** faz a projeção `R(q)·f_corpo` ficar errada → a gravidade não é
  cancelada corretamente → aparece uma **aceleração espúria** → integra em **erro de
  velocidade** → o **GPS** (vel. NE/D) vê esse erro e o corrige → e, pela
  **correlação cruzada na covariância `P`**, a correção "volta" e ajusta o tilt.

Ou seja: roll/pitch são corrigidos **indiretamente pela fusão de velocidade do GPS**,
não por uma medida direta de gravidade. O accel "puro" só domina o tilt quando a
aceleração linear é ~0 (repouso/hover steady) — e é por isso que durante manobras a
atitude do EKF **diverge** do `atan2` do acelerômetro. **Isso é a hipótese H3.**

### Yaw — vem do **magnetômetro (compass 2)**, com o GSF de backup
O yaw é o estado menos observável. No teu drone ele vem do **magnetômetro externo
(compass 2)**: o EKF projeta o campo magnético terrestre estimado no corpo e compara
com a leitura do mag (inovação `XKF3.IMX/IMY/IMZ`), corrigindo o heading. Em paralelo,
o EKF3 **sempre** roda o **EKF-GSF (Gaussian Sum Filter)** — um banco de mini-filtros
que deriva yaw da relação entre **velocidade do GPS** e **aceleração no corpo** durante
manobras (`XKY0/XKY1`); ele serve de *backup* e para re-alinhar o yaw se o mag ficar
inconsistente, mas não é a fonte primária aqui.

Mesmo com mag, o yaw tende a ser o canal mais fraco:
- roll/pitch são amarrados **com força** pela velocidade do GPS (1σ de 0.3 m/s);
- yaw depende de um mag **ruidoso** (inovações de 40–66 mGauss) com 1σ de medida de
  yaw de `EK3_YAW_M_NSE=0.5 rad` (~28°) → correção "mole";
- em hover puro o GSF quase não ajuda (sem velocidade) → sobra só o mag.

Isso casa com "yaw é o canal fraco" do teu modelo (`project_quad_sysid`). **H4** olha
`XKF3`/`XKY` e compara `ATT.Yaw` com o heading do mag (tilt-compensado) e com o GSF.

---

## 8. Resumo dos "pesos" — Q e R reais do teu drone

```
RUÍDO DE PROCESSO Q (confiança no modelo/predição — IMU):
   gyro      EK3_GYRO_P_NSE  = 0.015  rad/s
   accel     EK3_ACC_P_NSE   = 0.35   m/s²
   viés gyro EK3_GBIAS_P_NSE = 0.001  rad/s   (random walk)
   viés acc  EK3_ABIAS_P_NSE = 0.02   m/s²
   vento     EK3_WIND_P_NSE  = 0.2

RUÍDO DE MEDIDA R (confiança em cada sensor — maior = confia menos):
   GPS vel NE  EK3_VELNE_M_NSE = 0.3  m/s
   GPS vel D   EK3_VELD_M_NSE  = 0.5  m/s
   GPS pos NE  EK3_POSNE_M_NSE = 0.5  m
   baro alt    EK3_ALT_M_NSE   = 2.0  m
   mag         EK3_MAG_M_NSE   = 0.05 Gauss   (compass 2, usado p/ yaw)
   yaw         EK3_YAW_M_NSE   = 0.5  rad
```

A razão `Q/R` define o equilíbrio: o tilt segue o giro em alta frequência e é
"puxado" de volta pelo GPS em baixa frequência. Dá pra estimar uma **frequência de
crossover equivalente** de um filtro complementar a partir desses números (H5).

---

## 9. Estimação de viés online

O EKF não só usa o giro/accel — ele **aprende os viéses deles em tempo real**
(estados 11–16). Por isso a atitude do EKF não deriva como a integração pura:
o termo `(ω − b_g)` tem o viés continuamente removido.

- `XKF1.GX/GY/GZ` = viés do giro estimado (provável unidade: deg/s — **confirmar
  escala na fase 2**).
- `XKF2.AX/AY/AZ` = viés do accel estimado.

Comparar isso com o teu `estimate_bias.m` é um teste de consistência forte (H2).

---

## 10. Redundância: 2 cores + lane switching

`EK3_IMU_MASK=3` → o EKF3 roda **duas instâncias completas em paralelo**, uma com
IMU0 e outra com IMU1 (daí `XKF*_0` e `XKF*_1`). Cada core tem um "score" de saúde
(baseado em inovações/variâncias). O **core primário** (logado em `XKF4.PI`) é o que
alimenta o `ATT`. Se o primário degrada (ex.: clipping de IMU, vibração), o sistema
**troca de lane**. Vale checar se houve troca durante teus voos (afeta continuidade
do `ATT`).

---

## 11. Hipóteses para a Fase 2 (código) — mapeadas aos dados

Cada hipótese vira um script em `8_atitude_ekf_test/`. Tudo é verificável **com os
dados que você já tem** no `.bin` (precisamos re-extrair alguns campos além do que o
`test_flight.m` pega hoje — ver nota no fim).

| # | Hipótese | Como testar | Dados |
|---|---|---|---|
| **H1** | `ATT` == Euler do quatérnio do core primário | converter `XKQ` → euler e comparar com `ATT`; comparar core 0 vs core 1 | `ATT`, `XKQ`, `XKF1`, `XKF4.PI` |
| **H2** | Integrar `(giro − viés_EKF)` reduz drift vs giro puro | repetir `analyze_gyro_integration.m` subtraindo `XKF1.GX/GY/GZ`; comparar com `estimate_bias.m` | `IMU`, `XKF1`, `ATT` |
| **H3** | Roll/Pitch ≈ `atan2` do accel **só** em janelas quase-estáticas; diverge sob aceleração | calcular tilt do accel e comparar com `ATT`, segmentando por `‖a‖≈g` | `IMU`, `ATT` |
| **H4** | Yaw vem do **mag (compass 2)**; GSF é backup | comparar `ATT.Yaw` com o heading do mag (tilt-compensado) e com `XKY` (GSF); ver inovações `XKF3.IMX/IMY/IMZ` e ratio `SM` | `MAG`, `XKF3`, `XKY0/1`, `ATT` |
| **H5** | Reconstruir Q/R → ganho de Kalman / freq. de crossover do complementar equivalente | montar filtro complementar 1ª/2ª ordem com os `*_NSE` e comparar com `ATT` | `PARM`, `IMU`, `ATT` |
| **H6** | Inovações e test ratios contam quando cada sensor corrige/rejeita | plotar `XKF3` (inovações) e `XKF4` (`SV/SP/SH/SM`) ao longo do voo | `XKF3`, `XKF4` |
| **H7** (esticado) | Um mini-EKF/complementar em MATLAB alimentado com IMU+GPS reproduz o `ATT` | implementar versão reduzida e comparar | `IMU`, `GPS`, `BARO`, `ATT` |

### Nota sobre extração de dados
O `load_log_data.m`/`test_flight.m` hoje só pega `IMU/ATT/RCOU/RCIN/MAG/GPS(tempo)`.
Para a Fase 2 vamos precisar **estender a extração** para incluir `XKQ`, `XKF1`,
`XKF2`, `XKF3`, `XKF4`, `XKY0/1`, `GPS_0` (com vel/curso) e `BARO_0`. Faremos um
loader dedicado em `8_atitude_ekf_test/` para não mexer no pipeline de sysid.

---

## 12. Conclusão da fase teórica

1. **Sim** — o `Roll/Pitch/Yaw` do teu `ATT` é a saída do **EKF3** (quatérnio do core
   primário), não integração de giro nem `atan2` de accel.
2. O EKF3 é um filtro de **24 estados**, **quaterniônico**, com **predição strapdown
   INS** (giro+accel como entradas) e **correção multi-sensor** (GPS + baro + yaw).
3. Os "pesos" são `Q` (`*_P_NSE`) e `R` (`*_M_NSE`); o ganho de Kalman os equilibra
   adaptativamente via a covariância `P`.
4. **Tilt** (roll/pitch) é observável via **GPS-velocidade**; **yaw** vem do
   **magnetômetro (compass 2, externo)** — com o **GSF** rodando de backup. O mag é
   ruidoso, o que ajuda a explicar yaw ser o canal mais fraco.
5. Temos no log **tudo** para validar: o estado interno (`XKF1/2`), o quatérnio
   (`XKQ`), as inovações (`XKF3`), os test ratios (`XKF4`) e o estimador de yaw
   (`XKY`).

Próximo passo: dizer "bora codar" e eu começo pela extensão do loader + **H1** (a
sanity check `ATT == euler(XKQ)`).
