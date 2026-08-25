# 8_atitude_ekf_test/

Estudo: **de onde vem a atitude (`ATT`) do log do Pixhawk** e como o EKF3 do
ArduPilot a calcula a partir das medidas brutas (giro, accel, GPS, baro, mag).

## Fase 1 — Teoria ✅
- [`teoria_ekf3_atitude.md`](teoria_ekf3_atitude.md) — análise completa: arquitetura
  do EKF3, os 24 estados, predição/correção, os pesos (Q/R) **reais do drone**, e por
  que tilt vem do GPS e yaw vem do GSF (compass desligado).

## Scripts prontos (rodáveis agora)

Trabalham direto no `../1_data/logs_concat.mat`, janela **1–110 s**. Standalone.

| Script | O que faz | Resultado |
|---|---|---|
| [`ekf_atitude.m`](ekf_atitude.m) | **Recria o EKF** (quatérnio + viés de giro; predição pelo giro + correção de tilt pelo accel + correção de yaw pelo mag, com gating). Compara com `ATT`. | roll/pitch RMS ≈ **1.8°**; yaw RMS ≈ **12°** (limitado pelo mag, sem deriva). `USE_MAG=false` → yaw deriva ~24° |
| [`integra_giro_euler.m`](integra_giro_euler.m) | **Integração pura** da cinemática de Euler com p,q,r, realimentando a própria atitude (sem correção). Compara com `ATT`. | deriva crescente (roll ~ -74° = viés do giro × t) |

> ℹ️ `logs_concat.mat` agora inclui **MAG (compass 2) + velocidade do GPS** (via
> `test_flight.m`). O EKF usa giro+accel+mag → reproduz os **3 eixos**. O yaw fica
> *limitado* (~12°), não perfeito, porque o mag deste drone é ruidoso (~20° de
> espalhamento de heading). A velocidade do GPS ainda **não** é fundida — exigiria um
> INS-EKF com estados de velocidade (próximo nível); em hover ela ajudaria pouco.

## Fase 2 — Código (resto, a fazer)
Cada hipótese da teoria vira um script. Ordem sugerida:

| Script (planejado) | Hipótese | O que valida |
|---|---|---|
| `load_ekf_data.m` | — | loader que extrai `XKQ/XKF1..4/XKY/GPS/BARO` (além do que o `load_log_data.m` pega) |
| `h1_att_vs_quaternion.m` | H1 | `ATT` == euler(`XKQ`)? core 0 vs core 1? |
| `h2_gyro_bias.m` | H2 | integrar `(giro − viés_EKF)` vs giro puro; vs `estimate_bias.m` |
| `h3_accel_tilt.m` | H3 | roll/pitch do accel (`atan2`) vs `ATT`, segmentado por `‖a‖≈g` |
| `h4_yaw_source.m` | H4 | yaw vem do mag (compass 2) ou do GSF? quão ruidoso? |
| `h5_complementary.m` | H5 | filtro complementar com os pesos Q/R reais vs `ATT` |
| `h6_innovations.m` | H6 | inovações (`XKF3`) e test ratios (`XKF4`) ao longo do voo |
| `h7_mini_ekf.m` | H7 | mini-EKF/complementar reproduz `ATT`? (esticado) |

## Configuração do drone (extraída do log)
- ArduPilot **4.6.3** QuadPlane, `AHRS_EKF_TYPE=3` (EKF3), 2 cores (IMU0+IMU1).
- Fontes: GPS (pos/vel horiz + vel vert), Baro (alt), **mag compass 2 (externo)** → yaw; GSF de backup.
- Pesos `Q`/`R`: ver tabela na seção 8 da teoria.

## Dados
Os campos do EKF interno (`XKF*`, `XKQ`, `XKY*`) estão nos `.bin` brutos em
`../1_data/` (ex.: `5 ...bin-155264.mat`), **não** no `logs_concat.mat` atual.
A Fase 2 começa estendendo a extração.
