# 9_FEM_identification — Identificação por Filter-Error Method (FEM)

Identificação do quad **do zero**, baseada em:
- **Jategaonkar (2015)**, *Flight Vehicle System Identification: A Time-Domain Methodology* — Cap. 5 (FEM) e Cap. 9 (Unstable / closed-loop).
- **Tischler & Remple (2006)**, *Aircraft and Rotorcraft System Identification* — coerência (Eq. 7.54), domínio da frequência, malha fechada.

Mantém os **dados e janelas atuais**; é uma pasta separada (não toca em `3_identification`).

---

## Estrutura Quad-M (Jategaonkar, Fig. 1.5)

### Decisão-raiz: open-loop plant ID com entrada exógena (Cap. 9, §9.2 / Fig. 9.1)
Identificamos a **planta H** (`u → y`): `u` = PWM dos motores (RCOU) → forças; `y` = saídas da IMU.
Não é o modelo "lumped" `d_p → y`. O **RCIN (doublet do piloto)** é a **entrada exógena** (independente da saída) que valida a condição de malha fechada do Cap. 9.

### M1 — Maneuver
- Reusa as janelas atuais (ver abaixo). Doublets de roll/pitch já presentes (RCIN), `Δt≈0.45 s` ↔ `f_C≈0.67 Hz` (Eq. 2.10/2.12 do Jategaonkar; bate com a PSD).
- Yaw é fraco (excitação estreita) — fica documentado; voo dedicado é fase futura.

### M2 — Measurements
- IMU (gyro/accel), ATT (EKF), RCOU (`u`), RCIN (`d_p`), MAG, GPS.
- Amostragem 10 Hz (Nyquist ok p/ ≤2 Hz; log rápido fica p/ voo futuro).
- **Data Compatibility Check (Cap. 10)** como pré-passo (escala/bias/atraso de sensor) — a confirmar.

### M3 — Methods: FILTER-ERROR METHOD (Cap. 5 + §9.8)
Estimador de estado preditor+corretor com ganho de filtro K (Eqs. 9.33–9.36):
```
ỹ(tk)   = g[x̃(tk), u(tk), b]                 (obs. predita)
x̂(tk)   = x̃(tk) + K[z(tk) − ỹ(tk)]           (correção pela medida)
x̃(tk+1) = x̂(tk) + ∫ f[x,u,b] dt              (predição/propagação)
```
- O feedback `K[z−ỹ]` **estabiliza a integração** → roda em malha fechada/instável sem divergir (§9.8, p.350). K **não** se relaciona ao feedback do controlador real.
- Estima **ruído de processo (rajada/turbulência) + ruído de medida + bias** simultaneamente → o resíduo colorido de ~0.65 Hz pode ser tratado como **process noise**, reduzindo viés de parâmetro.
- Ressalva do livro (§9.8): *"good response match tends to mask the modeling discrepancies"* → validação rigorosa é essencial.

### M4 — Models (Cap. 3)
- **Estado x = [p, q, r, u, v, w]** (6). Atitude = **entrada medida** (projeta gravidade), NÃO estado/saída.
- **Observação y = [p, q, r, aₓ, a_y, a_z]** (6): taxas (estados) + acelerações (saída algébrica do `accelerometer_model`).
- Entrada: PWM (→ T/Mx/My/Mz via fT/fQ) + atitude medida.
- Parâmetros: reusa `parameters.m` (inércia CAD, k_T bancada) + bias `b_x,b_y` (§3.5) + covariâncias de ruído (FEM).
- Reusa `vtol_dynamics` e `accelerometer_model`.

### Validation (Cap. 11)
- Janela held-out `[605,625]`.
- **Validar em: p, q, r, aₓ, a_y, a_z.** Atitude NÃO é validada (é entrada medida — ver nota abaixo).
- Métricas: TIC/Theil, coerência (Tischler Eq. 7.54 γ²≥0.6), brancura (ACF), CRB.

### Por que NÃO validar atitude
No FEM/OEM `y` é o que se mede e quer casar. Atitude:
1. vem do **EKF** (não sensor bruto) e tem **yaw fracamente observável**;
2. integrada em malha aberta **deriva**;
3. usada como **entrada** (gravidade) elimina o drift e isola a física.
→ Observação = IMU `[p,q,r,aₓ,a_y,a_z]`, atitude = entrada exógena medida.

---

## Janelas (reaproveitadas do `3_identification/identify_plant.m`)
- **Treino (5):** `[4,24] [25,41] [42,62] [63,99] [100,125]` s  (~120 s)
- **Validação:** `[605,625]` s

---

## Roadmap de scripts (a construir nesta pasta)
| fase | script | papel |
|---|---|---|
| 0 | `show_windows.m` | imprime + plota as janelas de treino/validação (sanity) |
| 1 | `fem_config.m` ou header | janelas, dt, SG, ruídos, bounds, fonte de P |
| 2 | `identify_fem.m` | FEM: predict/correct (K), estima θ + bias + ruído de processo |
| 3 | `validate_fem.m` | validação em p,q,r,acc (TIC, coerência, brancura, CRB) |
| 4 | (futuro) voo dedicado de sysid (yaw) | M1 fase 2 |

Tudo lê dados via `setup_paths()` (logs em `../1_data`), sem duplicar arquivos.
