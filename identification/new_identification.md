# new_identification.m

Script principal de identificação de parâmetros do VTOL (quadcóptero H-frame, modo quadricóptero). Utiliza dados de voo (log ArduPilot) para identificar 24 parâmetros do modelo dinâmico.

## Pipeline de Execução

```
┌─────────────────────────────────────────────────┐
│  1. Configuração de janelas (treino + validação) │
│  2. Carregamento e interpolação dos dados        │
│  3. Modelos de referência dos motores            │
│  4. Extração dos segmentos de treino/validação   │
│  5. Chute inicial P0 e bounds                    │
│  6. Validação com P0 (referência)                │
│  7. FASE A: EEM rotacional → P_eem              │
│  8. FASE B: OEM progressivo → P_final           │
│  9. Validação final + diagnósticos               │
│ 10. Resumo comparativo P0 vs P_final             │
└─────────────────────────────────────────────────┘
```

## Dados de Entrada

- **`log_data.mat`**: Log de voo ArduPilot contendo:
  - `IMU`: giroscópio (GyrX/Y/Z) e acelerômetro (AccX/Y/Z) a ~100 Hz
  - `ATT`: atitude estimada pelo EKF (Roll, Pitch, Yaw) em graus
  - `RCOU`: sinais PWM dos 4 motores (C1..C4) em µs
  - `GPS`: dados GPS (não usado diretamente na identificação)

- **Interpolação**: Todos os sinais são interpolados para uma base de tempo comum com `dt = 0.1s` (10 Hz).

## Vetor de Parâmetros (24 elementos)

```
P(1:4)   = [Jx, Jy, Jz, Jxz]        Inércias
P(5:8)   = [k_T1, k_T2, k_T3, k_T4]  Fatores de escala empuxo
P(9:12)  = [k_Q1, k_Q2, k_Q3, k_Q4]  Fatores de escala torque reativo
P(13:15) = [Dp, Dq, Dr]               Amortecimento aerodinâmico
P(16:18) = [Bp, Bq, Br]               Biases (torques residuais)
P(19:20) = [dx_cg, dy_cg]             Offset do CG vs CAD
P(21:24) = [Xu_m, Yv_m, Zw_m, Bz]    Parâmetros translacionais
```

## Fase A: EEM (Equation Error Method)

- **Otimiza**: P(1:20) — todos os parâmetros rotacionais + CG offsets
- **Método**: Algébrico (sem integração ODE). Compara derivadas numéricas dos dados medidos com derivadas previstas pelo modelo.
- **Dados**: Concatenação de todos os segmentos de treino.
- **Pesos**: Inverso da variância de cada derivada (`p_dot`, `q_dot`, `r_dot`).
- **Solver**: `lsqnonlin` com `trust-region-reflective`, 1500 iterações máximas.
- **Arquivo**: `eem_cost_function.m`

**Vantagens do EEM**: Rápido, convexo, bom ponto de partida. **Desvantagens**: Sensível a ruído nas derivadas numéricas, não minimiza o erro de predição.

## Fase B: OEM (Output Error Method) Progressivo

- **Otimiza**: Todos os 24 parâmetros simultaneamente.
- **Método**: Integração temporal (compara saída simulada vs medida).
- **Janelas progressivas**: `[1s, 2s, 3s]` — começa com janelas curtas e aumenta progressivamente.
- **Seleção**: O estágio com melhor R² médio rotacional na validação é mantido como `P_final`.
- **Pesos**: Inverso da variância de cada canal (p, q, r, AccX, AccY, AccZ).
- **Solver**: `lsqnonlin` com `trust-region-reflective`, 500 iterações por estágio.

### Integração no Custo OEM

| Componente | Método | Sub-steps | Precisão |
|------------|--------|-----------|----------|
| **Rotacional** (p,q,r) | RK4 sub-stepping | n_sub=5, dt_sub=0.02s | O(dt⁴) ≈ 10⁻⁷ |
| **Translacional** (u,v,w) | Euler forward | dt=0.1s, sem sub-step | O(dt) ≈ 10⁻¹ |

### Custo Combinado

```
e = [e_rot; e_acc]

e_rot = [w_p*(p_med - p_sim); w_q*(q_med - q_sim); w_r*(r_med - r_sim)]
e_acc = [w_ax*(accX_med - accX_modelo); w_ay*(...); w_az*(...)]
```

Onde `accX_modelo` é a **força específica** (sem gravidade):
```
accX_m = r*v - q*w + Xu*u          (não inclui gx)
accY_m = p*w - r*u + Yv*v
accZ_m = q*u - p*v - T/m + Zw*w + Bz
```

## Validação (simulate_full)

Função interna que simula o modelo completo para treino e validação:

1. **Rotacional**: `ode45` com 3 estados `[p, q, r]`, usando `vtol_dynamics.m`.
2. **Translacional**: `ode45` com 3 estados `[u, v, w]`, usando `vtol_dynamics.m` no modo `translational`.
   - Utiliza **p, q, r simulados** (do modelo rotacional) como entrada — modelo acoplado.
   - Atitude (phi, theta, psi) vem dos dados **medidos** (EKF do ArduPilot).
3. **Força específica**: Calculada subtraindo a gravidade da aceleração inercial:
   ```
   accX_s = u_dot - gx    onde gx = -g*sin(theta)
   accY_s = v_dot - gy    onde gy =  g*cos(theta)*sin(phi)
   accZ_s = w_dot - gz    onde gz =  g*cos(theta)*cos(phi)
   ```

## Funções Internas (Local Functions)

| Função | Descrição |
|--------|-----------|
| `oem_multi_seg_cost` | Wrapper: acumula resíduos de múltiplos segmentos |
| `oem_ms_cost_func` | Core: calcula resíduos rotacionais (RK4) + translacionais (Euler) |
| `simulate_full` | Simulação completa com ode45 para validação |
| `print_R2` | Imprime R² de treino e validação |
| `plot_all_results` | Gera gráficos p/q/r e AccX/Y/Z (treino + validação) |
| `plot_torques` | Diagnóstico: momentos Mx, My, Mz vs entradas PWM |
| `plot_forces` | Diagnóstico: forças translacionais + AccX/Y/Z simulado vs medido |

## Diagnósticos Gerados

### Gráficos (salvos em `images/`)
- `pqr_val_*.png` — Velocidades angulares: medido vs simulado
- `acc_val_*.png` — Acelerações (força específica): medido vs simulado
- `torques_validacao.png` — Torques Mx, My, Mz e entradas por motor
- `forcas_validacao.png` — Forças translacionais com modelo simulado
- `pwm_analise.png` — Análise de PWMs, empuxo por motor, e T_ref cru

### Análise de Assimetria PWM
Calcula coeficientes de variação e assimetrias (roll, pitch, yaw) para distinguir entre:
- **Assimetria constante** → diferença de CG ou desempenho entre motores
- **Assimetria variável** → perturbação externa (vento)

## Métricas

- **R²** (coeficiente de determinação): `1 - SS_res / SS_tot`
  - R² = 1.0: ajuste perfeito
  - R² = 0.0: modelo equivalente à média
  - R² < 0: modelo pior que a média

## Resultados Típicos (Estado Atual)

| Canal | R² Validação |
|-------|-------------|
| p (roll rate) | 0.325 |
| q (pitch rate) | 0.413 |
| r (yaw rate) | 0.706 |
| AccX | 0.173 |
| AccY | -0.057 |
| AccZ | -0.429 |

## Dependências

### Arquivos Ativos (usados pelo pipeline)

| Arquivo | Descrição |
|---------|-----------|
| `new_identification.m` | Script principal |
| `vtol_dynamics.m` | Função ODE da dinâmica do VTOL |
| `eem_cost_function.m` | Custo da Fase A (EEM algébrico) |
| `create_thrust_model.m` | Cria modelo polinomial PWM → empuxo (N) |
| `create_torque_model.m` | Cria modelo polinomial PWM → torque (N·m) |
| `log_data.mat` | Dados de voo (ArduPilot log) |

### Arquivos Legados (NÃO usados pelo pipeline atual)

| Arquivo | Descrição |
|---------|-----------|
| `identification.m` | Versão anterior do script de identificação |
| `main.m` | Script antigo (Live Script convertido) |
| `model.m` | Script antigo de simulação de atitude |
| `cost_function_vtol.m` | Função de custo antiga (single-shooting) |
| `oem_cost_function.m` | OEM cost antigo (multiple-shooting separado) |
| `thrustFromPWM.m` | Modelo de empuxo hard-coded (substituído por `create_thrust_model`) |
| `torqueFromPWM.m` | Modelo de torque hard-coded (substituído por `create_torque_model`) |
| `motor_mapping_check.m` | Script auxiliar para verificar mapeamento de motores |
| `quad_model_v3.slxc` | Cache de modelo Simulink (não usado) |
