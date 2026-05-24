# identification/

Identificação de parâmetros do modelo da planta VTOL híbrido (DH2) a partir do log de voo experimental.

## Quickstart

```matlab
% No MATLAB:
cd /Users/graest/ita-master/artigo/artigo-dh-quad/identification
setup_paths();        % bootstrap: adiciona todas as subpastas ao path

identify_plant        % roda identificação completa (gera outputs/P_identified.mat)
compare_models('all') % roda 3 comparações (gera figuras em outputs/images/)
linearize             % linearização (gera outputs/linear_model.mat)
```

## Estrutura (por fase do workflow)

```
identification/
├── README.md                  ← este arquivo
├── setup_paths.m              ← resolve paths e adiciona ao MATLAB path
│
├── 1_data/                    ← Insumos experimentais
│   └── log_data.mat           (log de voo do drone)
│
├── 2_model/                   ← Equações da planta
│   ├── vtol_dynamics.m        (modelo ODE — 3/9/17 estados)
│   ├── motor_models.m         (spline PWM→T, PWM→Q)
│   └── parameters.m           (constantes, P0, bench data centralizados)
│
├── 3_identification/          ← Otimização de parâmetros
│   ├── identify_plant.m       (script principal — EEM + OEM multi-segmento)
│   └── eem_cost_function.m    (custo Fase A)
│
├── 4_simulink/                ← Modelo Simulink e setup
│   ├── quad_model_v4.slx      (modelo Simulink atual)
│   ├── setup_quad_v4.m        (prepara workspace pra simular)
│   └── P_J_to_simulink.m      (converte P_J (24) → P_estimated (23))
│
├── 5_validation/              ← Comparações vs medido / vs cenários / vs linear
│   └── compare_models.m       (3 modos: log, scenarios, linear)
│
├── 6_linear/                  ← Linearização
│   ├── linearize.m            (Jacobianos numéricos)
│   └── vtol_dynamics_linearized.m
│
├── outputs/                   ← Tudo gerado por scripts
│   ├── P_identified.mat       (parâmetros identificados — output de identify_plant)
│   ├── linear_model.mat       (output de linearize)
│   └── images/                (figuras de validação)
│
├── docs/                      ← Documentação técnica
│   ├── new_identification.md
│   └── vtol_dynamics.md
│
├── reference/                 ← Material externo (não modificar)
│   └── mirko/                 (modelo de comparação do Mirko)
│
└── legacy/                    ← Arquivado, sem uso atual
    ├── create_thrust_model.m  (substituído por motor_models.m)
    ├── create_torque_model.m  (substituído por motor_models.m)
    ├── simulate.m             (substituído por uso direto vtol_dynamics + ode45)
    ├── motor_mapping_check.m  (diagnóstico antigo)
    └── ... (identification.m, main.m, model.m, etc — versões iniciais)
```

## Workflow típico

```
        log_data.mat (1_data/)
                │
                ▼
      identify_plant.m ◄────── parameters.m, motor_models.m, vtol_dynamics.m
        (3_identification/)
                │
                ▼
      P_identified.mat (outputs/)
                │
                ├──► setup_quad_v4 → simular quad_model_v4.slx
                │
                ├──► compare_models('log'/'scenarios'/'linear')
                │
                └──► linearize → linear_model.mat (outputs/)
```

## Dependências entre pastas

| De | Chama / Carrega | Em |
|---|---|---|
| `3_identification/identify_plant.m` | `log_data.mat`, motor_models, vtol_dynamics, parameters, eem_cost_function | 1_data, 2_model, 3_identification |
| `4_simulink/setup_quad_v4.m` | `log_data.mat`, `P_identified.mat`, P_J_to_simulink, parameters | 1_data, outputs, 4_simulink, 2_model |
| `5_validation/compare_models.m` | quad_model_v4.slx, vtol_dynamics, motor_models, P_J_to_simulink, linear_model.mat | 4_simulink, 2_model, outputs |
| `6_linear/linearize.m` | vtol_dynamics, motor_models | 2_model |

Todas resolvidas via `setup_paths()` (sem hardcoded paths).

## Parâmetros do drone (fonte oficial: docs/arquitetura_pa.pdf, tabela 1.0)

| Parâmetro | Valor | Status no código |
|---|---|---|
| Massa total | 2.20 kg | ✅ atualizado em parameters.m (era 1.6011) |
| Ix (roll) | 0.14410 kg·m² | ⚠️ código tem 0.063 — task #73 pendente |
| Iy (pitch) | 0.11550 kg·m² | ⚠️ código tem 0.250 |
| Iz (yaw) | 0.25716 kg·m² | ⚠️ código tem 0.116 |
| Asa S | 0.27 m² | ⚪ documentado, ainda não usado |
| Asa b | 1.20 m | ⚪ documentado |
| Asa c | 0.226 m | ⚪ documentado |
| Airfoil | USA-35B | ⚪ documentado |
| Z_CG | -0.05 m | ⚠️ task #74 pendente (efeito de pêndulo) |

## Tarefas pendentes (ver task list)

- #65: Adicionar M_y_wing (momento de pitching da asa em hover)
- #66: Re-identificar com priors corretos (k_T fixo, regularização)
- #67: Estender vtol_dynamics com aerodinâmica de asa fixa
- #68: Linearizar em múltiplos pontos de operação (transição)
- #69-70: Controlador + simulação de takeoff/transição
- #73: Reconciliar valores de inércias (slide XFLR5 vs código)
- #74: Adicionar Z_CG ao modelo
