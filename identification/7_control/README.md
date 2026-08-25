# 7_control — Projeto de controle (cascata PID)

Pega o modelo linear (`outputs/linear_model.mat`, entrada `[T, Mx, My, Mz]`,
9 estados) e aplica um **controle em cascata** com setpoints de **altitude**,
**velocidade forward** e **heading**. Primeiro controle: foco em **validar a
estrutura**, não em atender requisitos de desempenho.

## Estrutura de controle (forward via pitch; roll regulado a zero)

```
altitude h     → T        (PD + I + feedforward m·g)
vel. forward u → θc → My  (PI externo de velocidade + PD interno de pitch)
heading ψ      → Mz       (PD com wrap ±180°)
roll φ = 0     → Mx       (PD: mantém nivelado)
```

> Obs.: como a velocidade forward é gerada inclinando em **pitch**, o **roll**
> é que fica regulado em zero (não o pitch). A velocidade lateral fica ~0.

## Arquivos

| arquivo | papel |
|---|---|
| `design_control.m` | projeta os ganhos por alocação de polos → `control_gains.mat` |
| `control_law.m` | a lei de controle em cascata (usada nas simulações `.m`) |
| `forces_to_pwm.m` | **a ponte** (control allocation) `[T,M] → PWM`: modo `linear` (M_mix) ou `nonlinear` (alocação em empuxo + fT⁻¹) |
| `sim_control_linear.m` | malha fechada no **modelo linear** (Diagrama A) → figura |
| `sim_control_nonlinear.m` | mesmo controlador na **planta não-linear** via ponte (Diagrama B) vs linear → figura |
| `quad_control_loops.slx` | **modelo Simulink**: UM BLOCO POR CONTROLADOR (5 subsistemas) + planta State-Space linear |
| `run_control_demo.m` | roda tudo em sequência |

## Modelo Simulink (`quad_control_loops.slx`)

Cinco subsistemas de controlador **separados** no nível de topo, cada um com
blocos visíveis dentro (Soma / Ganho / Integrador):

| bloco | entradas | saída |
|---|---|---|
| `Ctrl_Altitude` | h_sp, h, w | δT |
| `Ctrl_VelFwd`   | u_sp, u | θc |
| `Ctrl_Pitch`    | θc, θ, q | My |
| `Ctrl_Roll`     | φ, p | Mx |
| `Ctrl_Heading`  | ψ_sp, ψ, r | Mz |

Saídas → `Mux [δT, Mx, My, Mz]` → `Planta_Linear` (bloco State-Space, modelo
linear aumentado com altitude h: `ḣ = -w`, 10 estados) → realimentação.

## Como rodar

```matlab
>> run_control_demo        % projeto → sim linear → sim NL → Simulink
```

ou passo a passo:

```matlab
>> design_control                 % ganhos
>> sim_control_linear             % valida no modelo linear (.m)
>> sim_control_nonlinear          % valida na planta não-linear (.m) vs linear
>> open_system('quad_control_loops'); sim('quad_control_loops')   % Simulink
```

## Cenário de teste (degraus)

| t (s) | setpoint |
|---|---|
| 1  | altitude → 2 m |
| 8  | vel. forward → 3 m/s |
| 16 | heading → 30° |

## Resultado (consistência)

`.m linear` ≡ `.m não-linear` ≡ `Simulink`: altitude 2.00 m, vel 3.00 m/s,
heading 30.0°, erro de regime ~0. O `.m` não-linear mostra acoplamentos físicos
(perda de altitude ao inclinar — efeito cos θ) que a malha de altitude corrige.
