# vtol_dynamics.m

Função ODE que implementa a dinâmica do VTOL (quadcóptero H-frame) para uso com integradores como `ode45`.

## Assinatura

```matlab
dydt = vtol_dynamics(t, y, P, pwm_time, pwm_signals, func_T_ref, func_Q_ref, constants)
```

## Modos de Operação

A função opera em três modos, selecionados pelo tamanho do vetor de estados `y` e pelo campo `constants.mode`:

| Modo | Estados | Vetor `y` | Uso |
|------|---------|-----------|-----|
| **Rotacional** | 3 | `[p; q; r]` | Identificação rotacional (EEM/OEM) e validação |
| **9 estados** | 9 | `[p; q; r; phi; theta; psi; u; v; w]` | Simulação completa (não usado atualmente) |
| **Translacional** | 3 | `[u; v; w]` | Validação translacional com `p,q,r` e atitude de fonte externa |

O modo translacional é ativado quando `constants.mode = 'translational'` e requer `constants.measured_data`.

## Vetor de Parâmetros P (24 elementos)

| Índice | Parâmetro | Descrição | Unidade |
|--------|-----------|-----------|---------|
| P(1) | Jx | Momento de inércia em roll | kg·m² |
| P(2) | Jy | Momento de inércia em pitch | kg·m² |
| P(3) | Jz | Momento de inércia em yaw | kg·m² |
| P(4) | Jxz | Produto de inércia cruzado | kg·m² |
| P(5:8) | k_T1..k_T4 | Fatores de escala de empuxo por motor | adim. |
| P(9:12) | k_Q1..k_Q4 | Fatores de escala de torque reativo por motor | adim. |
| P(13) | Dp | Amortecimento aerodinâmico em roll | 1/s |
| P(14) | Dq | Amortecimento aerodinâmico em pitch | 1/s |
| P(15) | Dr | Amortecimento aerodinâmico em yaw | 1/s |
| P(16) | Bp | Bias em roll (torque residual) | rad/s² |
| P(17) | Bq | Bias em pitch (torque residual) | rad/s² |
| P(18) | Br | Bias em yaw (torque residual) | rad/s² |
| P(19) | dx_cg | Offset CG longitudinal (frente +) | m |
| P(20) | dy_cg | Offset CG lateral (direita +) | m |
| P(21) | Xu_m | Derivada de arrasto em u (força/massa) | 1/s |
| P(22) | Yv_m | Derivada de arrasto em v (força/massa) | 1/s |
| P(23) | Zw_m | Derivada de arrasto em w (força/massa) | 1/s |
| P(24) | Bz | Bias de força em z (por massa) | m/s² |

## Geometria do H-Frame

```
         Frente
    M3(FL,CCW)  M1(FR,CW)
         \      /
          \    /
           ----
          /    \
         /      \
    M2(RL,CW)  M4(RR,CCW)
         Traseira
```

**Braços nominais (do CAD):**
- Roll (lateral): Ly = 0.232 m (simétrico)
- Pitch frente: Lx_f = 0.311185 m
- Pitch traseira: Lx_r = 0.342865 m (assimétrico)

**Braços efetivos com offset do CG:**
```matlab
Lx_r = 0.232 - dy_cg     % direita (motores 1,4)
Lx_l = 0.232 + dy_cg     % esquerda (motores 2,3)
Ly_f = 0.311185 - dx_cg   % frente (motores 1,3)
Ly_r = 0.342865 + dx_cg   % traseira (motores 2,4)
```

## Equações Rotacionais

Constantes G derivadas das inércias (corpo rígido com produto cruzado Jxz):

```
gamma0 = Jx*Jz - Jxz²
G1 = Jxz*(Jx - Jy + Jz) / gamma0
G2 = (Jz*(Jz - Jy) + Jxz²) / gamma0
G3 = Jz / gamma0
G4 = Jxz / gamma0
G5 = (Jz - Jx) / Jy
G6 = Jxz / Jy
G7 = (Jx*(Jx - Jy) + Jxz²) / gamma0
G8 = Jx / gamma0
```

**Momentos:**
```
Mx = -(Lx_r*T1 - Lx_l*T2 - Lx_l*T3 + Lx_r*T4)   [roll]
My =  Ly_f*T1 - Ly_r*T2 + Ly_f*T3 - Ly_r*T4      [pitch]
Mz =  Q1 + Q2 - Q3 - Q4                            [yaw: CW+ CCW-]
```

Onde `Ti = k_Ti * T_ref(PWMi)` e `Qi = k_Qi * Q_ref(PWMi)`.

**Acelerações angulares:**
```
p_dot = G1*p*q - G2*q*r + G3*Mx + G4*Mz - Dp*p + Bp
q_dot = G5*p*r - G6*(p² - r²) + (1/Jy)*My - Dq*q + Bq
r_dot = G7*p*q - G1*q*r + G4*Mx + G8*Mz - Dr*r + Br
```

## Equações Translacionais

Usa a matriz de rotação completa R_nb (NED → Body) para projetar a gravidade:

```
G_body = R_nb' * [0; 0; m*g]
Fx = G_body(1)                    % gravidade em x
Fy = G_body(2)                    % gravidade em y
Fz = -T_total + G_body(3)        % empuxo (eixo -z body) + gravidade em z

u_dot = r*v - q*w + Fx/m + Xu*u
v_dot = p*w - r*u + Fy/m + Yv*v
w_dot = q*u - p*v + Fz/m + Zw*w + Bz
```

## Modo Translacional (constants.mode = 'translational')

Integra apenas `[u; v; w]`, interpolando `p, q, r, phi, theta, psi` de dados medidos (ou simulados) fornecidos em `constants.measured_data`:

```matlab
constants.mode = 'translational';
constants.measured_data.time  = time_vector;
constants.measured_data.p     = p_data;     % pode ser medido ou simulado
constants.measured_data.q     = q_data;
constants.measured_data.r     = r_data;
constants.measured_data.phi   = phi_data;   % em radianos
constants.measured_data.theta = theta_data;
constants.measured_data.psi   = psi_data;
```

## Sub-funções Internas

| Função | Descrição |
|--------|-----------|
| `translational_mode()` | Implementa o modo translacional de 3 estados |
| `euler_kinematics()` | Cinemática de Euler: `[p,q,r]` → `[phi_dot, theta_dot, psi_dot]` |
| `translational_eqs()` | Equações dinâmicas translacionais (gravidade via R_nb completa) |

## Entradas de Motor

As funções `func_T_ref` e `func_Q_ref` são handles criados por `create_thrust_model` e `create_torque_model`, que mapeiam PWM → empuxo (N) e PWM → torque (N·m) usando polinômios de grau 3 ajustados a dados de bancada.
