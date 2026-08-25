# Plano de Reorganização da Dissertação
*(comentário da banca sobre organização + notas da orientadora + alinhamento com as pendências)*
*Atualizado em 25/08/2026*

---

## 1. O diagnóstico da banca

> Estrutura foge do formato tradicional (Introdução, Fundamentação/Revisão, Metodologia,
> Resultados e Discussão, Conclusões). Especificamente:
> - **Cap. 4** mistura revisão bibliográfica, metodologia e critérios de validação;
> - **Cap. 5** mistura resultados de identificação com o desenvolvimento do controlador
>   (que é uma segunda contribuição);
> - falta separação clara entre **"como foi feito"** e **"o que foi obtido"**.

## 2. Estrutura ATUAL (arquivo → capítulo → conteúdo)

| Cap. | Arquivo | Conteúdo |
|---|---|---|
| 1 Introdução | `Cap1/cap1.tex` | Motivação; **revisão dos 3 VTOLs (azul, feita)**; objetivos; organização |
| 2 Modelagem Matemática | `Cap2/cap2.tex` | Coordenadas; transformações; propulsão; alocação; **composição F/M (azul)**; **aero da estrutura (azul)**; rotacional; translacional; acelerômetro; Euler; modelo não linear |
| 3 Plataforma | `Plataforma/plataforma.tex` | Plataforma; propulsão (motor/ESC); aviônica; **Campanha de voo** (`sec:campanha`) |
| 4 Identificação | `Cap3/cap3.tex` | 4.1 sistemas instáveis (caixas branca/cinza/preta; malha fechada×aberta — Figs 4.1/4.2 refeitas); 4.2 M4V: manobras, medições, OEM (Fig 4.10 refeita), modelos, critérios de validação |
| 5 Resultados | `Cap5/cap5.tex` | 5.1 coerência dos dados; 5.2 identificação OEM; 5.3 validação; 5.4 **influência aero (azul, nova)**; 5.5–5.7 controle (`\input{Cap4/cap4.tex}`: arquitetura, projeto PID, simulação linear, limitações); 5.8 comparação linear×NL |
| 6 Conclusões | `Cap6/Cap6.tex` | Conclusão; trabalhos futuros |
| Ap. A | `ApeA/apendiceA.tex` | Bancada do motor |
| Ap. B | `ApeB/apendiceB.tex` | **Linearização (JÁ é apêndice — ver §6.4)** |

## 3. Estrutura ALVO (notas da orientadora, interpretadas)

| Novo cap. | Título provisório | Conteúdo |
|---|---|---|
| 1 | Introdução | **Funil de 3 estágios**: (i) VTOLs para uso civil *(refs novas)* → (ii) controle de VTOL e a necessidade do modelo identificado *(refs novas)* → (iii) identificação aplicada (Græsdal/Salahudden/Nguyen — **pronto em azul**) + onde o trabalho inova. Mais referências no geral |
| 2 | Plataforma e Metodologia de Identificação ⚠️ *(nome/escopo "ainda não está certo")* | Plataforma **sem** a campanha de voo; formulação do problema (foco caixa cinza); panorama metodológico (malha fechada×aberta, EEM→OEM); ponte: "com o modelo, desenvolvem-se as malhas de controle". **⬅ destino natural da revisão metodológica crítica do comentário 1** |
| 3 | Modelagem Matemática | Atual Cap. 2 + apresentação dos coeficientes + fecho novo: **a necessidade da identificação** |
| 4 | Identificação | Atual Cap. 4 com correções; **Campanha de voo entra na "Aquisição de dados de saída"**; prosa das manobras/persistência/espectro (comentário 3) |
| 5 | Resultados | 5.1 coerência (**⚠️ verificar se permanece**); identificação; validação; influência aero; controle; comparação; **+ SEÇÃO DE DISCUSSÃO** |
| 6 | Conclusões | inalterado |

## 4. Matriz de movimentação

| Bloco | De → Para | Azul? |
|---|---|---|
| Plataforma (sem campanha) | Cap. 3 atual → Cap. 2 novo | não (movido) |
| Campanha de voo (`sec:campanha`) | Plataforma → Cap. 4 "Aquisição de dados de saída" | não (movido) |
| §4.1 (caixas, malha fechada×aberta, Figs 4.1/4.2) | Cap. 4 → Cap. 2 novo (metodologia) ⚠️ decidir quanto sobe | não |
| Modelagem completa | Cap. 2 → Cap. 3 novo | não |
| M4V detalhado (manobras/medições/OEM/modelos/critérios) | permanece no Cap. 4 | não |
| Projeto+simulação de controle (`Cap4/cap4.tex`) | dentro de Resultados → **decidir**: cap. próprio? (§6.2) | não |
| Linearização | já é Apêndice B — nada a mover (ver §6.4) | — |

## 5. Alinhamento com as PENDÊNCIAS (escrever uma vez, já no lugar certo)

| Pendência (comentário) | Destino na estrutura NOVA | Status |
|---|---|---|
| Revisão metodológica crítica: EEM/OEM/FEM, Kalman, frequência (CIFER/Wei/Ivler), ML (Bauersfeld), black-box; parâmetros típicos; vantagens/limitações (**com. 1**) | Cap. 2 novo (metodologia) | a escrever; Wei2017/Ivler2019 a entrar no bib |
| Funil da introdução + mais referências (**com. 1 / orientadora**) | Cap. 1 | a pesquisar (VTOL civil; controle de VTOL) |
| Revisão dos 3 VTOLs | Cap. 1, estágio (iii) do funil | ✅ feita (azul) |
| Prosa de manobras, Δt a priori, persistência, espectro (**com. 3**) | Cap. 4 (Manobras/Aquisição) | redigida em chat; falta inserir |
| Campanha de voo | Cap. 4 Aquisição | mover |
| "Necessidade da identificação" | fecho do Cap. 3 novo | a escrever (curto) |
| Robustez Monte Carlo (**com. 13**) | Resultados/Discussão (ou cap. de controle) + ajustar "Limitações" | ✅ rodado; texto proposto aguardando ok |
| Justificativa ωn/ζ + Limitações (**com. 13**) | acompanham o bloco de controle | ✅ feitas + regeneradas |
| Influência aero + separabilidade (**com. 2**) | Resultados/Discussão | ✅ feita (azul) |
| Apêndice identificação asa fixa | novo Apêndice C | a escrever |
| Discussão consolidada (orientadora) | seção nova em Resultados (ou cap. próprio) | a escrever |
| Siglas (CR, AVL, COCTA, GPS; duplicatas OEM/EEM/ESC/IMU/CG; M4V invertida) | independente | fazer APÓS a reorganização |
| Micro-edits Motivação/Objetivos | Cap. 1 | aprovados, não aplicados |

## 6. Decisões em aberto (levar à orientadora)

1. **Escopo do Cap. 2 novo**: quanto do atual §4.1 sobe para a "metodologia geral" e quanto
   fica no Cap. 4 detalhado? Proposta: Cap. 2 = problema + caixas + malha aberta×fechada +
   taxonomia crítica (revisão do com. 1) + visão geral M4V/OEM; Cap. 4 = execução detalhada.
2. **Controle: capítulo próprio?** A banca chama o controlador de "segunda contribuição".
   Um "Cap. 6 — Projeto e Avaliação do Controlador" (projeto, simulação, comparação,
   robustez MC, limitações) resolveria a crítica "resultados misturados" de frente.
   Alternativa mínima: manter em Resultados com subseções bem separadas + Discussão.
3. **§5.1 (coerência dos dados)**: é verificação de METODOLOGIA (qualidade da campanha),
   não resultado científico. Opções: mover para o fim do Cap. 4 (fecho da aquisição) ou
   manter como primeira seção dos Resultados. Tendência: Cap. 4.
4. **"Linearização → apêndice"**: a dedução JÁ está no Apêndice B. O que resta no corpo é
   o *projeto sobre o modelo linear* e a *simulação linear* — confirmar se a nota se refere
   a algo além disso.
5. **Nome do Cap. 2**: "Plataforma e Metodologia de Identificação" vs separar em dois
   capítulos curtos ("Fundamentação/Revisão Metodológica" + "Plataforma").

## 7. Riscos e cuidados

- **Azul**: trechos movidos NÃO recebem azul (convenção com o orientador); só texto novo.
- **23 referências textuais** a "Capítulo~N" no corpo — revisar todas após renumerar
  (labels simbólicos se ajustam sozinhos; menções por extenso não).
- Dependência cruzada: Plataforma cita `eq:alocacao` (da Modelagem). Na ordem nova
  (Plataforma ANTES da Modelagem) essa referência fica *para frente* — reescrever a
  passagem da matriz B₀ numérica (movê-la para o Cap. 3 novo ou citar "apresentada no
  Capítulo 3").
- Recompilar e comparar sumário a cada fase; conferir Figuras/Tabelas renumeradas nos
  textos que as citam por número em prosa.
- As figuras TikZ refeitas (4.1, 4.2, 4.3, 4.4, 4.10, 5.7, 2.6, 2.7) migram sem retrabalho.

## 8. Ordem de execução proposta

- **F0 — Decisões** (§6) com a orientadora. *(bloqueia F1)*
- **F1 — Mecânica**: reordenar `\input` no `tese.tex`; mover campanha p/ Cap. 4; separar
  controle (se cap. próprio); consertar `eq:alocacao` na Plataforma; varrer os 23
  "Capítulo~N"; recompilar e conferir sumário. *Sem texto novo.*
- **F2 — Textos novos nos lugares certos**: funil do Cap. 1; revisão metodológica (Cap. 2);
  necessidade da identificação (Cap. 3); prosa das manobras (Cap. 4); Discussão; robustez MC.
- **F3 — Apêndice asa fixa; siglas; varredura final** (números regenerados, refs, sumário azul).
