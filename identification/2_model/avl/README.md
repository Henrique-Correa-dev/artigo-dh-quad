# Modelo aerodinâmico do DH em AVL e XFLR5

## Para que serve (e para que não serve)

**Não serve** para obter momentos de inércia. Nem o AVL nem o XFLR5 calculam inércia
a partir da geometria: nos dois o usuário fornece as massas concentradas e o
programa apenas soma pelo teorema dos eixos paralelos, exatamente o que
`2_model/inertia_lumped.m` já faz. Para o a priori de inércia as saídas reais são
completar o CAD do SolidWorks (que hoje pesa 1,072 kg contra 1,993 kg reais) ou
medir com pêndulo bifilar.

**Serve** para as derivadas aerodinâmicas de asa fixa (Cl_p, Cm_q, Cn_r, Cm_α,
Cn_β, Cl_β). No modo multirrotor essas derivadas não aparecem sozinhas, mas
entram agregadas em `cp`, `cq`, `cr` conforme a velocidade de avanço, e por isso
delimitam o envelope de validade do modelo identificado (ver `prior_damping.m`).

## Estado do trabalho anterior

O relatório FINEP (`infos_Rela_FINEP_2.pdf`, seção 4.1) registra que o modelo AVL
do DH **não convergiu**: "durante as simulações o software não apresentou bons
resultados, devido a erros de convergência das derivadas de estabilidade". Os
valores da coluna "AVL / valor inicial" nas tabelas 4.0, 4.2 e 4.7 daquele
relatório **não são do DH**: são do *Baby Shark* de Græsdal (2021), adotado como
alternativa. O XFLR5, esse sim, rodou para o DH (Cm_α = −0,783, margem estática
positiva, short-period 1,078 Hz, fugoide 0,171 Hz, dutch-roll 0,542 Hz).

Ou seja, reproduzir "o AVL dele" significa refazer algo que ficou pendente. Os
arquivos deste diretório são uma reconstrução para fechar essa lacuna.

## Arquivos

| arquivo | conteúdo |
|---|---|
| `dh.avl` | geometria (asa, empenagem em H, superfícies de comando) |
| `dh.mass` | massas concentradas, para a análise de modos |
| `usa35b.dat` | **falta**: baixar em airfoiltools.com (perfil USA-35B) |

## Como rodar o AVL

1. Instalar (macOS): baixar o binário em <https://web.mit.edu/drela/Public/web/avl/>
   e colocar em `/usr/local/bin/avl`.
2. Baixar o perfil USA-35B em airfoiltools.com como `usa35b.dat` neste diretório.
   Sem ele, trocar as duas linhas `AFIL / usa35b.dat` por `NACA / 4412`, que é uma
   aproximação razoável do arqueamento do USA-35B.
3. Rodar:

```
avl dh.avl
```

4. Dentro do AVL:

```
MASS dh.mass        carrega as massas
MSET 0              aplica a todos os casos
OPER                menu de operação
   A A 5            ângulo de ataque = 5 graus
   X                executa
   ST               derivadas de estabilidade (é a tela que interessa)
   W dh_st.txt      salva
   Q
MODE                análise de modos (precisa do .mass)
   N                calcula autovalores
```

`ST` imprime Cl_p, Cm_q, Cn_r, Cn_β, Cl_β em forma adimensional, prontos para
comparar com as estimativas de `prior_damping.m`.

## Por que o AVL do relatório provavelmente falhou

O AVL é um método de malha de vórtices linear: ele não "deixa de convergir" no
cálculo aerodinâmico. Quando o relatório fala em erro de convergência, quase
sempre a causa está na análise de **modos**, que precisa do arquivo de massa. As
causas comuns, em ordem de frequência:

1. Arquivo `.mass` com massa total ou inércia inconsistente (é o caso aqui: o CAD
   dava 1,072 kg, e o peso real é 1,993 kg).
2. CG (`Xref`) fora da faixa entre o bordo de ataque e o ponto neutro, o que gera
   modo instável de raiz muito grande e trava a interpretação.
3. `YDUPLICATE` declarado em superfície que já foi definida nos dois lados,
   duplicando a área e dobrando as derivadas.
4. Empenagem em H modelada como uma única superfície vertical no plano de simetria
   (o DH tem duas derivas fora do plano, com `YDUPLICATE`).

## XFLR5

O XFLR5 é mais simples de reproduzir e já funcionou para o DH:

1. `Plane Design → Define a New Plane`
2. Asa: envergadura 1,20 m, corda de raiz 0,25 m, corda de ponta 0,20 m, perfil
   USA-35B (importar o `.dat` em `Direct Foil Design` antes).
3. Estabilizador horizontal: 0,45 m de envergadura, corda 0,17 m, a 0,52 m do
   bordo de ataque da asa, perfil simétrico NACA 0009.
4. Duas derivas de 0,12 m nas pontas do estabilizador.
5. `Define Inertia`: inserir as massas concentradas (as mesmas de `dh.mass`).
6. Análise: `Stability Analysis`, tipo 7, velocidade de cruzeiro 15 m/s.
7. A aba de resultados traz as derivadas e os autovalores dos modos.

Cuidado com o número de Reynolds: a 15 m/s e corda 0,226 m dá Re ≈ 230 000, então
a polar do perfil precisa ser calculada nessa faixa, não nos 506 000 usados no
artigo (que corresponde a 35 m/s).
