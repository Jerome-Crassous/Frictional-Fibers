# FrictionalFibers — Documentation (Fiber_v9.6)

**Jérôme Crassous**  
Univ Rennes, CNRS, IPR (Institut de Physique de Rennes) – UMR 6251, F-35000 Rennes, France  
PMMH, CNRS, ESPCI Paris, Université PSL, Sorbonne Université, Université de Paris, 75005 Paris, France  
jerome.crassous@univ-rennes1.fr — [https://jerome-crassous.github.io/](https://jerome-crassous.github.io/)

> *Cette documentation a été rédigée avec l'aide d'une intelligence artificielle (Claude, Anthropic), à partir du code source et de l'article de référence.*

<table align="center">
  <tr>
    <td align="center"><img src="images/curved_fabric.png" height="160" alt="Tissu tricoté courbé"><br><sub>Tissu tricoté courbé</sub></td>
    <td align="center"><img src="images/bowline003.png" height="160" alt="Nœud de chaise (bowline)"><br><sub>Nœud de chaise (<i>bowline</i>)</sub></td>
    <td align="center"><img src="images/multiple_coils.png" height="160" alt="Assemblée de fibres hélicoïdales"><br><sub>Assemblée de fibres hélicoïdales</sub></td>
  </tr>
  <tr>
    <td align="center"><img src="images/cordage.png" height="160" alt="Filet tressé impacté par une sphère"><br><sub>Filet tressé impacté par une sphère</sub></td>
    <td align="center"><img src="images/bobine.png" height="160" alt="Enroulement sur une bobine"><br><sub>Enroulement sur une bobine</sub></td>
    <td align="center"><img src="images/hair.png" height="160" alt="Fibres enchevêtrées (cheveux)"><br><sub>Fibres enchevêtrées (cheveux)</sub></td>
  </tr>
</table>

Simulation par éléments discrets (DEM) de tiges élastiques (fibres) en contact frictionnel, sur GPU (OpenCL). Le modèle mécanique et sa validation sont publiés dans :

> J. Crassous, *Discrete-element-method model for frictional fibers*, Phys. Rev. E **107**, 025003 (2023). DOI: 10.1103/PhysRevE.107.025003

Cette documentation ne redémontre pas ce modèle en détail (voir la Partie III et l'article ci-dessus) ; elle explique comment **utiliser** le programme (Partie I), comment il est **implémenté en parallèle** (Partie II), et rappelle en annexe le **modèle physique** sous-jacent (Partie III).

## Sommaire

- [Partie I — Guide utilisateur](#partie-i--guide-utilisateur)
- [Partie II — Algorithme et kernels](#partie-ii--algorithme-et-kernels)
- [Partie III — Annexe : modèle physique](#partie-iii--annexe--modèle-physique)

---

# Partie I — Guide utilisateur

Cette partie s'adresse à qui veut **utiliser** le programme sur son propre cas, en s'appuyant sur les exemples fournis, sans nécessairement comprendre le détail des kernels GPU.

## I.1 Vue d'ensemble

Le programme simule un ensemble de fibres élastiques discrètes (chaînes de segments cylindriques) soumises à des forces internes (étirement, flexion, torsion) et à des forces de contact frictionnelles entre elles. Le cœur du calcul (forces, contacts, intégration temporelle) s'exécute sur GPU via des kernels OpenCL ; le CPU ("host") ne fait que préparer la géométrie initiale, piloter la boucle, et périodiquement relire l'état pour le sauvegarder.

Chaque **cas d'étude** (une géométrie particulière : poutre encastrée, hélice, tissu, filet de tennis, etc.) est un petit projet C autonome qui réutilise une bibliothèque commune (`Librairie_v9.6`) et un jeu de kernels communs, et qui ne fournit que :

- la géométrie et les paramètres initiaux (fonctions host `SetFibersParameters` / `SetFibersInitialPositions`),
- trois kernels "spécifiques" (`kernel_specific.cl`) qui appliquent les forces externes et les contraintes géométriques propres au cas étudié (encastrement, traction imposée, torsion imposée, etc.).

Tout le reste (calcul des forces élastiques, détection et traitement des contacts, intégration de Verlet, transport du repère matériel) est **générique** et partagé par tous les exemples.

## I.2 Organisation des fichiers sources

```
sources/
├── fiberLib_Common_Macros_v9.6.h      # constantes/macros partagées (WG, statuts, ...)
├── fiberLib_OpenCL_v9.6.h             # structures host (fiber, parameter, contact) + prototypes
├── functions_v9.6.cl                  # petites fonctions OpenCL utilitaires (rotation, transport parallèle, ...)
├── kernel_le_force_v9.6.cl            # forces élastiques : longueurs/e_i, étirement, flexion, torsion
├── kernel_poss_contact_v9.6.cl        # détection des contacts potentiels (broad-phase)
├── kernel_contact_v9.6.cl             # calcul géométrique + forces de contact (narrow-phase)
├── kernel_integrate_shift_v9.6.cl     # intégration de Verlet + transport du repère matériel
├── kernel_max_v9.6.cl                 # réductions (déplacement max, comptage de contacts, dissipation, ...)
│
├── Librairie_v9.6/                    # bibliothèque host (C) partagée par tous les exemples
│   ├── fiber.c                        # allocation d'une fibre (AllocateOneFiberLib)
│   ├── geometry.c                     # repère matériel initial, courbure, longueurs (host)
│   ├── energy.c                       # énergies (host, pour diagnostic/plots)
│   ├── math.c                         # algèbre vectorielle 3D, rotations, RNG
│   ├── read_save_config.c             # sauvegarde/lecture de configurations (config.dat, reduced_Config.tmp)
│   ├── device_init.c                  # initialisation OpenCL, compilation des kernels, Execute_Kernel_*
│   ├── device_host_device.c           # allocation des buffers, transferts host <-> device
│   ├── device_err_code.c, times.c, print.c
│
├── Visu_v9.6/                         # petit visualiseur OpenGL/GLUT + export POV-Ray
│
└── Fiber_v9.6_<NomDuCas>/              # un dossier par cas d'étude, ex. SimpleBending, Twist, Coil, ...
    ├── main.c / main.h                # point d'entrée, paramètres globaux par défaut
    ├── core.c                         # DoOneIteration() : orchestre les kernels à chaque pas de temps
    ├── <nom_du_cas>.c                 # SetFibersParameters() + SetFibersInitialPositions() (spécifique)
    └── kernel_specific.cl             # kernel_Specific_OneTime / _Force / _Position (spécifique)
```

Les exemples fournis dans ce dépôt sont : `SimpleBending` (flexion d'une poutre encastrée), `Twist` (flambage en torsion d'une tige), `Coil` (fibre unique enroulée en hélice), `Multiple_coils` (assemblée de fibres hélicoïdales, type mèche/fil), `Jersey` (tricot, mailles jersey sur une grille), et `Tennis_v2_square` (filet tressé, type raquette de tennis).

## I.3 Compilation et prérequis

Le code est écrit en C, cible Windows (`windows.h`, `omp.h`) et utilise l'API **OpenCL** (`CL/cl.h`) pour piloter le GPU. Chaque cas d'étude est un exécutable indépendant qui doit être compilé avec :

- accès au SDK OpenCL (headers + lib de votre fournisseur GPU, ou un ICD générique),
- les fichiers de `Librairie_v9.6` compilés/liés au projet,
- exécution **depuis le dossier du cas d'étude** : `Create_Kernels()` charge les sources `.cl` par chemin relatif (`"..\functions_v9.6.cl"`, `"kernel_specific.cl"`, etc.), donc le répertoire de travail courant doit être ce dossier.

À l'exécution, le programme choisit une plateforme OpenCL (il préfère une plateforme dont le nom commence par "NV", sinon la première trouvée), crée un contexte et une file de commandes, puis compile tous les `.cl` en un seul programme.

## I.4 Discrétisation d'une tige élastique : notions de base

Une fibre est représentée par une **chaîne de N nœuds** $r_0, \dots, r_{N-1}$ reliés par **N−1 segments cylindriques** rectilignes (le segment $i$ relie $r_i$ à $r_{i+1}$). C'est la discrétisation "tiges élastiques discrètes" (*discrete elastic rods*, Bergou et al.) reprise et adaptée au contact frictionnel par l'article de référence.

Pour chaque segment $i$ on définit :

- le **vecteur directeur unitaire** $e^i = (r_{i+1}-r_i)/\lVert r_{i+1}-r_i\rVert$,
- la **longueur courante** $l^i = \lVert r_{i+1}-r_i\rVert$, censée rester proche de la longueur au repos $l_0$ (le calcul de flexion/torsion suppose $l^i \simeq l_0$, l'étirement réel étant traité séparément par un ressort longitudinal),
- deux **vecteurs matériaux unitaires** $m^i_{(1)}, m^i_{(2)}$ tels que $(e^i, m^i_{(1)}, m^i_{(2)})$ forme un trièdre direct. Ils portent l'orientation "torsion" du segment (comme les repères matériels d'une poutre de Kirchhoff).
- un **angle de torsion** $\theta_i$ : $m^i_{(1)}$ est obtenu par rotation d'un vecteur de référence sans torsion $\bar m^i_{(1)}$ ("m1_bar") d'un angle $\theta_i$ autour de $e^i$. $\bar m^i_{(1)}$ est mis à jour à chaque pas de temps par **transport parallèle** (sans torsion ajoutée) le long de la courbe : c'est ce qui permet de séparer proprement la courbure (flexion) de la torsion.

L'état cinématique complet d'une fibre à un instant donné est donc entièrement décrit par :

- les positions des nœuds $r_i$, $0 \le i \le N-1$,
- les angles de torsion $\theta_i$, $0 \le i \le N-2$ (un angle par segment).

Chaque nœud/segment porte aussi une **courbure naturelle** ("au repos") $\bar\kappa_{(1),i}, \bar\kappa_{(2),i}$ : c'est la courbure que la fibre aurait spontanément sans contrainte extérieure — nulle pour une fibre droite, non nulle si l'on veut modéliser une fibre naturellement bouclée/frisée. C'est un champ (`kappa1_bar[i]`, `kappa2_bar[i]`) fourni par l'utilisateur segment par segment.

Enfin, chaque **fibre** (et non chaque segment) porte des propriétés matériaux uniformes : rayon `radius`, longueur de segment au repos `l0`, raideur d'étirement `k0`, masse par nœud `masse`, moment d'inertie de torsion `j`, module de flexion `bending` ($B$) et module de torsion `c` ($C$).

## I.5 Les structures de données côté utilisateur

Deux structures C (définies dans `fiberLib_OpenCL_v9.6.h`) sont manipulées côté host pour décrire une simulation.

### La structure `fiber` (une par fibre)

| Champ | Rôle |
|---|---|
| `n` | nombre de nœuds N de la fibre (donc N−1 segments) |
| `status` | `STATUS_FREE` (fibre dynamique), `STATUS_FIXED` (immobile), `STATUS_VIRTUAL` (fibre "bouche-trou", voir I.7) |
| `radius`, `l0`, `k0`, `masse`, `j` | rayon du cylindre, longueur de segment au repos, raideur d'étirement, masse par nœud, moment d'inertie |
| `bending` ($B$), `c` ($C$) | module de flexion et module de torsion |
| `xt[i][3]`, `xtm[i][3]` | position du nœud $i$ à l'instant $t$ et à l'instant précédent $t-\Delta t$ (nécessaire à l'intégrateur de Verlet) |
| `thetat[i][3]`, `thetatm[i][3]` | seule la composante `[2]` (l'axe $e^i$, stocké comme composante "z" locale) sert : c'est $\theta_i$ à $t$ et $t-\Delta t$ |
| `e[i][3]` | vecteur directeur $e^i$ du segment $i$ (calculé, pas à fournir) |
| `m1_bar[i][3]`, `m1[i][3]` | repère matériel de référence (transporté sans torsion) et repère matériel réel |
| `kappa1_bar[i]`, `kappa2_bar[i]` | courbure naturelle du segment $i$ |
| `f[i][3]`, `moment[i][3]` | force et moment accumulés sur le nœud/segment (calculés par les kernels) |
| `flag[i]` | entier libre, utilisable par le kernel spécifique pour marquer certains nœuds (ex. bords d'un tissu) |

### La structure `parameter` (paramètres globaux de la simulation)

Les champs les plus utilisés en pratique :

| Champ | Rôle |
|---|---|
| `dt` | pas de temps $\Delta t$ |
| `kn`, `kt` | raideurs normale et tangentielle de contact |
| `lambda` | amortissement visqueux global (transverse), $f^{visc} = -\lambda \dot r$ |
| `lambda_internal` | amortissement visqueux de l'étirement (longitudinal) |
| `lambda_contact_n`, `lambda_contact_t` | amortissement visqueux du contact (normal/tangentiel) |
| `mu` | coefficient de frottement de Coulomb |
| `nFiber`, `nSegment` | nombre de fibres, nombre total de segments (calculé) |
| `nContactMax` | taille de la liste de contacts (doit être un multiple de `WG`, voir Partie II) |
| `epsStar` | seuil de distance pour la recherche de contacts potentiels (auto-ajusté, voir II.5) |
| `periodic` | `PERIODIC_NO`, `PERIODIC_XY`, ... conditions aux limites périodiques |
| `bending`, `c` | valeurs par défaut de $B$ et $C$ (souvent recopiées dans chaque `fiber`) |
| `iter`, `iterStop` | compteur d'itérations et itération d'arrêt |
| `lx, ly, lz`, `R`, `RHelix`, `pitch`, ... | paramètres géométriques libres, utilisés à la convenance de chaque cas d'étude (boîte périodique, rayon d'un cylindre d'enroulement, pas d'une hélice, etc.) |

Les constantes utiles sont dans `fiberLib_Common_Macros_v9.6.h`, notamment :

- `WG = 1024` : taille de work-group utilisée par la plupart des kernels 1D (voir Partie II),
- `STATUS_FREE=0`, `STATUS_FIXED=1`, `STATUS_VIRTUAL=2`,
- `NOCONTACT`, `OLDCONTACT`, `STILLCONTACT`, `NEWCONTACT` : états d'un contact.

## I.6 Cycle de vie d'une simulation (`main.c`)

Le `main()` d'un cas d'étude fait toujours, dans l'ordre :

1. Allouer `gParameterPtr` et fixer les paramètres par défaut (dt, kn, kt, lambda, bending, mu, nContactMax, ...).
2. Allouer `gFiber` (tableau de `fiber`) et `gContact`.
3. Appeler les deux fonctions **spécifiques au cas d'étude** :
   - `SetFibersParameters()` : alloue chaque fibre (`AllocateOneFiberLib`) et fixe ses propriétés matériaux ;
   - `SetFibersInitialPositions()` : fixe la géométrie initiale $r_i$, $\theta_i$, $\bar\kappa_{(1,2),i}$, et initialise le repère matériel via `Calculate_Relaxed_M1_Bar` puis `Calculate_M1_From_M1Bar_And_Theta`.
4. Initialiser le device : `InitDevice()`, `Create_Device_Ptrs()`, `Create_Host_Ptrs()`, `Create_Kernels()` (compile les `.cl` et fixe les arguments de chaque kernel une fois pour toutes).
5. Copier l'état initial host → device (`Copy_Utils_host2dev`, `Copy_Fiber_host2dev`, `Copy_Param_host2dev`).
6. Boucler indéfiniment sur `DoOneIteration()` (définie dans `core.c` du cas d'étude, identique dans son squelette d'un exemple à l'autre — voir Partie II.3), qui périodiquement rapatrie l'état sur le host (`Copy_Fiber_dev2host`) pour l'écrire sur disque (`SaveReducedConfigLib`) et afficher un diagnostic (ex. énergie de flexion).

## I.7 Écrire sa propre géométrie : les 4 fonctions à fournir

Pour créer un nouveau cas d'étude, il faut fournir :

**1) `SetFibersParameters(void)`** — alloue chaque fibre réelle avec `AllocateOneFiberLib(gFiber, iFiber, n)` et fixe `status`, `radius`, `l0`, `k0`, `masse`, `bending`, `c`, `j`.

**2) `SetFibersInitialPositions(void)`** — fixe `xt[i]`, `thetat[i][2]`, `kappa1_bar[i]`, `kappa2_bar[i]`, initialise `m1_bar[0]` avec un vecteur non parallèle à $e^0$ puis appelle `Calculate_Relaxed_M1_Bar` (propage le repère sans torsion le long de la fibre par transport de Darboux) et `Calculate_M1_From_M1Bar_And_Theta` (applique la rotation $\theta_i$). Termine en recopiant `xt → xtm` et `thetat → thetatm` (état "au repos", vitesse initiale nulle pour l'intégrateur de Verlet).

**Fibre virtuelle.** Les kernels 1D tournent avec une taille de work-group fixe `WG=1024` ; le nombre total de segments doit donc être un multiple de `WG`. Tous les exemples ajoutent en conséquence, après les fibres "réelles", une dernière fibre de statut `STATUS_VIRTUAL` dont la seule fonction est de compléter `nSegment` au multiple de `WG` immédiatement supérieur :

```c
missingSegment = ((gParameterPtr->nSegment + WG - 1) / WG) * WG - gParameterPtr->nSegment;
iFiber = gParameterPtr->nFiber - 1;
AllocateOneFiberLib(gFiber, iFiber, missingSegment);
gFiber[iFiber].status = STATUS_VIRTUAL;
```
Les segments virtuels sont ignorés par les kernels de force et de contact (tests `if (status[iFiber]==STATUS_VIRTUAL) return;`).

**3) Trois kernels "spécifiques" dans `kernel_specific.cl`** — c'est le point d'extension principal, appelé à chaque itération :

- `kernel_Specific_OneTime(param_Uint, param_Double)` : exécuté **une seule fois** par itération (taille de grille = 1), utile pour mettre à jour un compteur ou un paramètre global (ex. incrémenter `param_Uint[0]`, faire évoluer une force cible dans le temps).
- `kernel_Specific_Force(...)` : exécuté **une fois par segment**, après le calcul des forces élastiques ; il ajoute les **forces externes** propres au cas étudié (traction imposée en bout de fibre, gravité, etc.). Dans l'exemple `SimpleBending` :
  ```c
  __kernel void kernel_Specific_Force(...) {
      iSegment = (uint)(get_global_id(0));
      iFiber = iFiberFromSegment[iSegment];
      i = iFromSegment[iSegment];
      if ((iFiber == 0) && (i + 1 == n[iFiber]))
          f[3 * iSegment + 2] -= 1.e-4;   // force ponctuelle au bout de la fibre 0
  }
  ```
- `kernel_Specific_Position(...)` : exécuté **une fois par segment**, après l'intégration temporelle ; il impose les **contraintes géométriques** (encastrement, déplacement ou torsion imposés) en écrasant `xt`/`thetat` par une valeur prescrite. Dans `SimpleBending`, les deux premiers nœuds de la fibre sont ramenés à leur position précédente, ce qui réalise l'encastrement :
  ```c
  if ((iFiber == 0) && (i <= 1)) {
      for (k = 0; k < 3; k++) {
          xt[3*iSegment+k] = xtm[3*iSegment+k];
          thetat[3*iSegment+k] = thetatm[3*iSegment+k];
      }
  }
  ```
  Le fichier `functions_v9.6.cl` fournit une fonction utilitaire `FastCheck(fiber1,i1,fiber2,i2)` qu'on peut redéfinir dans `kernel_specific.cl` pour exclure *a priori* certaines paires fibre/nœud de la détection de contact (retourner `SKIP_CONTACT_DETECTION`), ce qui évite de tester des contacts nécessairement absents et accélère la simulation.

**4) `main.c`/`main.h`** — copiés depuis un exemple existant et adaptés : nombre de fibres, `nContactMax`, valeurs par défaut de `dt`, `kn`, `kt`, `bending`, `mu`, etc.

## I.8 Exemple minimal commenté : `SimpleBending`

C'est l'exemple le plus simple : la flexion statique d'une poutre encastrée sous une force ponctuelle en bout, comparable au premier test de validation de l'article (déflexion d'une tige encastrée, Fig. 5).

**Paramètres** (`main.c`) : `nFiber = 2` (une fibre réelle + une fibre virtuelle), `ls = 1.5` (longueur de segment $l_0$), `bending = 0.1` ($B$), `nContactMax = 32·WG`, `dt = 0.1`, `kn = 1`, `kt = 0.5`, `lambda_internal = 2.8`, `mu = 0.5`.

**Géométrie** (`simple_bending.c`) : une fibre de 20 segments, alignée sur l'axe $x$ :
```c
AllocateOneFiberLib(gFiber, iFiber, 20);
gFiber[iFiber].status = STATUS_FREE;
gFiber[iFiber].l0 = gParameterPtr->ls;
...
gFiber[iFiber].xt[i][0] = ((double)i) * gFiber[iFiber].l0;   // ligne droite selon x
gFiber[iFiber].thetat[i][1] = 0.;                            // pas de torsion
gFiber[iFiber].kappa1_bar[i] = 0.;                           // pas de courbure naturelle
```

**Kernels spécifiques** (`kernel_specific.cl`, voir I.7) : une force `f_z -= 1e-4` est appliquée en permanence au dernier nœud de la fibre, et les deux premiers nœuds sont gelés à leur position précédente à chaque pas de temps — ce qui réalise l'encastrement.

**Déroulement** (`core.c`, fonction `DoOneIteration`) : à chaque itération, forces élastiques → force spécifique → détection/calcul des contacts (ici sans effet, la fibre ne se touche pas elle-même) → intégration → contrainte géométrique spécifique → mesure du déplacement max. Toutes les `10 000` itérations, l'état est rapatrié sur le host, sauvegardé, et la hauteur du bout de la poutre ainsi que l'énergie de flexion sont affichées :
```c
printf("iter = %u h = %e enrgy = %e\n", iter, gFiber[0].xt[gFiber[0].n-1][2], eBending);
```
À `iter == iterStop`, le profil final de la poutre est écrit dans `profile.txt` (colonnes $x_i, z_i$) et le programme s'arrête.

C'est le squelette à reproduire pour n'importe quel nouveau cas : seuls la géométrie initiale et le contenu des trois kernels spécifiques changent.

## I.9 Tour des autres exemples fournis

- **`Twist`** — une fibre unique alignée selon $-z$, avec un repère matériel initial `m1_bar[0] = (0,1,0)`. Utilisé pour étudier le flambage en torsion d'une tige comprimée/tordue à ses extrémités (comparable au test de flambage de l'article, Fig. 6) : les contraintes aux extrémités (couple imposé, blocage transverse) sont réalisées dans `kernel_specific.cl` (non détaillé ici, à adapter selon le couple cible).

- **`Coil`** — une fibre unique positionnée initialement sur une **hélice paramétrique** :
  ```c
  x = RHelix*sin(t);  y = RHelix*cos(t);  z = pitch*t;
  ```
  (fonction `S2PositionHelix`, paramètres `gParameterPtr->RHelix`, `gParameterPtr->pitch`), les nœuds étant placés à abscisse curviligne constante $l_0$ le long de cette courbe. La courbure naturelle `kappa1_bar`/`kappa2_bar` est calculée nœud par nœud à partir de la variation de $e^i$ le long de l'hélice, de sorte que l'hélice soit une configuration **au repos** (sans contrainte) de la fibre — typiquement pour étudier une mèche/ressort déjà enroulé(e).

- **`Multiple_coils`** — même principe géométrique que `Coil`, mais appliqué à *plusieurs* fibres (`SetOneFiberInitialPositions` appelée pour chaque fibre), chacune décalée aléatoirement (`shift` tiré avec `MyRand()`) dans une boîte `lx × ly × lz` : utile pour étudier un paquet de fibres/ressorts en interaction frictionnelle.

- **`Jersey`** — modélise une **maille tricotée** (jersey) : la forme d'une boucle de tricot est lue depuis un fichier de points externe (`x_theta_lw70_lc70.txt`), rééchantillonnée (`loop_raw` → `loop`, interpolation linéaire), puis reproduite fibre par fibre sur une grille `GN_CELL_X × GN_CELL_Y` de cellules (fibres traversantes horizontales, fibres basses/hautes, coins), avec un traitement particulier des fibres de bord. C'est l'exemple le plus proche d'une géométrie "industrielle" complexe construite par composition de blocs.

- **`Tennis_v2_square`** — un filet croisé (façon cordage de raquette) : une moitié des fibres suit une position `S2Position` en cosinus (`amp*cos(2π t /(2L))`) puis est symétrisée et répétée avec alternance de signe (`pow(-1, iFiber-nHalfHalf)`) pour créer un tressage, l'autre moitié des fibres étant droites et perpendiculaires. Un fichier `test_energy.c` additionnel permet de vérifier numériquement les énergies calculées par les kernels contre les formules host (`BendingEnergyLib`, etc.).

Dans tous les cas, comparez le fichier `<nom>.c` du cas qui vous intéresse à celui de `SimpleBending` : la structure est identique (allocation, positions, `m1_bar`, mise au repos), seule la formule géométrique change.

## I.10 Sorties, sauvegarde et visualisation

Deux formats de sauvegarde sont fournis par `read_save_config.c` :

- **Configuration complète** (`SaveConfigLib` / `ReadConfigLib`, fichier `config.dat` ou `config_<line>_<n>.tmp`) : tous les champs de `fiber` en double précision — permet de reprendre une simulation à l'identique.
- **Configuration réduite** (`SaveReducedConfigLib` / `ReadReducedConfigLib`, fichier `reduced_Config.tmp` ou `reduced_Config_<line>_<n>.tmp`) : uniquement `xt`, `e`, `thetat`, `m1` en simple précision — c'est ce format que `core.c` écrit périodiquement pendant la simulation, et que lit le visualiseur.

Le dossier **`Visu_v9.6`** est un petit programme OpenGL/GLUT séparé qui relit les fichiers `reduced_Config_<line>_<n>.tmp` produits par une simulation et les affiche comme un "film" : au démarrage il lit un fichier `movie.txt` contenant, dans l'ordre, le chemin des données, le numéro de `line`, l'index de départ, le pas entre deux images et une temporisation, puis avance automatiquement (`glutIdleFunc(DoOneIteration)`). Il peut aussi exporter chaque image en scène **POV-Ray** (`povray.c`) pour un rendu photo-réaliste hors-ligne, avec un jeu de couleurs aléatoires par fibre (`SetColors`).

Pour un diagnostic rapide sans visualiseur, `PrintReducedConfigLib` écrit un fichier texte lisible `check_Reduced_Config.txt` (une ligne par nœud : position, $e$, $\theta$, $m_1$), et les fonctions de `energy.c` (`BendingEnergyLib`, `TwistEnergyLib`, `TractionEnergyLib`, `KineticEnergyLib`) permettent de suivre l'évolution des différentes énergies au cours du temps depuis `core.c`.

## I.11 Conseils pratiques de choix de paramètres

Le modèle est habituellement utilisé en unités **adimensionnées** (voir Partie III.5 pour la justification physique complète). Quelques repères pratiques, cohérents avec les valeurs utilisées dans les exemples et dans l'article :

- Pas de temps $\Delta t^* \approx 0.1$ dès lors que $k_n^* \sim 1$ (sinon $\Delta t^* \sim 0.1\, \min(1, (k_n^*)^{-1/2})$).
- Raideurs de contact $k_n^* = k_t^* \sim 1$.
- Amortissement interne $\lambda_{internal}^* \sim 1$–$3$ pour supprimer rapidement les ondes de compression parasites le long des segments (les exemples utilisent des valeurs entre 1 et 2.8).
- Amortissement de contact $\lambda_{n}^* \sim 1$.
- Amortissement visqueux global $\lambda^*$ doit rester petit (souvent $10^{-3}$–$10^{-4}$) : il sert uniquement à stabiliser le mouvement transverse hors contact, pas à modéliser une vraie dissipation physique.
- Toujours vérifier que la force adimensionnée typique reste petite ($f^* \sim 10^{-5}$–$10^{-3}$) pour rester dans l'hypothèse de faible élongation utilisée par le calcul discret de flexion/torsion — et vérifier la convergence en doublant `N` (le nombre de segments) et en divisant `dt` par 2.
- `nContactMax` doit être choisi assez grand pour ne jamais saturer (le programme affiche `maxPossibleContacts` à chaque recalcul de la liste de contacts potentiels, voir II.5) et doit rester un multiple de `WG`.

---

# Partie II — Algorithme et kernels

Cette partie détaille l'implémentation parallèle : organisation mémoire sur le device, séquence exacte des kernels exécutés à chaque itération, et fonctionnement interne du pipeline de contact — le point le plus spécifique de ce code, adapté à des objets non sphériques (cylindres allongés) sur GPU.

## II.1 Principe général

Toutes les grandeurs par fibre, par segment (nœud) et par contact sont stockées en **structure de tableaux** (SoA) sur le device — un tableau par champ (`fiber_xt_devPtr`, `fiber_e_devPtr`, `fiber_moment_devPtr`, ...) plutôt qu'un tableau de structures — pour un accès mémoire coalescent dans les kernels.

Chaque **segment global** est indexé par un entier `iSegment` (0 ≤ `iSegment` < `nSegment`), obtenu en concaténant les segments de toutes les fibres. Deux tableaux utilitaires permettent de revenir à l'identité "fibre + rang local" :

```c
iFiber = iFiberFromSegment[iSegment];
i      = iFromSegment[iSegment];       // 0 <= i < n[iFiber]
```
(côté host, la fonction équivalente est `NSegment(fiber*, parameter*, iFiber, i)`, qui fait l'opération inverse).

Deux constantes de taille de work-group structurent tous les kernels :

- `WG = 1024` : taille de work-group pour les kernels 1D (un thread par segment, ou par contact).
- `WG2D = 32` : taille de work-group (32×32) pour le kernel 2D de détection de contacts.

`NFPS = 128` est le nombre maximal de contributions de force/moment qu'un même segment peut recevoir en un pas de temps (accumulation via une liste, voir II.5), et `NPARAM = 100` la taille des tableaux de paramètres scalaires partagés `param_UInt`/`param_Double`.

## II.2 Table des paramètres device (`param_UInt` / `param_Double`)

Plutôt que de passer individuellement chaque scalaire à chaque kernel, deux tableaux globaux `param_UInt[NPARAM]` et `param_Double[NPARAM]` sont copiés sur le device et relus par tous les kernels via un emplacement ("slot") fixe. `Copy_Param_host2dev` (dans `device_host_device.c`) fait la correspondance :

**`param_UInt`** (extraits) :

| Slot | Contenu |
|---|---|
| 0 | `iter` (compteur d'itération) |
| 2 | `nFiber` |
| 3 | `nSegment` |
| 4 | `nContactMax` |
| 6 | `periodic` |
| 7 | nombre de contacts effectifs (écrit par `kernel_Count_Contact_Step2`) |
| 20–29 | compteurs/itérations libres (`iter0`, `iter1`, `low_High_Eta`, `iter3`, `rRate`, `low_High_Mu`, `iter6..8`, `iterStop`) |
| 30–39 | drapeaux libres (`flag0..flag7`, `flag_Phase`), et **39 = `flag_Reset_Contact`** (force la reconstruction de la liste de contacts) |

**`param_Double`** (extraits) :

| Slot | Contenu |
|---|---|
| 0 | `dt` |
| 1, 2 | `kn`, `kt` |
| 3, 4 | `lambda`, `lambda_internal` |
| 5, 6 | `lambda_contact_n`, `lambda_contact_t` |
| 7 | `mu` |
| 30 | `epsStar` (seuil courant de détection de contact, auto-ajusté) |
| 31 | déplacement max de l'itération (sortie de `kernel_Max_Displacement`) |
| 32 | déplacement intégré depuis le dernier recalcul de la liste de contacts (`Δ` de l'algorithme II.5) |
| 33–37 | accumulateurs de dissipation (visqueuse globale, étirement, contact normal/tangentiel, "W" opérateur) |
| 50–52 | `lx, ly, lz` (utilisés par `ComputeShift` pour les images périodiques) |

Les slots libres (40 et au-delà pour les doubles, 30 et au-delà pour les entiers hors ceux listés) sont à la disposition de chaque cas d'étude pour ses propres besoins (force cible, angle cible, etc.), lus dans `kernel_Specific_*` via `param_Double[...]`/`param_UInt[...]`.

## II.3 La boucle d'itération

`DoOneIteration()` (identique dans son squelette pour tous les exemples) enchaîne, à chaque pas de temps :

```c
Execute_Kernel_Calculate_le(p);          // e_i, l_i ; force visqueuse + étirement
Execute_Kernel_Bending_Force(p);         // force + moment de flexion
Execute_Kernel_Twist_Force(p);           // force + moment de torsion

Execute_Kernel_Specific_Force(p);        // [utilisateur] forces externes

Execute_Kernel_PossibleContact(p);       // (re)construction éventuelle de la liste de contacts potentiels
Execute_Kernel_Contact(p);               // géométrie exacte des contacts potentiels
Execute_Kernel_RemoveDoubleContact(p);   // dédoublonnage
Execute_Kernel_CalculateContactForce(p); // forces normale + tangentielle (Cundall-Strack)
Execute_Kernel_AddContactForce(p);       // réduction des contributions par segment
Execute_Kernel_ShiftContact(p);          // vieillissement de l'état des contacts (t-dt <- t)

Execute_kernel_Integrate_and_Shift(p);   // Verlet : x(t+dt), theta(t+dt)

Execute_kernel_Compute_m1bar(p);         // transport parallèle du repère de référence
Execute_kernel_Specific_Position(p);     // [utilisateur] contraintes géométriques
Execute_kernel_Compute_m1(p);            // applique theta pour obtenir m1 réel

Execute_Kernel_Max_Displacement(p);      // reduction : déplacement max -> decide du refresh contact
```
Ce squelette correspond très directement à l'algorithme général de l'article de référence (forces internes → forces externes → contacts → intégration → transport → contrainte géométrique). Chaque kernel n'exige de l'utilisateur aucune intervention **sauf** les deux marqués `[utilisateur]`, qui appellent directement `kernel_Specific_Force`/`kernel_Specific_Position` définis dans `kernel_specific.cl` du cas d'étude (voir I.7).

## II.4 Détail des kernels de forces internes

**`kernel_Calculate_le`** (un thread par segment) — calcule pour chaque segment $e^i$ et $l^i$ (à $t$ et $t-\Delta t$), remet à zéro la force/moment du segment en y plaçant : la force visqueuse globale $-\lambda(x_t - x_{tm})/\Delta t$, la force et l'amortissement d'étirement longitudinal (ressort $k_0$ + amortisseur $\lambda_{internal}$ entre nœuds voisins), et accumule les dissipations correspondantes dans `dissip_visc_global`/`dissip_visc_stretch` (utilisées pour un bilan d'énergie, cf. `Execute_Kernel_Dissipation`).

**`kernel_BendingForce`** (un thread par segment, appelé avec `local_size = WG/4`) — implémente la discrétisation de la courbure et de la force de flexion donnée en Partie III.2 : pour chaque segment $i$, on relit une fenêtre de 5 nœuds ($x_{i-2}\dots x_{i+2}$) et de $e$, $m_1$ voisins, on calcule la dérivée seconde discrète $d^2r/ds^2$ à trois positions consécutives ($i-1$, $i$, $i+1$), on en déduit la courbure réduite $\kappa_{(1,2)} - \bar\kappa_{(1,2)}$ à chacune, et on cumule la contribution de force sur le segment courant ainsi que le moment de flexion résultant (projeté sur l'axe $e^i$, cf. la décomposition force/moment de la Partie III.2).

**`kernel_TwistForce`** (un thread par segment, `local_size = WG/2`) — implémente l'énergie de torsion discrète $\propto \beta_i^2$ (Partie III.2, formule du "produit mixte" $\beta_i = (e^{i-1}+e^i)\cdot(m^i\times m^{i-1})$) : calcule $\beta$ et le moment symétrisé sur les segments $i-1$, $i$, $i+1$, en déduit la force perpendiculaire équivalente (couple de forces) et la met à jour.

Ces trois kernels lisent/écrivent en simple précision (`float`, `float3`) pour les grandeurs orientées (e, m1, forces, moments) mais gardent les positions/angles en double précision (`double`, `double3`) — compromis précision/performance classique en DEM sur GPU, la position cumulée devant rester précise sur de longues simulations alors que les forces sont recalculées à chaque pas.

## II.5 Détail du pipeline de contact

Le traitement des contacts est la partie la plus spécifique de ce code : les objets ne sont pas des sphères mais des **cylindres allongés reliés entre eux**, ce qui interdit d'utiliser telles quelles les méthodes classiques de recherche de voisins des codes DEM granulaires (voir la discussion de l'article, section III.C.3, sur le choix entre cellules de taille $\sim l_0$ ou $\sim r$).

**a) Maintien d'une liste de contacts potentiels (broad-phase), `kernel_Poss_Contact`.**
Recalculer à chaque pas de temps la distance entre **toutes** les paires de segments serait $O(N_{seg}^2)$ et inutilement coûteux : on ne le refait que lorsque c'est nécessaire. Un compteur $\Delta$ (`param_Double[32]`) accumule deux fois le déplacement maximal mesuré à chaque itération (`kernel_Max_Displacement`) ; tant que $\Delta < \epsilon^\star$ (`epsStar`, `param_Double[30]`), on est certain qu'aucun segment hors de la liste actuelle n'a pu entrer en contact, et on **saute entièrement** la reconstruction (`Execute_Kernel_PossibleContact` retourne immédiatement). Quand $\Delta \ge \epsilon^\star$, on relance la recherche complète :

- `kernel_Refresh_and_Store_Poss_Contact` sauvegarde l'ancienne liste (`oldPossContSeg1/2`) et vide la nouvelle ;
- `kernel_Poss_Contact` est lancé sur une grille **2D** `nSegment × nSegment` (work-groups `32×32`), avec `iSeg2 > iSeg1` uniquement (moitié de matrice), et copie d'abord les positions en mémoire locale pour limiter les accès à la mémoire globale ; pour chaque paire, il teste la distance minimale cylindre-cylindre / cylindre-sphère / sphère-sphère (formules de la Partie III.4) élargie d'une marge $\epsilon^\star/l$, et si elle est inférieure à $r_1+r_2+\epsilon^\star$, ajoute la paire à la liste ;
- pour équilibrer la charge entre threads malgré la forte hétérogénéité du nombre de contacts par segment le long d'une même fibre, chaque contact potentiel n'est pas simplement ajouté en fin de liste globale mais réparti dans un **groupe** `iGroup` déterminé par une fonction de `iSeg2` (repliée autour du milieu de la fibre, cf. formule ci-dessous) de sorte que les `NcontactMax / WG` groupes reçoivent une charge à peu près égale ;
- `kernel_Align_Poss_Contact` réordonne ensuite chaque groupe pour que les contacts déjà existants au pas précédent gardent le **même indice** `iCont` (nécessaire pour préserver l'historique du déplacement tangentiel $u_t$ d'un pas à l'autre, voir III.3) ;
- enfin `epsStar` est **auto-ajusté** : si un groupe dépasse `3·WG/4` contacts potentiels on le réduit (`epsStar -= 0.02`), sinon on l'augmente légèrement (`epsStar += 0.01`, plancher `0.01`) — ce qui maintient un bon compromis entre fréquence de reconstruction et remplissage des listes.

La fonction de répartition en groupe (implémentée dans `kernel_Poss_Contact`, et rappelée dans la documentation d'origine du code) est :

```
iAux  = (2*(iSeg2+1) > nSegment) ? nSegment - iSeg2 - 1 : iSeg2
iGroup = (2*iAux * (nContactMax / (WG*32))) / (nSegment/32)
```

**b) Détermination géométrique + forces (narrow-phase), un thread par contact potentiel.**
- `kernel_Contact` recalcule, pour chaque entrée de la liste, la géométrie exacte du contact (distance minimale entre les deux segments considérés comme couple cylindre+sphères terminales, cf. III.4), en choisissant le cas (sphère-sphère / cylindre-sphère / sphère-cylindre / cylindre-cylindre) donnant l'interpénétration $\delta$ la plus grande ; il en déduit la normale $n$, le point de contact $C$, et marque l'état du contact `NEWCONTACT` ou `STILLCONTACT` (selon qu'il existait déjà au pas précédent, état `OLDCONTACT`).
- `kernel_RemoveDoubleContact` élimine les doublons géométriques qui peuvent apparaître à la jonction entre deux segments consécutifs d'une même fibre (même point de contact détecté deux fois), en ne gardant que celui de plus grande interpénétration.
- `kernel_CalculateContactForce` calcule la force normale (ressort-amortisseur), met à jour et seuille le déplacement tangentiel `u_t` (modèle de Cundall–Strack + Coulomb, voir III.3), calcule les composantes de force à répartir sur les deux nœuds de chaque segment en contact (pondération par l'abscisse curviligne $s_1, s_2$ du point de contact), et les **accumule par segment** de façon atomique dans deux tableaux `force_per_Segment[NSegment × NFPS]` / `moment_per_Segment[...]` (chaque segment ne pouvant recevoir plus de `NFPS=128` contributions par pas de temps).
- `kernel_AddContactForce` (un thread par segment) fait la somme finale des contributions de contact et les ajoute à `fiber_f`/`fiber_moment`.
- `kernel_ShiftContact` fait vieillir l'état : `STILLCONTACT`/`NEWCONTACT → OLDCONTACT` (et `deltatm ← deltat`, `utm ← ut`, pour le pas suivant), sinon `→ NOCONTACT`.

**c) Machine à états d'un contact.** Un emplacement de la liste passe par les états `NOCONTACT → NEWCONTACT → STILLCONTACT (→ STILLCONTACT ...) → NOCONTACT`, via l'intermédiaire `OLDCONTACT` entre deux pas de temps ; c'est ce qui permet de savoir, à chaque pas, si le déplacement tangentiel accumulé `u_t` doit être remis à zéro (nouveau contact) ou prolongé (contact persistant), condition nécessaire pour un frottement de Coulomb physiquement correct.

## II.6 Intégration temporelle et transport du repère matériel

**`kernel_Integrate_and_Shift`** applique un schéma de **Verlet sans vitesse explicite** (aussi appelé Störmer–Verlet) séparément sur les positions et les angles :
```
x(t+dt) = 2x(t) - x(t-dt) + f(t)·dt²/masse
theta(t+dt) = 2·theta(t) - theta(t-dt) + moment(t)·dt²/j
```
(et décale `xtm ← xt`, `thetatm ← thetat`), en ignorant les fibres qui ne sont pas `STATUS_FREE`.

**`kernel_Compute_m1bar`** transporte ensuite le repère de référence $\bar m_{(1)}$ **sans torsion ajoutée** le long du nouveau segment : transport parallèle de $e_{tm}$ vers $e_t$ (fonction `parallel_transport` de `functions_v9.6.cl`, $x \mapsto x + (e_{tm}\times e_t)\times x$), projection pour rester perpendiculaire à $e_t$, puis renormalisation. C'est cette étape qui garantit que $\theta_i$ mesure bien la torsion "physique" ajoutée par rapport au transport géométrique naturel, indépendamment de la façon dont la fibre se courbe.

**`kernel_Compute_m1`** applique enfin la rotation d'angle $\theta_i$ (stocké dans `thetat[...][2]`) autour de $e^i$ pour obtenir le repère matériel réel $m_1$ (fonction `rotate` de `functions_v9.6.cl`, formule de Rodrigues plane) — sauf pour les fibres réduites à une seule sphère (`n==1`), traitées comme un corps rigide dont on tourne directement `m1` et `e` de $\Delta\theta$.

`kernel_Unwarp` (utilitaire, appelé séparément si besoin) détecte et corrige les sauts de $\pi$ dans $\theta_i$ le long d'une fibre, utile pour "dérouler" un historique de torsion continu à des fins d'analyse.

## II.7 Conditions périodiques

Lorsque `periodic = PERIODIC_XY` (ou `PERIODIC_XYZ`/`PERIODIC_Z`), la fonction `ComputeShift(periodic, iShift, dx, dy, dz, &shift)` (`functions_v9.6.cl`) renvoie l'un des 9 vecteurs de translation d'image périodique dans le plan $(x,y)$ (centre + 8 voisins), à partir de la taille de boîte `lx, ly` (`param_Double[50..51]`). `kernel_Poss_Contact` est alors lancé avec une troisième dimension de grille `global_Size3D[2] = 9` (une par image), chaque contact potentiel retenu mémorisant l'image `iShift` utilisée (`iShiftCont`), reproduite ensuite par `kernel_Contact`/`kernel_CalculateContactForce` (lignes actuellement commentées dans le code fourni — l'ossature est en place mais le calcul du shift y est désactivé par défaut ; à réactiver si votre cas d'étude a besoin de conditions périodiques effectives).

## II.8 Compilation des kernels et transferts host ↔ device

`Create_Kernels()` (dans `device_init.c`) concatène en mémoire, dans cet ordre précis, les sources suivantes avant de les compiler en un seul programme OpenCL (`clBuildProgram`) :

1. `fiberLib_Common_Macros_v9.6.h` (macros/constantes),
2. `functions_v9.6.cl` (fonctions utilitaires),
3. `kernel_specific.cl` **du cas d'étude courant**,
4. `kernel_max_v9.6.cl`, `kernel_le_force_v9.6.cl`, `kernel_poss_contact_v9.6.cl`, `kernel_contact_v9.6.cl`, `kernel_integrate_shift_v9.6.cl`.

Chaque kernel a ensuite tous ses arguments fixés **une fois pour toutes** par `clSetKernelArg` (pas à chaque itération), ce qui impose que la signature de `kernel_Specific_Force`/`kernel_Specific_Position`/`kernel_Specific_OneTime` dans `kernel_specific.cl` corresponde exactement à ce que `Create_Kernels()` attend (voir la liste d'arguments dans `device_init.c` si vous ajoutez un buffer à un kernel spécifique).

Côté transferts, `device_host_device.c` fournit les fonctions symétriques `Create_Host_Ptrs`/`Create_Device_Ptrs` (allocation), `Copy_Fiber_host2dev`/`Copy_Fiber_dev2host` (l'intégralité de l'état d'une fibre), `Copy_FiberPosition_dev2host` (positions seules, plus léger), `Copy_Param_host2dev`/`Copy_Param_dev2host`, et `Copy_Contact_dev2host` (relit la liste de contacts effective dans `gContact[]`, utile pour du post-traitement sur le host). Ces transferts complets ne sont typiquement effectués que périodiquement (voir `core.c`, toutes les 10⁴ itérations dans les exemples), le calcul lui-même restant intégralement sur le device d'une itération à l'autre.

---

# Partie III — Annexe : modèle physique

Cette annexe résume le modèle mécanique implémenté par les kernels, tel que publié dans l'article de référence (Crassous, *Phys. Rev. E* **107**, 025003, 2023) ; s'y reporter pour les démonstrations complètes, les tests de validation (poutre encastrée, flambage en torsion, capstan, nœuds élastiques, chaîne tombante, modèle de fil retors) et les comparaisons expérimentales.

## III.1 Description mécanique de la fibre

Une fibre est modélisée (Bergou et al., *discrete elastic rods*) comme un ensemble de $N$ points connectés $r_i$, $0\le i \le N-1$, reliés par $N-1$ segments de vecteur directeur $e^i = (r_{i+1}-r_i)/\lVert r_{i+1}-r_i\rVert$ et de longueur $l^i$. Le segment $i$ est géométriquement un cylindre de diamètre $d$ ; chaque nœud $r_i$ porte en outre une sphère de même diamètre — l'enveloppe complète de la fibre pour la détection de contact est donc l'union de ces cylindres et sphères (nécessaire notamment pour bien traiter la jonction entre deux segments consécutifs). Une masse $m_0$ est associée à chaque nœud, un moment d'inertie $J$ à chaque cylindre.

Les nœuds peuvent se translater (flexion, étirement) ; les segments restent rectilignes mais peuvent tourner autour de leur propre axe (torsion). L'état cinématique complet est donc $\{r_i\}_{0\le i\le N-1}$ et $\{\theta_i\}_{0\le i\le N-2}$. Toute force ou moment agissant sur un segment (élastique ou de contact) est décomposé en un couple axial (le long de $e^i$) et des forces appliquées aux deux nœuds extrémités — cette décomposition est détaillée pour la torsion en III.2 et pour le contact en III.3.

## III.2 Énergies élastiques discrètes

**Étirement.** Pour une raideur $k_0$ et un amortissement $\lambda$, la force d'étirement exercée par le point $i+1$ sur le nœud $i$ est
$$f^{(e)}_{i+1;i} = \big[k_0(l^i - l_0) + \lambda\,\dot l^i\big]\,e^i ,$$
chaque nœud intérieur recevant la somme des contributions de ses deux segments voisins.

**Flexion.** À partir de l'énergie continue $E^{(b)} = \int_s (B/2)\kappa^2\,ds$, l'énergie discrète s'écrit, pour $B$ constant le long de la fibre,
$$E^{(b)} = \frac{Bl_0}{2}\sum_{i=1}^{N-2}\kappa_i^2 ,\qquad
\kappa_i^2 = \frac{4\,(l^{i-1})^2(l^i)^2 - (l^{i-1}\cdot l^i)^2}{(l^{i-1})^2 (l^i)^2 (l^{i-1}+l^i)^2}$$
où $\kappa_i$ est la courbure au nœud $i$ (cercle passant par $r_{i-1}, r_i, r_{i+1}$, formule de Héron). Dans la limite d'une fibre faiblement étirée ($l^i \simeq l_0$) et faiblement courbée ($l^i \simeq l^{i-1}$), utilisée par le kernel `kernel_BendingForce`, la force de flexion se réduit à un opérateur linéaire à 5 points :
$$f^{(b)}_i = -\frac{B}{l_0^3}\big[r_{i-2} - 4r_{i-1} + 6r_i - 4r_{i+1} + r_{i+2}\big]$$
(avec des expressions particulières aux deux premiers/derniers nœuds). Dans l'implémentation retenue par ce code, la courbure est en outre **projetée** sur les deux directions matérielles $m_{(1)}, m_{(2)}$ (plutôt que sur la seule normale de Frenet, non définie pour une fibre rectiligne), avec une **courbure naturelle** $\bar\kappa_{(j),i}$ soustraite :
$$\kappa_{(j),i} = \Big(\frac{dt}{ds}\Big)_i \cdot m_{(j),i} = \frac{1}{2l_0^2}(r_{i-1}-2r_i+r_{i+1})\cdot(m^{i-1}_{(j)}+m^i_{(j)}),$$
$$E^{(b)} = \sum_{j=1,2}\frac{B_j l_0}{2}\sum_{i=1}^{N-2}\big(\kappa_{(j),i}-\bar\kappa_{(j),i}\big)^2 .$$
Chaque terme d'énergie $E^{(b)}_i$ produit des forces sur les trois nœuds $i-1,i,i+1$ (Éq. 18–21 de l'article) ; comme $\partial m_{(j)}/\partial r_i$ n'est pas rigoureusement nul, ces trois forces seules ne conservent pas le moment cinétique — un moment correctif $M$ est calculé (Éq. 22–24) puis redécomposé en une paire de forces perpendiculaires à $e^{i-1}$/$e^i$ (Éq. 30–32), c'est ce calcul complet que réalise `kernel_BendingForce`.

**Torsion.** À partir de l'énergie continue $E^{(t)} = \int_s (C/2)\,\tau^2\,ds$ ($\tau$ = torsion totale, somme d'une torsion "interne" $(\theta_i-\theta_{i-1})/l_0$ et d'une torsion géométrique de la ligne centrale $\tau_s$), l'énergie discrète, en utilisant $m_{(1)}\cdot m_{(2)}=0$, se réécrit en fonction du produit mixte
$$\beta_i = (e^{i-1}+e^i)\cdot\big(m^i \times m^{i-1}\big), \qquad
E^{(t)} = \frac{C}{8l_0}\sum_{i=1}^{N-2}\beta_i^2 .$$
Le moment de torsion élastique sur le segment $(i,i{+}1)$ est
$$m^{(t)}_i = \frac{C}{l_0}\Big[(\theta_{i+1}-\theta_i+l_0\tau_{s,i+1})\,e^{i+1} - (\theta_i - \theta_{i-1} - l_0\tau_{s,i})\,e^{i-1}\Big],$$
dont la composante axiale $m_i^{(t)} = m^{(t)}_i\cdot e^i$ est le couple retenu, et dont la composante perpendiculaire est réexprimée comme un couple de forces $f^{(t)}_i = -f^{(t)}_{i+1} = (m_i^{(t)}/l^i)\times e^i$ appliqué aux deux nœuds extrémités du segment. C'est ce calcul (moments $M^i_i$, $M^{i-1}_i$ sur les segments voisins, symétrisation pour assurer $M^i_{s,i}+M^{i-1}_{s,i}=0$, puis décomposition longitudinale/transverse) que met en œuvre `kernel_TwistForce`.

## III.3 Modèle de contact frictionnel (Cundall–Strack)

Le contact entre deux segments (ou entre segment et sphère terminale, ou entre deux sphères) est modélisé comme un ressort-amortisseur normal plus un frottement de Coulomb tangentiel régularisé, à la manière des codes DEM granulaires classiques (modèle de Cundall & Strack, 1979) :
$$f^{(c)}_n = -\big[k_n\,\delta + \lambda_n\,\dot\delta\big]\,n, \qquad
f^{(c)}_t = -\min\!\big(k_t\lVert u_t\rVert,\ \mu k_n\delta\big)\,\frac{u_t}{\lVert u_t\rVert}$$
où $\delta = (d_1+d_2)/2 - d(t)$ est l'interpénétration (distance minimale $d(t)$ entre axes des cylindres, cf. III.4), $n$ la normale, $\mu$ le coefficient de frottement, et $u_t$ le déplacement tangentiel relatif accumulé **depuis la formation du contact**. Ce déplacement est mis à jour à chaque pas de temps par :

1. rotation de $u_t(t-dt)$ pour suivre la rotation de la normale entre $t-dt$ et $t$ (transport de la mémoire tangentielle lorsqu'un point de contact glisse sur la surface d'un cylindre) ;
2. ajout de la vitesse relative des points matériels coïncidant avec le point de contact, $u_t(t) = u_t^{rot}(t-dt) + (v^{(2)}-v^{(1)})\,dt$, où chaque vitesse combine la translation du nœud et la rotation (axiale $\dot\theta_i e^i$ + transverse $\Omega_\perp = \frac1{l^i}e^i\times(\dot r_{i+1}-\dot r_i)$) du segment porteur ;
3. projection pour rester perpendiculaire à $n$, puis écrêtage à $u_t \leftarrow \mu k_n\delta/k_t$ si $k_t\lVert u_t\rVert > \mu k_n\delta$ (loi de Coulomb).

La résultante de contact $f^{(c)} = f^{(c)}_i + f^{(c)}_{i+1}$ et le moment $(r_{i+1}-r_i)\times f_{i+1} + m_i e^i = (r_C-r_i)\times f^{(c)}$ ne déterminent pas de façon unique la répartition entre les deux nœuds extrémités du segment ; le choix retenu (justifié en détail dans l'article, Appendix A3) est
$$m^{(c)}_i = \big[(r_C-r_i)\times f^{(c)}\big]\cdot e^i,\qquad
f^{(c)}_i = (1-s)\,f^{(c)} + \frac{R}{l}(f^{(c)}\cdot e^i)\,n,\qquad
f^{(c)}_{i+1} = s\,f^{(c)} - \frac{R}{l}(f^{(c)}\cdot e^i)\,n,$$
où $s\in[0,1]$ est l'abscisse curviligne du point de contact sur le segment. C'est directement cette formule (pondération par `s1`/`s2`, composante normale $R/l\,(f\cdot e)$) qu'implémente `kernel_CalculateContactForce`.

Une force de **viscosité globale** $f^{(v)}_i=-\lambda_v \dot r_i$ peut être ajoutée pour amortir les mouvements transverses hors contact (le modèle élastique de fibre seul n'a sinon aucune dissipation perpendiculaire à son axe), ainsi que des forces volumiques (gravité) ou des forces externes propres à chaque géométrie (traction aux extrémités, etc.) — c'est précisément le rôle du kernel spécifique `kernel_Specific_Force` (Partie I.7 / II.3).

## III.4 Distance entre segments / recherche de contacts

Pour deux segments d'axes portant chacun une sphère terminale (rayon $r$) et un tronçon cylindrique ($0\le s\le 1$), on paramètre chaque point de l'axe par son abscisse $s$ ; la distance au carré entre deux points d'abscisses $s_1, s_2$ sur les deux axes s'écrit $d^2(s_1,s_2) = \lVert a + s_1 b - s_2 c\rVert^2$ où $a$ relie les deux origines et $b,c$ sont les vecteurs des deux segments. La distance minimale sur $s_1,s_2\in\mathbb R$ s'obtient en annulant les deux dérivées partielles (système linéaire $2\times2$) :

$$s_1^\star = \frac{cc\cdot ab - ac\cdot bc}{bb\cdot cc - bc^2}, \qquad
s_2^\star = -\frac{bb\cdot ac - ab\cdot bc}{bb\cdot cc-bc^2}$$
(notations $aa=a\cdot a$, $ab=a\cdot b$, etc.) Si $(s_1^\star,s_2^\star)\in[0,1]^2$, le contact **cylindre–cylindre** est retenu ; sinon on cherche, dans l'ordre, un contact **cylindre–sphère** ($s_1=ab/bb$, avec $s_1\in[0,1]$), **sphère–cylindre** ($s_2=-ac/cc$), puis en dernier recours **sphère–sphère** ($s_1=s_2=0$). C'est exactement cette cascade de quatre cas (`deltaCC`, `deltaCS`, `deltaSC`, `deltaSS`, en retenant l'interpénétration $\delta = r_1+r_2-d$ la plus grande) qu'implémentent `kernel_Poss_Contact` (avec une marge $\epsilon^\star$ élargissant le domaine $[0,1]$ pour anticiper le mouvement entre deux reconstructions de liste, cf. II.5) et `kernel_Contact` (calcul exact, sans marge, une fois la paire retenue).

Pour un système de fibres fortement anisotropes ($l_0 \gg r$), les méthodes classiques de partitionnement en cellules (*linked cells*) des codes DEM granulaires sont inefficaces : une cellule dimensionnée sur $l_0$ contient trop de segments, une cellule dimensionnée sur $r$ en dénombre trop. Le choix retenu ici — exploiter le fait que les segments d'une même fibre sont connectés pour parcourir les candidats par fibre voisine croissante (cf. Fig. 4c de l'article) — est ce qui motive l'organisation par groupes équilibrés de la Partie II.5, plutôt qu'un maillage spatial classique.

## III.5 Échelles physiques et choix des paramètres numériques

Le code travaille en unités adimensionnées : échelle de masse $m_0$ (masse d'un nœud), échelle de longueur $l_0$ (longueur de segment au repos), échelle de raideur $k_0$ (raideur d'étirement). Pour une fibre élastique réelle de rayon $r$, module d'Young $E$, coefficient de Poisson $\nu$, masse volumique $\rho$ :
$$k_0 = \frac{E\pi r^2}{l_0}, \qquad m_0 = \rho\pi r^2 l_0, \qquad t_0 = \sqrt{m_0/k_0} = l_0\sqrt{\rho/E}, \qquad f_0 = k_0 l_0 = E\pi r^2 .$$
$t_0$ est le temps de propagation d'une onde de compression sur un segment ; $f_0$ la force qui étirerait de 100% une fibre parfaitement élastique. Pour toute grandeur $x$ d'échelle $x_0$, on note $x^\* = x/x_0$.

En pratique, on choisit la force typique adimensionnée $f^\*\sim 10^{-5}$–$10^{-3}$ pour rester dans l'hypothèse de faible élongation (si $f^\*$ est trop petit en revanche, la propagation des ondes transverses, de vitesse $v_t^\* = (f^\*)^{1/2}$ en l'absence de flexion, devient très lente). Les modules adimensionnés de flexion et de torsion valent, pour une fibre cylindrique homogène,
$$B^\* = \frac{B}{k_0 l_0^3} = \frac{(r^\*)^2}{4}, \qquad C^\* = \frac{C}{k_0 l_0^3} = \frac{(r^\*)^2}{2(1+\nu)} .$$
L'amortissement longitudinal est pris $\lambda \sim \sqrt{k_0 m_0}$ (soit $\lambda^\*\sim 1$) pour supprimer rapidement les ondes de compression parasites.

Pour le contact, la raideur normale est estimée en linéarisant le contact hertzien entre deux cylindres (équivalent à une sphère de rayon $r$ contre un plan, $f_n = \frac{4}{3}E_{eff}\sqrt r\,\delta^{3/2}$) autour de la force de traction typique appliquée, ce qui conduit à
$$k_n^\* \approx \frac{(f^\*)^{1/3}}{r^\*}$$
(à un facteur numérique près) ; en pratique $f^\*\sim10^{-5}$–$10^{-3}$ et $r^\*\sim10^{-1}$ donnent $k_n^\*\sim 1$, d'où le choix usuel $k_t^\*=k_n^\*=1$ et un amortissement de contact $\lambda_n^\*\sim1$ pour une relaxation rapide de l'oscillation de contact. Le pas de temps est enfin choisi pour résoudre correctement à la fois la relaxation de longueur de segment (temps $t_0$) et l'établissement du contact (temps $t_0\sqrt{k_0/k_n}$) :
$$dt^\* = \frac1{10}\min\!\big(1,\ (k_n^\*)^{-1/2}\big) \ \Longrightarrow\ dt^\*=0.1 \text{ lorsque } k_n^\*\sim1,$$
valeur effectivement utilisée par défaut dans tous les exemples de ce dépôt. On vérifiera toujours, pour un jeu de paramètres donné, que les résultats sont inchangés en divisant `dt` par 2 (convergence temporelle) et en doublant le nombre de segments (convergence spatiale, la convergence de la flexion discrète étant en $N^{-2}$, cf. Fig. 5c de l'article).

## III.6 Référence

J. Crassous, *Discrete-element-method model for frictional fibers*, **Phys. Rev. E 107, 025003 (2023)**, DOI: [10.1103/PhysRevE.107.025003](https://doi.org/10.1103/PhysRevE.107.025003). L'article contient en outre quatre illustrations complètes (tige encastrée / flambage en torsion, capstan sans flexion, nœuds élastiques ouverts avec et sans frottement, impact d'une chaîne tombante, modèle de fil retors à fibres multiples) qui constituent autant de cas de validation supplémentaires du modèle décrit ci-dessus.
