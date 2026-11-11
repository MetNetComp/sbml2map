# sbml2map

sbml2map draws metabolic maps from SBML models. It can export Escher maps,
SVG images, and SBML files with compound and reaction positions and connecting lines.
To export SBML with positions, use `--export-sbml-layout`.

The structural Escher view paints red above orange above gray and offsets coincident routes of different reactions.
Biological nodes stay fixed; same-reaction trunks remain bundled. Primary backbone
centerlines stay fixed during parallel-lane separation, but may detour around unrelated red/orange nodes.
Parallel red/orange segments from different reactions that run together for at
least one grid unit are separated by approximately 0.3 grid units using rectangular
detours in both SVG and Escher structural views. Node collisions are avoided;
unresolved routes are recorded in `parallel_routing.json` (SVG) and map metadata
(Escher). Display detours add length without changing the original layout objective.
FBA uses a dedicated complete map with neutral orthogonal routes, retaining the saved pathway coordinates.

`sbml2map.py` draws metabolic networks with nodes on a grid.
Cofactors may also be placed off the grid. Inorganic species are ignored.
Exchange reactions are placed on the outer boundary. Key pathways, such as
glucose to pyruvate, are drawn as straight lines and highlighted in red.

When flux sampling is available, reactions with zero median flux are not drawn.
Orange paths connect the key pathways to exchange reactions. With flux sampling,
the layout favors short paths through reactions with larger absolute median flux.
Use `--orange-flux-weight` to control this preference.

The criteria for selecting key pathways are provisional. The drawing criteria
are also provisional, and we welcome all ideas for improvement.

## Two useful options

| Option | What it does |
| --- | --- |
| `--export-sbml-layout` | Saves compound and reaction positions and connecting lines in `layout.sbml`. Requires SBML Level 3 input and the optional libSBML dependency. |
| `--orange-flux-weight NUMBER` | Controls the penalty for low-flux orange routes. The default is `1`. Use `0` for length only, or a larger number to favor higher-flux reactions more strongly. |

The orange penalty uses flux sampling results. It has no effect with
`--no-sampling`. SBML Layout stores positions and lines; colors are not included.
See [positioned SBML export](#export-positioned-sbml-optional) and
[orange route penalty](#orange-route-penalty) for examples.

## Browse saved results and publish with GitHub Pages

Open [`results/e_coli_core/index.html`](results/e_coli_core/index.html) for the SVG
viewer, [`escher_preview.html`](results/e_coli_core/escher_preview.html) for the
structural Escher view, or [`escher_fba.html`](results/e_coli_core/escher_fba.html)
for the FBA result. The structural view starts without a flux overlay; FBA uses
only simulation-driven colors and widths.

The ready-to-publish [`docs/`](docs/README.md) directory contains the E. coli core
views. After updating results, run `python3 prepare_github_pages.py` and commit/push
`docs/`. An enabled Pages site can use branch `main`, directory `/docs`.

The iCH360 structural drawing is also available in the [live SVG viewer](https://metnetcomp.github.io/sbml2escher/eco_ich_360/)
and [live Escher viewer](https://metnetcomp.github.io/sbml2escher/eco_ich_360/escher_preview.html).
It shows 326 of 349 reactions without flux sampling; no new model files are included.

Run the packaged checks with `python -m unittest discover`.

## Escher labels

Use `--hide-metabolite-labels` and/or `--hide-reaction-labels` to initially hide
labels in the Escher preview. Use `--allow-label-overlap` to show all enabled
labels, including overlapping ones:

```bash
python sbml2map.py --model models/e_coli_core.xml --no-sampling \
  --hide-metabolite-labels --hide-reaction-labels --output-dir /tmp/ecore_no_labels
```

By default, overlapping labels are filtered together, with priority given to
EX reactions and their connected compounds, red pathways, orange connectors,
other reactions and primary compounds, then secondary compounds. Ties use IDs
for stable ordering. Cofactor copies remain low priority except at EX reactions.
Priorities come from structural metadata and stay the same when flux colors change.

The generated HTML provides `Metabolite labels`, `Reaction labels`, and
`Allow label overlap` checkboxes. Hidden label types stay hidden even when overlap
is allowed. Bounds are measured again after zooming, resizing, and label edits.
Names, IDs, chemistry and coordinates are preserved. These options apply to
interactive SVG and Escher views, including FBA. Static PNG/SVG downloads retain
the complete annotation set. Important SVG/PNG labels have white halos and are
painted above lines; Escher labels also use white halos without changing colors.

## Installation

Run these commands from the project root:

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements-layout.txt -r requirements-escher.txt
```

On Ubuntu, if virtual environment creation fails because `ensurepip` is
unavailable, install the venv package and then rerun the commands above:

```bash
sudo apt update
sudo apt install python3-venv
```

If the error message names a version-specific package, install that package
instead (for example, `python3.12-venv` for Python 3.12).

The required Python packages are NumPy, SciPy, Matplotlib, and jsonschema.
Optimization uses HiGHS through SciPy; no separate solver executable or MATLAB
installation is required. A browser is needed to view the interactive output.
The Escher JavaScript viewer is bundled under `vendor/escher/`.

Run all examples from the project root with the virtual environment active.
Always specify `--model` to select an input file.

## Export positioned SBML (optional)

SBML files are **not generated by default**. To explicitly export a positioned
SBML Level 3 model, install the optional dependency and pass the export flag:

```bash
python -m pip install -r requirements-sbml-layout.txt
python sbml2map.py --model models/e_coli_core.xml --no-sampling \
  --export-sbml-layout --output-dir results/runs/e_coli_core_layout
```

This writes `layout.sbml` and `sbml_layout_report.json`. Without the flag, these
files are neither created nor updated and libSBML is not required. JSON and SBML
Level 2 inputs are not supported by this exporter. `--discover-only` cannot be
combined with layout export.

To export an existing saved map without recalculating layout or samples:

```bash
python export_sbml_layout.py --model models/e_coli_core.xml \
  --layout-dir results/e_coli_core --output results/e_coli_core/layout.sbml --verify-fba
```

The original full SBML model is preserved, including FBC objectives, bounds,
stoichiometry, genes, notes, and annotations. Only displayed entities receive
glyphs: the E. coli core example retains 95 source reactions and positions 76.
Local cofactor copies become multiple SpeciesGlyphs referencing one species.
Curves preserve the saved horizontal/vertical segments and lane offsets, using
the original SBML species/reaction IDs. Existing layouts are preserved and a
unique new layout ID is chosen. `--verify-fba` additionally checks the full-model
objective before and after the export. Geometry uses the standard Layout
package; Render colors and arrow styling are not exported in this first version.

## Directory structure

```text
./
  README.md
  sbml2map.py                  Main command-line entry point
  *.py                           Supporting Python modules
  requirements-layout.txt        Numerical and plotting dependencies
  requirements-escher.txt        Escher schema-validation dependency
  requirements-layout-lock.txt   Recorded layout dependency versions
  data/bacteria/                 Chemical-role reference catalog
  major_pathway.json             Canonical backbone reference
  pathway_output/                Additional canonical backbone references
  vendor/escher/                 Bundled viewer, schema, license, and provenance
  models/                        Input SBML models
  prepare_github_pages.py        Refresh the E. coli core Pages files
  docs/                          E. coli core website and downloads
  results/                       Saved drawing results
    manifest.json                Model-to-result mapping and verification record
```

## Input models: `models/`

This folder contains seven SBML models. The program also accepts BiGG/COBRA
model JSON through `--model`. Escher map JSON is an output format and cannot
be used as an input model.

| File | Description |
| --- | --- |
| [e_coli_core.xml](models/e_coli_core.xml) | E. coli core metabolic model |
| [iML1515.xml](models/iML1515.xml) | E. coli K-12 MG1655 model from BiGG |
| [iMM904.xml](models/iMM904.xml) | Saccharomyces cerevisiae model from BiGG |
| [Human-GEM.xml](models/Human-GEM.xml) | Human1, Human-GEM v1.19.0 |
| [iNESC2DN_CMap_Annotated.xml](models/iNESC2DN_CMap_Annotated.xml) | NeuroMap CMap with original CellDesigner annotations |
| [iNESC2DN_CTMap_Annotated.xml](models/iNESC2DN_CTMap_Annotated.xml) | NeuroMap CTMap with original CellDesigner annotations |
| [iNESC2DN_ComboMap_Annotated.xml](models/iNESC2DN_ComboMap_Annotated.xml) | NeuroMap ComboMap with original CellDesigner annotations |

The SBML reader targets metabolic models with fixed stoichiometry. Dynamic
rules, events, initial assignments, and variable stoichiometry are unsupported.
Models keep their original attribution and usage conditions.

Chemical-role and reaction-identifier matching also accept `isRelatedTo`
annotations. Explicit `is` annotations and recognized model IDs take precedence over related identifiers.
Conflicting related identifiers remain ambiguous and can be resolved with
`--identity-map`; their evidence is recorded in `id_mapping.json`.

MetaNetX chemical IDs are read from annotations or species IDs such as
`MNXM23@MNXC3`. Identifier namespaces and ChEBI prefixes are normalized.
Structured notes `SOURCE: M_f6p_c` provide supplementary exact-ID evidence.
To use an additional local MNXref release, pass `--mnxref-dir /path/to/release`
containing `chem_xref.tsv` and `chem_depr.tsv`. The release headers and SHA-256
hashes are saved in `pathway_extraction.json`; no automatic download is made.
Deprecated IDs with multiple replacements stay ambiguous unless stronger
identity evidence or an explicit override resolves them.

When standard chemical identifiers are missing, the reader checks structured
SBML notes (names, abbreviations, formulas, and VMH links) and CellDesigner
species names. Matches use exact reference identifiers or names, including
explicit legacy D/L aliases. Ambiguous matches stay unresolved and can be set
with `--identity-map`. The evidence is saved in `id_mapping.json`.

Exchange evidence includes SBO:0000627, BiGG EX annotations, and notes SOURCE
labels. An explicit model-boundary compartment (including MetaNetX `BOUNDARY`)
can supply the second participant of an exchange. A matching chemical,
opposite coefficient signs, and boundary evidence on every use are required.
Verified boundary participants and CellDesigner source/sink symbols are hidden
only in the drawing; the source model keeps every participant and bound.
They trigger structural mode when sampling is requested.

Reactions with one participant in the source model are included as boundary
reactions even without an EX label. DM/SK demand and sink subtypes are recorded;
unlabeled single-participant reactions are marked `exchange_candidate`.
Counting occurs before display filtering. Visible boundary reactions occupy
the perimeter and participate in orange connector discovery. Disconnected
ones are reported, with no fabricated connection. CellDesigner coordinates
and colors are not used to select pathways or position nodes.

## Example commands

### Draw a structural map without flux sampling

```bash
python sbml2map.py \
  --model models/e_coli_core.xml \
  --no-sampling \
  --coverage 100 \
  --output-dir results/runs/e_coli_core_structural
```

This displays all eligible reactions. Reactions with no visible participants
after display filtering are excluded. Arrows indicate directions permitted by
the source bounds: reversible reactions have arrows in both directions and
zero-bound reactions have no arrows. No flux or activity ranking is calculated.

A larger-model example:

```bash
python sbml2map.py \
  --model models/iML1515.xml \
  --no-sampling \
  --output-dir results/runs/iML1515_structural
```

To calculate structural maps for every included model in Bash:

```bash
for model in models/*.xml models/*.sbml; do
  [ -f "$model" ] || continue
  name=$(basename "$model")
  python sbml2map.py \
    --model "$model" \
    --no-sampling \
    --coverage 100 \
    --output-dir "results/runs/${name%.*}_structural"
done
```

These examples use `results/runs/` so the supplied result folders remain
available for comparison.

### View a drawing as an image

Open a PNG file to view the complete drawing.
Use the SVG version to zoom in without losing image quality.

| Model | PNG | SVG |
| --- | --- | --- |
| E. coli core | [PNG](results/e_coli_core/layout.png) | [SVG](results/e_coli_core/layout.svg) |
| iML1515 | [PNG](results/iML1515/layout.png) | [SVG](results/iML1515/layout.svg) |
| iMM904 | [PNG](results/iMM904/layout.png) | [SVG](results/iMM904/layout.svg) |
| Human1 | [PNG](results/Human-GEM/layout.png) | [SVG](results/Human-GEM/layout.svg) |
| iNESC2DN CMap | [PNG](results/iNESC2DN_CMap/layout.png) | [SVG](results/iNESC2DN_CMap/layout.svg) |
| iNESC2DN CTMap | [PNG](results/iNESC2DN_CTMap/layout.png) | [SVG](results/iNESC2DN_CTMap/layout.svg) |
| iNESC2DN ComboMap | [PNG](results/iNESC2DN_ComboMap/layout.png) | [SVG](results/iNESC2DN_ComboMap/layout.svg) |

The iNESC2DN drawings use SBML notes and CellDesigner identity information to
identify key pathways and cofactors. Red shows the matched pathways; orange
shows verified connections to exchange reactions. Unmatched compounds and
unreachable exchanges are listed in the result reports.

### Use a drawing in Escher

The initial Escher view has no reaction flux data loaded: sampled values, flux
widths, and the flux color scale are hidden. Arrows show model-permitted
directions; red and orange describe the pathway structure. Loading reaction
data automatically disables those fixed pathway colors, and clearing the data
restores them. FBA pages start with simulation values and a blue-to-red
absolute flux scale (gray for zero), with fixed pathway colors off.

Escher connections are exported as horizontal and vertical segments. Diagonal
incidences receive an intermediate `multimarker` at the corner; already aligned
connections remain one segment. These routing nodes do not add compounds or
reactions. Original node coordinates, stoichiometry, bounds, and flux values are
preserved. Arrowheads stay on terminal metabolite segments, and pathway colors
and purple reaction highlights apply to both parts of a routed connection.
Moving nodes manually in Escher does not automatically reroute them.

To show a biomass-maximizing FBA solution on the saved E. coli core layout:

```bash
python run_fba_escher.py --model models/e_coli_core.xml --layout-dir results/e_coli_core
```

Open `results/e_coli_core/escher_fba.html`. FBA uses the full source model and its
original bounds. Its dedicated `escher_fba_map.json` includes all 95 reactions
and all 360 reaction–metabolite connections, including water, protons and
cofactors as reaction-local copies. Zero-flux reactions remain visible in thin
gray. Existing pathway coordinates and sampled medians are preserved. All 95
source reaction values are available in `fba_fluxes.csv` and `fba_fluxes.json`;
`escher_fluxes.json` contains sign-adjusted values for all 95 reactions.
Label controls only change label visibility; the full FBA map has no coverage filter.
The bundled result has a biomass objective of approximately 0.873922 h⁻¹.
This is one optimal FBA solution; other optimal flux distributions may exist.

Both HTML views include a purple-highlight panel. Choose compounds or reactions,
search by ID or name, and click a result to select it. You can also add exact
model IDs in bulk, separated by commas or newlines. Compound selection colors
all displayed copies; reaction selection colors its lines, arrows, labels, and
grid-view reaction node. Purple overrides pathway colors. Clear selections to
restore the original palette.

In `index.html`, **Always show purple targets** overrides the display filters for
selected objects and shows the compounds connected to selected reactions in
their own colors. Objects omitted when the map was generated cannot be added;
their status is shown in the panel. Selection does not change layout or fluxes.
Save/load selections as JSON, or save the currently displayed SVG. Browser
changes do not update PNG files. To generate purple highlights in PNG, SVG, and
the initial HTML selection, run:

```bash
python sbml2map.py --model models/e_coli_core.xml --no-sampling \
  --highlight-metabolites pyr_c accoa_c --highlight-reactions PDH CS \
  --output-dir results/runs/e_coli_core_highlight
```

Alternatively, use `--highlight-file highlights.json`, with the same format as
the panel export: `{"metabolites": ["pyr_c"], "reactions": ["PDH"]}`.
Unknown IDs cause a CLI error. Known targets omitted from the drawing are
reported in `report.json` under `highlights.unavailable`.

Open the [E. coli core Escher preview](results/e_coli_core/escher_preview.html)
in a browser to explore the map.

To load this map into another Escher session, download
[escher_map.json](results/e_coli_core/escher_map.json)
and select it as the map file.

To display the sampled reaction values, also download
[escher_flux_medians.json](results/e_coli_core/escher_flux_medians.json)
and load it as reaction data.

These values are sampling estimates, not measured fluxes.
See the [report](results/e_coli_core/report.json) for diagnostics.

GitHub file pages do not run the HTML previews. Download the repository,
then open `results/e_coli_core/escher_preview.html` in a browser.
To use the interactive layout instead, open `results/e_coli_core/index.html`.

### Discover pathways without drawing

```bash
python sbml2map.py \
  --model models/iMM904.xml \
  --discover-only \
  --output-dir results/runs/iMM904_discovery
```

This writes `pathway_extraction.json` and `id_mapping.json` without sampling or
layout calculation.

### Draw a map using flux sampling

Sampling is enabled when `--no-sampling` is omitted:

```bash
python sbml2map.py \
  --model models/e_coli_core.xml \
  --samples 2000 \
  --chains 4 \
  --thinning 100 \
  --burn-in 10000 \
  --seed 0 \
  --output-dir results/runs/e_coli_core_sampled
```

The model's original bounds are used without growth-floor or objective
conditioning. Reaction selection prioritizes identified backbones and then
absolute median flux. Zero-median reactions are excluded. Sampling diagnostics
can fail to meet convergence thresholds; such output is explicitly provisional.
If sampling is unavailable or fails, the program automatically falls back to
structural mode and records the reason in `report.json`.

To reuse samples generated by the preceding command:

```bash
python sbml2map.py \
  --model models/e_coli_core.xml \
  --reuse-samples results/runs/e_coli_core_sampled \
  --coverage 50 \
  --output-dir results/runs/e_coli_core_sampled_50
```

Reuse requires the saved raw samples and an exactly matching parsed-model hash.
It cannot be combined with `--no-sampling`. Use a separate output directory for
structural runs rather than one containing sampling files.

### Inspect available options

```bash
python sbml2map.py --help
```

Useful options include `--coverage` (0–100), `--time-limit` (the grid MILP time
limit, not a total runtime limit), `--max-milp-nodes`, `--grid-cols`, `--grid-rows`,
`--glucose-id`, `--identity-map`, and `--classification`. Identity-map files map
source metabolite IDs to reference roles; classification files map IDs to
`inorganic`, `coenzyme`, or `organic`.

## Saved drawings: `results/`

This folder contains drawings for the seven included models. Some results use
flux sampling; others show structure only. The table lists the mode and the
number of displayed reactions for each model.

| Result folder | Mode | Displayed / source reactions |
| --- | --- | ---: |
| [e_coli_core](results/e_coli_core/index.html) | Sampled; mixing checks passed | 76 / 95 |
| [iML1515](results/iML1515/index.html) | Sampled; convergence checks not met | 1,635 / 2,712 |
| [iMM904](results/iMM904/index.html) | Sampled; convergence checks not met | 848 / 1,577 |
| [Human-GEM](results/Human-GEM/index.html) | Structural; no sampling | 12,737 / 12,971 |
| [iNESC2DN_CMap](results/iNESC2DN_CMap/index.html) | Structural; no sampling | 1,834 / 1,915 |
| [iNESC2DN_CTMap](results/iNESC2DN_CTMap/index.html) | Structural; no sampling | 1,976 / 2,056 |
| [iNESC2DN_ComboMap](results/iNESC2DN_ComboMap/index.html) | Structural; no sampling | 2,222 / 2,305 |

[results/manifest.json](results/manifest.json) lists the model files, calculation
settings, validation results, and file hashes. Raw flux samples for E. coli core
are available in [flux_samples.npz](results/e_coli_core/flux_samples.npz).
Other supplied sampled results do not include raw samples. To use
`--reuse-samples`, use the E. coli core samples or generate samples with the
sampling command above.

Typical files in each result folder are:

| File | Purpose |
| --- | --- |
| `index.html` | Interactive layout; open in a browser |
| `layout.png`, `layout.svg` | Complete static drawing |
| `orange_routes.png`, `orange_routes.svg` | Backbone and EX-connector drawing, when provided by that result version |
| `escher_preview.html` | Browser-based Escher viewer |
| `escher_map.json` | Escher-format map |
| `source_model.json` | Parsed source model, including original reaction participants and bounds |
| `display_model.json` | Display projection after filtering |
| `coordinates.csv`, `edges.csv` | Saved node coordinates and reaction incidences |
| `report.json` | Layout, sampling, validation, and limitation report |
| `pathway_extraction.json`, `id_mapping.json` | Pathway matches and supporting chemical-role mappings |
| `reaction_ranking.csv` | Reaction selection and ordering |
| `exchange_connectors.json` | EX connections and unreachable exchanges, when available |

The display projection omits classified inorganic participants and verified
source/sink symbols. It must not be used as an FBA model. The original parsed model is retained in
`source_model.json`.

## Display conventions and limitations

- Red identifies recognized primary pathways. In results with connector
  highlighting, orange identifies shared undirected connections to exchange
  reactions; other branches are gray. Sampled results may also use flux colors.
- Ordinary metabolite and reaction nodes occupy distinct integer-grid points.
  Local cofactor copies normally sit off-grid near their reactions. Cofactors
  required to connect pathways to exchanges can instead become shared grid nodes.
- Exchange reactions are placed on the outer boundary. Layout constraints can
  differ between saved results; see each `report.json`.
- Large layouts use heuristic placement. Global optimality is not assumed, and
  edge crossings and label overlap are not minimized. Dense maps are best
  inspected through their interactive or SVG views.
- Connectivity and permitted directions do not establish steady-state
  feasibility or physiological pathway activity.
- The iNESC2DN files retain their CellDesigner annotations. Chemical-role matching
  uses a limited reference catalog, so some compounds remain unresolved. The
  reports list matching evidence and exchange reactions with no verified path.
  Restoring annotations does not prove physiological activity.
- Human1 and iMM904 also have partial pathway recognition. Consult the individual
  extraction reports rather than interpreting unmatched pathways as absent.

## Flux sampling with SBML boundary species

Species explicitly marked `boundaryCondition=true` are excluded from the
steady-state balance equations. All original reaction columns and bounds are
retained; internal water, protons and cofactors remain balanced even when hidden
in the drawing. Compartment names and drawing-only placeholder flags do not
exempt species from balance. JSON accepts boolean `boundaryCondition` values
(default: false). No growth constraint is added.

Use `--require-sampling` to fail if sampling cannot be acquired or validated,
instead of automatically producing a structural drawing. This cannot be combined
with `--no-sampling`. Mixing diagnostics remain separate: check
`sampling_report.json` for `mixing_checks_passed`, split R-hat and ESS. Saved-sample
reuse verifies the balance policy, constraints fingerprint, model/sample hashes,
actual balance residuals and bounds. Legacy samples without boundary species
remain supported. Unsupported additional constraints and drawing placeholders
without an explicit boundary flag remain excluded from sampling support.

Sampled Escher previews provide a `Show reactions` slider (0–100% of saved
reactions). Red pathways come first, followed by orange connectors, then other reactions;
absolute sampled median flux determines order within each group. Attached nodes and labels follow the filter. Changing the
percentage does not resample fluxes or recompute coordinates.

## Orange route penalty

This option makes reactions with lower absolute median flux cost more when
choosing orange paths. Higher values give a stronger preference to higher-flux
reactions, even if the path becomes longer. The default is `1`; use `0` to
consider length only. A fractional value such as `0.5` is also allowed.

```bash
python sbml2map.py --model models/e_coli_core.xml \
  --orange-flux-weight 0.5 --output-dir results/runs/e_coli_core_weighted
```

The value must be zero or positive. The option is ignored with `--no-sampling`.

Use `--orange-flux-weight 1` to favor larger absolute sampled medians during
orange route selection and placement. The default is 1; use 0 for length only. Candidate
non-red reactions receive midrank percentile values q in [0, 1], and each unique
orange edge costs `length * (1 + weight * (1 - q))`. Shared edges count once.
`exchange_connectors.json` records physical `total_length`, weighted `total_cost`,
normalization values and reaction coefficients separately. Sampling-free runs
use the original unweighted objective. Routing, heuristic moves, MILP and
incumbent comparisons use the same objective; global optimality is not assumed.


## Avoid nodes with whole-segment rectangular detours

In structural SVG and Escher views, red/orange segments that cross an unrelated
red/orange metabolite or reaction node can move sideways as three sides of a rectangle.
The direct U uses two bends to connect the segment endpoints to the shifted line;
node positions, reaction participants and flux values stay unchanged. This also
applies to primary pathway lines. Both the long side and the two short sides
must be free of unrelated red/orange nodes. Gray and other nodes may be crossed.
The same structural node classification drives SVG colors and both viewers’
obstacle lists; purple selections and display filters do not change it.
Gray cofactor copies on red reactions remain outside the obstacle list.
If a direct U is obstructed, the line first
steps approximately 0.15 grid units along the original direction at one endpoint,
then turns twice to reach the offset lane. Both endpoints can use this fallback.
On short segments, the step is capped at one quarter of the segment length.
Existing short coincident endpoint stems may remain, without increasing their
overlap; new coincident paths are rejected. If every candidate is obstructed,
the original segment is retained and the reason is recorded.

```bash
python sbml2map.py --model models/e_coli_core.xml --no-sampling \
  --route-offset 0.3 --output-dir results/runs/e_coli_core_node_detours
```

`--route-offset` specifies the displacement from the original line in grid units
(default `0.3`). Use `0.5` for a wider detour or `0` to disable node avoidance.
Negative or non-finite values are rejected. Existing parallel-lane separation
is independent and remains at 0.3 grid units. The setting and unresolved cases
are saved under `node_detours` in `parallel_routing.json` for SVG and
`overlap_routing.node_detours` in Escher map metadata. This is display routing;
layout objectives still use the original routes.

Accepted node detours record the chosen side and endpoint steps in `detours`.
`endpoint_step_fallback_segments` counts segments that needed this fallback.
For E. coli core, the saved SVG and Escher routes from `for_e` to `EX_for_e`
now avoid `succ_c` with the default 0.3-grid displacement.

Shared trunks are grouped only when their collinear intervals overlap with
positive length. Sharing an endpoint alone does not bundle branches. This keeps
collision-free short branches straight, including `glc__D_e`–`GLCpts` in the
saved E. coli core SVG; long branches still detour when needed.


## Reserve red parallel channels during layout

`--red-parallel-spacing` sets the minimum centerline spacing for red segments
of different reactions that run in parallel for at least one grid unit. The
default is `1.0`; use `0` to disable this layout condition. Same-reaction shared
trunks and endpoint-only contacts are excluded. Negative/non-finite values fail.

```bash
python sbml2map.py --model models/e_coli_core.xml --no-sampling \
  --red-parallel-spacing 1 --route-offset 0.3 \
  --output-dir results/runs/e_coli_core_red_spacing
```

The grid layout stage reserves compatible routing channels before rendering;
SVG and Escher consume the same planned red paths. Fixed biological node
positions are preserved. When fixed endpoints cannot be separated by moving
nodes, routing channels can extend the drawing grid; `routing_envelope` records
this expansion. The layout objective remains the original node-incidence cost;
channel displacement is optimized separately, without a joint global-optimality
claim. Both MILP-based and heuristic grid layouts use the same geometric checks.
A bounded greedy allocation is tried first; the MILP path can use a bounded
finite-candidate solver if necessary. Each channel-allocation search has a budget of at most 3 seconds (or the smaller
layout time limit), plus final geometry validation. Candidate grid layouts can
invoke this allocation more than once; the placement report records the total
additional time as `red_parallel_layout_total_seconds`.
Unresolved assignments are reported, never marked as satisfying the condition.

`--route-offset` still controls ordinary red/orange node-avoidance detours.
Verified red channels are kept fixed by subsequent display routing. Reports
record the requested spacing, minimum actual spacing, violations, elapsed time
and selected channels under `red_parallel_layout` / `red_parallel_spacing`.


## Prefer simple important routes; let gray paths yield

Structural red/orange routing selects feasible candidates by bend count first,
then path length, then crossings/interference. Protected red/orange nodes and
verified red/red channel spacing remain constraints. Gray coincidences cannot
force extra bends on an important path: conflicting gray paths are rerouted
separately while important geometry stays fixed. Shared reaction trunks remain
bundled. Candidates include straight, L/Z/U and endpoint-step forms; unresolved
gray paths are recorded rather than changing important routes to accommodate
them. Preallocated red channels are preserved.

The saved E. coli core `EX_ac_e`–`ac_e` branch uses a two-bend U. The conflicting
gray PGL route yields. Diagnostics are `important_simplification` and `gray_yield`
in SVG reports and Escher `overlap_routing` metadata.


## Label visibility and readability

SVG, structural Escher, and FBA viewers expose the same three controls:
`Metabolite labels`, `Reaction labels`, and `Allow label overlap`. Hiding a
reaction label also hides its displayed flux value. Disabled label types remain
hidden when overlap is allowed, when filters/zoom change, and when selected
objects are forced visible. Important labels take priority during overlap
filtering. Reaction-label copies in Escher’s priority layer follow the same
visibility state. Saved current-view SVG exports retain visibility and halos.

White outlines separate important labels from red/orange lines without moving
nodes, routes, label coordinates or changing reaction data. Labels are drawn in
front of the lines. This first readability pass uses halos and overlap filtering.

FBA colors run from blue through violet and magenta to red, with square-root
scaling to distinguish low and moderate absolute fluxes. Widths range from 3 to
8; zero flux stays gray at width 1. Legend ticks and signed labels retain actual
flux units. Numerical FBA results and branch visibility are unchanged.
