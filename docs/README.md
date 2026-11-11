# Interactive drawing results

The `docs/` directory is ready for GitHub Pages. In the repository's Settings → Pages,
select **Deploy from a branch**, branch **main**, directory **/docs**, and Save.
Pages must be enabled before the public site becomes available.

When Pages is enabled, the default viewer URLs are:

- MetNetComp: https://metnetcomp.github.io/sbml2escher/e_coli_core/
- ELIXIR: https://elixir2026-project22.github.io/sbml2map/e_coli_core/

Append `escher_preview.html` for the structural Escher view or `escher_fba.html`
for the FBA overlay. Initial structural views have no reaction flux data loaded;
FBA views use simulation values and no fixed pathway colors.

After generating new results, refresh these files from the package root:

```bash
python3 prepare_github_pages.py
```

Commit and push `docs/` to update an enabled Pages site. This folder publishes
the E. coli core and iCH360 views and related downloads, without raw samples or new model files.

The explicitly exported positioned SBML example is available as [layout.sbml](e_coli_core/layout.sbml),
with its [validation report](e_coli_core/sbml_layout_report.json). The CLI does not generate or update SBML by default.

## iCH360 (eco_ich_360)

- [Interactive SVG viewer](eco_ich_360/index.html)
- [Escher viewer](eco_ich_360/escher_preview.html)
- [PNG](eco_ich_360/layout.png)

Live site: https://metnetcomp.github.io/sbml2escher/eco_ich_360/

These saved structural views show 326 of the source model's 349 reactions,
without flux sampling. No new model files or source/display model snapshots
are included. `prepare_github_pages.py` refreshes E. coli core and preserves
this separately published iCH360 directory.

FBA color and width show absolute flux: gray at zero, then blue, light blue,
yellow, orange, and red as flux increases. Labels retain signed values.
User selections are highlighted in purple to distinguish them from low flux.

iCH360 publication excludes SBML/XML model files and links to them.
Run `python3 prepare_github_pages.py` before publishing to check this policy.

## Parallel red/orange routes

Highlighted segments from different reactions that run in parallel for at least
one grid unit and are less than 0.3 grid units apart receive rectangular detours.
Nodes stay fixed; primary backbone centerlines stay fixed for parallel-lane
separation and can detour around unrelated nodes. Same-reaction trunks remain bundled. SVG diagnostics are in each viewer folder’s `parallel_routing.json`;
Escher diagnostics are in `escher_map.json` under `overlap_routing.parallel_detours`.
Blocked lanes remain in place and are recorded. Display detours add length;
reported optimization costs still describe the original routes.


Node avoidance shifts a whole straight segment as three rectangle sides with
two bends by default, including on red primary pathways. The displacement is
`--route-offset` in grid units (default `0.3`; `0` disables node avoidance).
Diagnostics are under `node_detours` in the SVG report and under
`overlap_routing.node_detours` in Escher map metadata. Both requested sides
are checked; obstructed cases remain unchanged and are recorded.

If the direct U is obstructed, a 0.15-grid step near either or both endpoints
allows two right-angle turns into the offset lane. Existing short endpoint
stems may remain coincident without increasing their overlap. The chosen shape
is recorded in `detours`; `endpoint_step_fallback_segments` counts its uses.
The E. coli core `for_e`–`EX_for_e` branch now clears `succ_c` in both viewers.


Node avoidance protects only structurally red/orange metabolite and reaction
nodes. Gray/other nodes may be crossed. Cofactor copies on red reactions stay
gray unless explicitly used by an orange connector. SVG and Escher use the same
classification; selections and display filters do not alter routing. The policy
and protected source node IDs are recorded as `node_obstacle_policy` and
`obstacle_node_ids` in each SVG report and Escher `overlap_routing` metadata.

Shared trunks require positive-length overlap. Branches that only touch at an
endpoint are routed separately, keeping the E. coli core `glc__D_e`–`GLCpts`
branch straight while colliding long branches can detour.


Red/red runs of at least one grid unit now reserve channels during grid layout,
with `--red-parallel-spacing 1` by default (`0` disables the condition). SVG and
Escher share these paths. Fixed node positions remain unchanged; outer routing
channels may expand the drawing grid. Actual spacing and unresolved constraints
are recorded under `red_parallel_spacing`; layout planning and its elapsed time
are under `red_parallel_layout`. This option is independent of `--route-offset`.


Important red/orange routes now prefer fewer bends, then shorter length, among
feasible candidates. Gray lines yield instead of forcing extra important bends.
The E. coli core `EX_ac_e`–`ac_e` route is a simple two-bend U in both viewers;
the gray PGL route moves aside. Reports record simplification and unresolved
gray conflicts under `important_simplification` and `gray_yield`. Preallocated
one-grid red channels and protected-node avoidance remain in effect.


All SVG/Escher viewers, including the FBA page, have metabolite-label,
reaction-label and overlap controls. Important static labels and Escher labels
use white halos above the lines; label colors and coordinates are preserved.
Visibility is maintained across zoom/filter changes and selected-object forcing.
