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
