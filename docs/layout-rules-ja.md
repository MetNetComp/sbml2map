# Layout rules and priorities

[Latest circle and square views](test-tca/) · [Validation results](test-tca/validation.json)

The TCA cycle can be a circle or a square. Use `--tca-shape circle` (default)
or `--tca-shape square`, or choose **TCA shape** on the comparison page.
Only the TCA cycle follows clockwise order. Other pathways follow their own layout rules.

1. Keep the source metabolites, reactions, stoichiometry, and allowed directions.
2. Use the same size for nodes of the same type, including `pyr_c` and `accoa_c`.
3. Allow branches at a common node to overlap within 0.5 grid units of its center.
4. Place every node except TCA nodes and cofactors on integer grid intersections. Place cofactors near their reactions; they may be off the grid.
5. For all colors, forbid parallel or overlapping runs less than one grid unit apart that last for two or more grid units. Short auxiliary branches have no fixed count limit.
6. Keep the TCA cycle's order and chosen shape. The circular TCA has 18 nodes. Other nodes may be inside. Unrelated lines stay at least one grid unit from the outline. Real TCA connections may use short outward branches near their TCA node.
7. Place red main pathways first. Allow a red/orange point crossing only if avoiding it would at least double the length, add four grid units, or leave no valid route. Record the reason and candidate lengths.
8. Place confirmed exchange and transport chains in a straight line from outside inward, using real source species and reactions, including intermediate compartments.
9. Shorten right and bottom exchange tails where possible. Keep valid internal routes and recheck all saved routes.
10. Move static labels to avoid overlap. Fit margins around labels and arrows. SVG, PNG, and Escher use the same routes. Save static drawing checks in `render_validation.json`.

## Red-route priorities

Compare required rules first, then red length, changes to the saved red layout,
orange length, and margins. Initial placement uses Manhattan length for red
connections other than arcs. A MILP suggests TCA shifts. Route candidates,
including whole-red-layout shifts, and compare actual lengths with detours and Bezier arcs.

Check direct red-to-exchange routes, then save and fix the best red layout.
Try moving and rerouting other nodes, or expanding the boundary, first.
Reserve required straight routes before free routes. Only if that fails,
try equal-length red rerouting, equal-length red node moves, then longer red candidates.

`red_optimization.json` records initial and final lengths, moved nodes, changed
routes, fallback reasons, and failed candidates. Check red preservation after
margin trimming and in saved Escher data. The finite search does not prove a global minimum.

## E. coli core results

| TCA shape | Previous red length | Current red length | Red nodes moved | Red routes changed | Red length added after selection |
| --- | ---: | ---: | ---: | ---: | ---: |
| Circle | 181.53 | 151.53 | 0 | 0 | 0 |
| Square | 124.00 | 90.50 | 0 | 0 | 0 |

| TCA shape | Reactions | Long close parallel runs | Label overlaps | Boundary `[left,right,bottom,top]` |
| --- | ---: | ---: | ---: | --- |
| Circle | 84 | 0 | 0 | `[0, 33, 0, 46]` |
| Square | 84 | 0 | 0 | `[0, 30, 3, 46]` |

Each view has 220 labels. All 116 ordinary nodes outside the TCA and cofactor
groups are on integer grid intersections. Saved Escher routes and SVG node sizes
were checked. The 89 related tests passed.

The default is `--crossing-policy adaptive-grid`. Placement uses a MILP for some
steps; routing uses A* search and heuristics. The search does not prove the
shortest total routes or smallest area. Earlier FBA and SBML Layout examples
use separate layouts.

## iCH360 results (2026-11-13)

[Circle and square comparison](test-tca/?model=eco_ich_360) · [Validation results](test-tca/validation-eco_ich_360.json)

Each view shows 326 of 349 reactions, with 921 nodes and 1,073 connections.
Red length is 100.03 grid units for the circle and 96.50 for the square.
Both keep the selected red nodes and routes unchanged, with no long close
parallel runs and no label overlaps.

The layout uses real extracellular → periplasm → cytosol chains: 17 multi-step
transport groups and two simple exchange pairs, with 87 straight segments.
No intermediate reactions are added, and the input model is unchanged.
