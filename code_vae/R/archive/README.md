# Archive

Results and scripts that are **not** used for the current version of the
paper, kept for the record. All the figures of the paper are produced from the
colSBM fit `Rdata/save_sbmfit/iidcolsbm_trim_round8.rds` (selected model:
`$best_fit`, 6 groups), computed by
`07b_fit_colsbm_proxy_graphs_replicate_results.R`, and from the bootstrap fits
of `07c_fit_colsbm_proxy_graphs_resampling.R`.

The paths inside the archived scripts were not updated: they point to the
locations the files had before being archived.

## Rdata/save_sbmfit/

| File | Content |
|---|---|
| `iidcolsbm_trim.rds` | colSBM fit of the first version of the paper (`estim = 'saintclair'` in the former versions of `08_analyse_behaviour.R` and `10_analyze_group.R`), made by `scripts/07_fit_colsbm_proxy_graphs.R` |
| `iidcolsbm_trim_round2.rds` ... `round7.rds` | intermediate fits made while reproducing the first fit; the model with 5 blocks of `round7` was `estim = 'nousQ5'` |
| `sbm_*`, `mlvsbm_*`, `picolsbm_*`, `deltacolsbm_*`, `colsbm_year.rds`, `iidcolsbm_month.rds` | other models and time aggregations (year, quarter, month) fitted in 2023, not used in the paper |

## Rdata/save_colsbmfit_subsample/

Bootstrap fits of the first version of the paper (`subnets_80/90`: sub-sampled
networks, `res_iidcolsbm_subnets_80/90`: fits). **They are wrong**: the saved
sub-networks contain an object named `net_BCq8_trim` (the full-data networks)
instead of `sub_net_BCq8_trim`, and the estimation loop used the object
`sub_net_BCq8_trim` left in the R session, so the fits of a given fraction
were all estimated on the same network collection. They are replaced by the
fits of `07c_fit_colsbm_proxy_graphs_resampling.R`
(`Rdata/save_colsbm_bootstrap/`).

## scripts/

| Script | Content |
|---|---|
| `07_fit_colsbm_proxy_graphs.R` | fit of the first version (`iidcolsbm_trim.rds`) |
| `07b_compare_with_earlier_fits.R` | former version of `07b_...R`: fit of `round8`, then comparison of its mesoscale structure and clustering with the first fit |
| `11_analyze_results_secondrun.R` | income per group, duplicate of `10_analyze_group.R` |
| `prospective/` | exploratory description of the groups of the first fit |

## images/, images2/

- `images/Meso_saintclair.png`: mesoscale structure and groups of the first
  fit (5 groups)
- `images/iidcolsbm_group_trim.png`: earlier rendering of Figure 10 with the
  current fit (6 groups), replaced by `images2/iidcolsbm_group_trim.png`
- `images2/Ammount_per_trip_per_group.png`: income per trip for the three fits
  'saintclair', 'nous' and 'nousQ5'
