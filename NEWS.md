# richCluster 2.0.0 (2026-09-15)

## Important: results differ from earlier versions

* This release fixes several problems in the similarity and clustering core,
  listed below. Results may differ from those of richCluster 1.0.2 and
  earlier, and we recommend re-running affected analyses with 2.0.0.
* `cluster_result` objects saved with an earlier version should be recomputed
  from the input enrichment tables rather than re-plotted. The bundled example
  object has been regenerated.
* There is no compatibility mode that restores the earlier behaviour; see
  "Reproducing results from richCluster 1.0.2" below.

## Fixes to similarity and clustering

* Kappa similarity now counts the genes shared by two terms correctly.
* Corrected the cluster merge logic so that qualifying clusters are merged as
  documented.
* `linkage_method = "single"` and `linkage_method = "complete"` now apply
  single and complete linkage respectively; earlier versions had the two
  exchanged.
* `linkage_method = "ward"` now implements Ward's minimum-variance criterion;
  earlier versions used average linkage for this option.
* `distance_metric = "jaccard"` now returns the Jaccard index,
  |A intersect B| / |A union B|.
* Kappa now returns 1 for two terms whose identical gene sets span the whole
  gene universe.
* `david_cluster()` now repeats its merge stage until no pair of final
  clusters exceeds the multiple-linkage threshold.

## Other changes in behaviour

* In `cluster()` and `runRichCluster()`, negative kappa values are set to 0, so
  kappa similarity lies on a [0, 1] scale. This non-negative similarity is a
  deliberate choice for clustering and differs from an unmodified Cohen's
  kappa. `david_cluster()` uses unmodified kappa, as the DAVID algorithm does.
* `cluster()` now selects terms on the adjusted p-value by default, matching
  its documentation. The new `filter_on` argument selects the column
  (`"Padj"` or `"Pvalue"`).
* `distance_cutoff` and `linkage_cutoff` are strict thresholds throughout: a
  score must exceed the cutoff to link or merge.
* Results are now reproducible across macOS, Linux and Windows: cluster
  numbering, term order and tie-breaking follow a canonical term order, and
  `merge_enrichment_results()` orders its rows independently of the locale.
* Cluster numbers are now contiguous and consistent with those used by the
  plotting functions. The representative term of a cluster is the term with
  the smallest value in the requested significance column.
* `full_network()` no longer includes self-loops, and `cluster_hmap()` now
  produces unique row labels.
* The diagonal of the exported similarity matrix is 1.
* `term_hmap()` now returns a plotly object drawn by heatmaply, like
  `cluster_hmap()`, and the iheatmapr dependency has been removed.

## New features

* `distance_metric = "dice"` adds the Dice coefficient.
* `cluster()` and `runRichCluster()` gain `gene_delim`, the separator used to
  split gene lists (default `","`).
* `cluster()`, `david_cluster()` and `runRichCluster()` gain `verbose`
  (default `FALSE`); progress output is silent unless `verbose = TRUE`.

## Diagnostics

* Warnings and errors carry stable condition classes:
  `richCluster_filter_on_fallback`, `richCluster_colname_collision`,
  `richCluster_duplicate_terms` and `richCluster_missing_geneid`.
* `merge_enrichment_results()` adds a `DatasetCount` column and warns when an
  input dataset repeats a `Term` value.
* `david_cluster()` validates its parameters, and an input too large for the
  similarity matrix (above about 46,000 terms) stops early with an
  informative error.

## Reproducing results from richCluster 1.0.2

* To reproduce an existing figure or table, install 1.0.2 in a separate
  library with `remotes::install_version("richCluster", "1.0.2")`. Use 2.0.0
  for new analyses.
