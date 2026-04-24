# MarkerRepo

MarkerRepo is a curated repository and annotation toolkit for cell type marker genes in single-cell data. It combines three things:

1. **A database** of YAML marker lists with structured metadata (organism, tissue, source, tags, reference genome), searchable and combinable into custom selections.
2. **An annotation engine** that assigns cell types to clusters by matching each cluster's differentially ranked genes against marker lists, weighting matches by both differential-expression rank and marker specificity.
3. **Cross-organism transfer** of marker lists via BioMart (Ensembl) or HomoloGene, so lists curated for one organism can be applied to another.

The package supports both **gene expression** (scRNA-seq, matched by symbol) and **genomic regions** (scATAC-seq, matched by interval overlap). Inputs are an `AnnData` object with clusters and ranked genes; outputs are cell type assignments in `adata.obs` together with per-cluster score tables.

## Quickstart

```bash
conda env create -f environment.yaml
conda activate marker-repo
pip install .
python -m ipykernel install --user --name=marker-repo
```

Open `notebooks/markerrepo.ipynb` and select the `marker-repo` kernel.

## Notebook Sections

| Section | What it does | Prerequisite |
|---|---|---|
| 0. Setup & Imports | Load packages, set `repo_path` | none |
| 1. Search & Combine Marker Lists | Search the marker database, export combined lists | none |
| 2. Cell Type Annotation | Annotate a clustered AnnData object | `test_data/adata_annotation.h5ad` |
| 3. Create Marker Lists (YAML) | Build YAML marker lists from a TSV file | `test_data/marker_list_blood.tsv` / `marker_list_brain.tsv` |

Each section is independently runnable.

## How Cell Type Annotation Works

The core annotation routine (`markerrepo.annotation.annot_ct`) scores every candidate cell type against every cluster and assigns the top-scoring cell type to each cluster.

**Inputs per cluster.** For each cluster, the algorithm takes the ranked genes stored in `adata.uns[rank_genes_column]` — a list of genes ordered by their differential-expression score (e.g. the output of `scanpy.tl.rank_genes_groups`).

**Per-marker weighting.** Each marker in a cell type's list carries an optional *ubiquitousness index* `ub_i ∈ (0, 1]` describing how broadly that gene is expressed across cell types. The per-marker weight is

```
ubiquity_weight = round(sqrt(1 / ub_i))
```

This down-weights pan-expressed genes and up-weights specific markers. Two-column marker lists without an index use a flat weight of `1.0` (`annotation.py:382–399`).

**Score for a (cluster, cell type) pair.** For each marker in the cell type that also appears in the cluster's ranked genes:

```
weighted_score = rank_score × ubiquity_weight
```

These are summed across all matching markers and normalized by the square root of the cell type's marker list size, to avoid bias toward cell types with very long marker lists:

```
score = sum(weighted_score for matching markers) / sqrt(num_markers)
```

A cell type is only considered if at least `min_hits` markers match (default 4) — this suppresses spurious assignments from single-gene coincidences (`annotation.py:661–684`).

**Cell type ranking.** Within each cluster, cell types are sorted first by `score` (descending) and then by the ratio `match_count / num_markers` (descending). The second key acts as a specificity tiebreak: between two cell types with similar total scores, the one where a larger fraction of its markers matched is preferred (`annotation.py:696–699`).

**Genomic regions (ATAC).** For region markers (e.g. `chr1:1000-2000`), the same scoring formula applies, but matches are determined by interval overlap using an `IntervalTree` rather than exact string equality. Per-marker intervals can be padded with `upstream_offset` / `downstream_offset` to extend into promoter regions (`annotation.py:491, 532`).

**Output.** For each cluster, the top-ranked cell type is written to `annotation.txt`, and the full ranking table (`cell_type, score, hits, num_markers, mean_ubiquity`) is written to `ranks/cluster_{id}`. The top assignment is then added to `adata.obs`.

## Marker List Database

Each marker list is a YAML file under `lists/` with two sections:

- **`metadata`** — id (UID), name, organism + taxonomy ID, marker type (`Genes` / `Genomic regions`), submitter, source, tags (tissue, disease, life stage, etc.), and reference genome (for regions). The allowed fields and whitelists are defined in `keys.yaml`.
- **`marker_list`** — a list of cell types, each with a `name` and an array of markers. Markers may carry an Ensembl ID (`"SYMBOL ENSEMBL_ID"`) so that downstream steps can resolve either identifier.

The database is loaded in parallel via `ProcessPoolExecutor` (`parsing.py`) and flattened into pandas DataFrames for searching. Key APIs:

- `guided_search(repo_path)` — interactive column-by-column search with `+`/`-` include/exclude prefixes
- `search_df(df, keywords)` — programmatic search with the same semantics
- `export_marker_list(df, ..., marker_id="symbol"|"ensembl")` — export a selection

**Export formats:**

| Format | Columns | Use case |
|---|---|---|
| `two_column` | marker, cell type | Custom annotation, SCSA input |
| `score` | marker, cell type, score | Scored output where `score` is each marker's prevalence across the selected lists (0 = most specific, 1 = most common), min-max scaled (`scoring.py:compare_marker_lists`) |
| `ui` | marker, cell type, score | Same shape as `score`, but values come from the PanGlaO ubiquitousness index directly — transferred via homology for non-human/mouse organisms (`scoring.py:update_scores`) |
| `panglao` | PanGlaO six-column format | Direct drop-in for tools that expect the PanGlaO schema |

## Cross-Organism Transfer

`markerrepo.homology.prepare_gene_transfer` transfers marker lists between organisms using one of two backends:

- **BioMart / Ensembl** — online query via `apybiomart`; broader organism coverage
- **HomoloGene** — offline lookup against a local `homologene.data` file; faster, smaller organism set

Use this when marker lists exist for one organism (e.g. human) but you need them for another (e.g. zebrafish), or let `create_marker_lists(..., force_homology=True)` invoke it automatically when no lists are available for the target organism.

## Contributing Marker Lists

To create and submit your own marker lists, see `dev/notebooks/submit_lists.ipynb`. It walks through metadata entry, whitelist-based validation, and submission.

## Additional Notebooks

The `dev/` directory contains notebooks for less common workflows:

- `dev/notebooks/guided_annotation.ipynb` — interactive step-by-step annotation
- `dev/notebooks/scoring.ipynb` — marker ubiquitousness scoring
- `dev/notebooks/homology.ipynb` — step-by-step cross-organism transfer via BioMart or HomoloGene
- `dev/examples/` — example analyses on published datasets
