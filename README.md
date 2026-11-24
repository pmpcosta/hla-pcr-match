# HLA PCR Match

[![CI](https://github.com/pmpcosta/hla-pcr-match/actions/workflows/ci.yml/badge.svg)](https://github.com/pmpcosta/hla-pcr-match/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/License-Apache_2.0-blue.svg)](LICENSE)

A fast, standalone Rust desktop GUI for **in-silico PCR primer matching** against the official **IPD-IMGT/HLA** human allele database (and any custom FASTA or IPD-MHC file), modeled after the [IPD-MHC Primer Match](https://www.ebi.ac.uk/ipd/mhc/primer/match/) web tool.

## Features

- **Exact IPD-MHC Primer Match form & semantics:**
  - **Specific sequences required** — comma- or newline-separated loci/alleles (e.g. `HLA-A`, `A`, `DRB1*15:01`, `B*07`, `MICA`), or leave empty to search all loci
  - **Forward primer (5'→3') & Reverse primer** — full IUPAC degeneracy codes (`R, Y, S, W, K, M, B, D, H, V, N`), live length/degeneracy readout, and live reverse-complement preview
  - **Type** — switch between `CDS database` (`hla_nuc.fasta`) and `Genomic database` (`hla_gen.fasta`, including introns & UTRs)
  - **Options** — `End Identity` (strict 3'-end exact match, default `3`), `Max mismatches` (default `3`), and `Flip reverse` (`Yes`/`No`, default `Yes`)
  - **Product Size** — `Min` / `Max` amplicon size bounds (`-1` = unbounded)
- **One-click IMGT/HLA Downloader & Local Cache:**
  - Downloads `release_version.txt`, `hla_nuc.fasta`, and `hla_gen.fasta.zip` directly from [ANHIG/IMGTHLA](https://github.com/ANHIG/IMGTHLA) and caches them locally for offline use.
  - Also supports loading any local `.fasta`, `.fa`, `.fna`, or `.zip` file.
- **Analysis & Export Extras:**
  - Parallelized 4-bit IUPAC bitmask matching (`rayon`) with virtualized, sortable, filterable tables (`egui_extras`)
  - Forward, reverse, and total mismatch counts per amplicon
  - Interactive color-coded **Binding Site Alignment Inspector** (3'-end identity bases highlighted in green, mismatches in red) + one-click **Copy Amplicon DNA**
  - **Locus Coverage** breakdown (% of alleles amplified per locus)
  - **Missed Alleles** diagnostic tab showing why each allele failed (forward unbound, reverse unbound, neither bound, invalid orientation, or product size out of range)
  - Export **Matched Hits (`.csv` / `.tsv`)**, **Missed Alleles (`.csv`)**, **Locus Coverage (`.csv`)**, and **Amplicon Sequences (`.fasta`)**

## Download

Pre-built standalone binaries for Windows, macOS, and Linux are published on the [GitHub Releases](https://github.com/pmpcosta/hla-pcr-match/releases) page whenever a version tag (`v*`) is pushed.

## Build from Source

Requires Rust (edition 2024 / Rust 1.85+).

```bash
cargo build --release
```

On Linux (Debian/Ubuntu), install the native GUI/file-dialog development packages first:

```bash
sudo apt-get update && sudo apt-get install -y \
  libgtk-3-dev libxcb-render0-dev libxcb-shape0-dev libxcb-xfixes0-dev libxkbcommon-dev
```

Run tests:

```bash
cargo test
```

## Data Attribution

Human HLA sequence files downloaded by this tool come from the official [IPD-IMGT/HLA Database](https://www.ebi.ac.uk/ipd/imgt/hla/) GitHub mirror ([ANHIG/IMGTHLA](https://github.com/ANHIG/IMGTHLA)), maintained by the Anthony Nolan HLA Informatics Group and EMBL-EBI. Please cite IPD-IMGT/HLA when using its data in publications.

## License

Licensed under the [Apache License, Version 2.0](LICENSE).
