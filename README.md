# Blanes/Oneblue Project — DADA2 Amplicon Pipeline (16S & 18S, multi-run)

This README documents the full parameter-selection process, diagnostic tests, and findings from running the [dada2_guidelines](https://github.com/adriaaula/dada2_guidelines) pipeline on the 16S and 18S datasets. It's written so anyone in the group can rerun this analysis, understand *why* each decision was made, and avoid re-running experiments we've already done (some of which failed informatively).

Two amplicons were processed — 16S (prokaryotes) and 18S (eukaryotes) — each in its own project directory. The 16S side involved **three separate sequencing runs** (MetaLor1, MetaLor2, MetaLor4) later combined into one project (`oneblue_16S_project`). MetaLor3 was the corresponding 18S run.

---

## 1. Primers and expected amplicon length

| | Forward primer | Reverse primer | Amplicon length (incl. primers) |
|---|---|---|---|
| **16S** | 515F-Y: `GTGYCAGCMGCCGCGGTAA` (19 nt) | 926R: `CCGYCAATTYMTTTRAGTTT` (20 nt) | ~400–420 bp |
| **18S** | TAReuk454FWD1: `CCAGCASCYGCGGTAATTCC` (20 nt) | V4RB: `ACTTTCGTTCTTGATYRR` (18 nt) | ~400–430 bp |

Reads provided to DADA2 were already primer-trimmed, so primer lengths only matter for calculating the *expected biological insert length* (amplicon − primers) used to check read overlap during merging:

- 16S insert (post-primer): **~361–381 bp**
- 18S insert (post-primer): **~362–392 bp**

---

## 2. Step 0 — Quality inspection, and a costly early mistake

Before choosing `truncLen`, we inspected per-position quality plots for forward (R1) and reverse (R2) reads.

**Initial visual read of the plots** suggested truncLen ≈ 240/195 (18S) and 240/175 (16S).

**This was wrong**, and it cost a full failed run: `filterAndTrim()` discards any read *shorter* than `truncLen` — it does not skip or pad. Running with `truncLen=240,195` filtered out **100% of reads, in every sample**, including high-depth ones. The quality plot's x-axis maximum is not the same as the actual max read length present in the file (especially after upstream primer trimming, which shortens reads further than the plot suggests).

**Lesson:** always confirm real read lengths directly before trusting a quality-plot-based truncLen guess:
```bash
zcat sample_R1.fastq.gz | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c | tail -20
# or, for many files at once:
seqkit stats *.fastq.gz
```

---

## 3. Corrected Step 1 parameters (`01_run-dada2.sh`)

Using real read lengths and the overlap requirement (`truncLen(fwd) + truncLen(rev) ≥ insert length + ~20bp overlap`):

| Run | truncLen (fwd, rev) | maxEE | minOverlap |
|---|---|---|---|
| MetaLor2 (16S) | `225,175` | `2,6` | 15 |
| MetaLor4 (16S) | `220,190` | `2,6` | 15 |
| MetaLor1 (16S) | `220,190` | `2,6` | 15 |
| 18S | `225,195` | `2,6` | 15 |

No pooling was used in any run (independent sample inference), which also determined the chimera-removal method choice later (Section 6).

**Key realization (Section 7): `truncLen` does not need to match across runs.** `mergePairs()` reconstructs the true biological insert length regardless of how trimming is split between fwd/rev — the overlapping region is trimmed away during merging. As long as each run's own overlap margin clears `minOverlap` with a healthy buffer, different runs can (and often must) use different truncLen values.

---

## 4. Checking the error model plots

After Step 1's `learnErrors()`, inspect `errors_<name>_fwd.pdf` / `_rev.pdf`.

**How to read them:** each panel = one base transition (e.g. `A2C`). Black dots = observed error rate per quality score; black line = fitted model; red line = theoretical nominal-Q error rate. **Good sign:** smooth monotonic decline as Q increases, fitted line tracking the point cloud. **Bad sign:** erratic/flat curves — usually means too few reads survived filtering (check truncLen again).

All runs, once truncLen was corrected, showed clean, well-behaved error plots.

---

## 5. Step 2 — Chimera removal & length trimming, per run

### 5a. Don't set the length-trim range from raw ASV counts

`table(nchar(colnames(seqtab)))` counts **distinct ASVs per length**, not reads — a bin with few ASVs can still hold huge read counts (real, abundant biology), and vice versa. Always check the **read-weighted** distribution instead:

```r
fina <- readRDS("path/to/seqtab.rds")
tapply(colSums(fina), nchar(colnames(fina)), sum)
```

### 5b. What the read-weighted data showed, across all three 16S runs

All three runs — despite different truncLen splits — converged tightly on the same real biological insert length window:

| Run | Dominant peak(s) |
|---|---|
| MetaLor2 | ~370–374 bp |
| MetaLor4 | 374 bp (4.26M reads at that length alone) |
| MetaLor1 | 372–374 bp (1.5–2.6M reads) |

This convergence was itself a useful validation: it confirmed all three runs were measuring the same amplicon and would merge cleanly.

**Final length-trim range used for all three 16S runs: `362,380`** (18S: `320,401`, given genuine multi-modal length variation across eukaryotic lineages — see Section 5c).

### 5c. 18S needed a wider window — real biology, not noise

18S read-weighted data showed **multiple real peaks** (327, 333, 374, 379–386, 390, 394 bp) — expected, since different eukaryotic lineages (diatoms, dinoflagellates, ciliates, etc.) have genuinely different V4 lengths. A narrow cutoff would have wrongly excluded whole taxonomic groups. One isolated spike (35,751 reads at exactly 225bp, the truncLen floor, with no neighboring real peak) was identified as likely primer-dimer/non-specific product and excluded.

### 5d. Chimera removal method

`consensus` was used initially (correct default given no pooling was used in Step 1). Later, `per-sample` was tested for the merged 16S table — see Section 8.

---

## 6. Multi-run merging (16S: MetaLor1 + MetaLor2 + MetaLor4)

### 6a. `02_chimerarem_merge.R` handles multi-run merging internally

Passing a comma-separated list of `seqtab.rds` paths as argument [1] merges the tables first, then runs chimera detection + length filtering on the combined table in one pass — no separate step needed:

```bash
Rscript scripts/preprocessing/02_chimerarem_merge.R \
    data/dada2/01_errors-output/all_16S/rn_metalor1_seqtab.rds,data/dada2/01_errors-output/all_16S/rn_metalor2_seqtab.rds,data/dada2/01_errors-output/all_16S/rn_metalor4_seqtab.rds \
    data/dada2/ \
    oneblue_16S_project \
    362,380 \
    consensus
```

### 6b. Duplicate control names across runs had to be resolved before merging

`BLANK`, `BPCR5`, etc. were named identically across MetaLor1/2/4 (same control naming convention reused per run). Since `seqtab` rows are indexed by sample name, merging as-is would have silently collapsed distinct controls from different runs into one row. **Fix:** rename only the colliding rownames per run before merging:

```r
rn <- rownames(ml1)[rownames(ml1) %in% rownames(ml2)]
rownames(ml1)[rownames(ml1) %in% rn] <- paste0("MetaLor1_", rownames(ml1)[rownames(ml1) %in% rn])
rownames(ml2)[rownames(ml2) %in% rn] <- paste0("MetaLor2_", rownames(ml2)[rownames(ml2) %in% rn])
rownames(ml4)[rownames(ml4) %in% rn] <- paste0("MetaLor4_", rownames(ml4)[rownames(ml4) %in% rn])
saveRDS(ml1, "metalor1_seqtab.rds")
saveRDS(ml2, "metalor2_seqtab.rds")
saveRDS(ml4, "metalor4_seqtab.rds")
```
**Only overlapping names were touched** — samples with unique names across all three runs kept their original names, avoiding unnecessary renaming.

### 6c. Validating that combining runs before chimera removal doesn't distort results

Before trusting the merged output, we directly compared per-sample chimera-removal results **combined vs. run separately** for the same samples. Deltas were within ±0.1–2.4% (noise-level) across every sample checked, for both MetaLor1 and MetaLor2 subsets. **Conclusion: merging runs before chimera detection does not distort per-sample outcomes** — `removeBimeraDenovo()`'s consensus method is robust to being run jointly vs. separately here.

---

## 7. Investigating uneven read retention — the "40% loss" samples

A subset of samples consistently retained far fewer reads through merging/chimera removal than others (`diff.total` as low as 0.46–0.59 vs. 0.85–0.98 for the best samples). This was investigated at length, with both a hypothesis-testing phase and a parameter-tuning phase.

### 7a. Ruling out a length/community-shift explanation

Pulling the length distribution for one of the worst offenders (`OBM21LSDNA`, `diff3=0.54`) showed the **same** dominant 372–378bp peak as healthy samples — no shifted or secondary peak. This ruled out "different dominant taxon with a different amplicon length" as the explanation, and pointed toward something depressing merge success within the correct length window (e.g., quality/error-driven), not a length-filter or composition issue.

### 7b. Two truncLen/minOverlap experiments — both confirmed the parameters were already near-optimal

- **Test 1 — tighter margin:** `truncLen=210,175` (sum=385, ~4bp overlap margin, minOverlap=15) → **catastrophic, uniform collapse**: every sample's `diff.total` fell to 0.02–0.35, including previously-perfect samples. This confirmed the overlap margin had been cut below what `minOverlap=15` requires — global failure, not a targeted fix.
- **Test 2 — looser margin + relaxed minOverlap:** `truncLen=220,180`, `minOverlap=12` (sum=400, ~19bp margin) → **uniform +0.01 improvement across every sample**, regardless of whether it was a "bad" or "good" sample. A real fix for the low-diff3 samples would show a *disproportionate* gain for them; a flat, universal +0.01 instead confirmed this parameter space wasn't where the problem lived. **Not adopted** — added complexity for no meaningful benefit.

**Conclusion: `truncLen`/`minOverlap` tuning cannot fix this gap.** Kept `220,190`/15 (and `225,175`/15 for MetaLor2) as final.

### 7c. Finding the real driver: depth-associated community complexity

Metadata correlation revealed the actual pattern. Initially flagged as a **sediment Littoral (L) vs. Offshore (O) naming split** (all "O" samples: `diff3` 0.85–0.98; most "L" samples: 0.46–0.72) — but this didn't hold for non-oneblue campaign samples lacking the L/O suffix, so **depth itself** was tested directly:

- **Sediment:** Pearson r=0.343 (p=0.003) between depth and `diff.total`, but not a smooth trend (Spearman not significant) — instead, a **threshold effect**: every sample >800m depth scored consistently well (0.81–0.93), while shallow (<150m) samples spanned almost the entire range (0.46–0.93), depth alone not explaining the shallow-sample spread.
- **Water column confirmation:** a 4-depth series from one station (`OBA07`, water samples at 5m/13m/1000m/3016m) showed the same pattern cleanly: `diff.total` = 0.65 / 0.64 / 0.76 / 0.82 — shallow consistently worse, deep consistently better, in the same station, same day, ruling out cross-site/batch confounds.

**Biological interpretation:** shallow water and shallow sediment both host far more diverse, complex microbial communities (light, nutrients, mixing, biological activity) than the deep ocean. More diverse, closely-related co-amplifying 16S variants during PCR → more sequences flagged as explainable-as-chimeras by `removeBimeraDenovo()`, even when many are true, distinct organisms. This is a known limitation of *de novo* chimera detection in high-diversity samples, not a pipeline defect.

### 7d. Confirming the loss happens at chimera detection, not length filtering

Directly quantified where reads were actually being lost, per sample:

| Sample | Lost to chimera removal | Lost to length filter (362–380) |
|---|---|---|
| OBA07OaWDNA (5m) | 35.14% | 0.033% |
| EX251CA1t4DNA | 49.16% | 0.019% |
| OBA36LaWDNA | 38.81% | 0.010% |
| NADEL302WDNA | 34.03% | 0.066% |

**The length-filter step (362–380bp) removes at most ~0.07% of reads in every case checked** — it is not the bottleneck. Virtually all loss happens earlier, at `removeBimeraDenovo()`. **This means widening the length filter (e.g., to 340,400) would not help** — it would only rescue that negligible ~0.03–0.07%, not the 20–50% actually lost at the chimera step. Tested and confirmed via direct calculation rather than assumed.

### 7e. Bottom line on read-yield variance

This is a real, reproducible, depth/diversity-driven feature of the dataset — validated across sediment and water samples, across all three runs, and across multiple parameter combinations (truncLen, minOverlap, length-filter width). **No further parameter tuning in Step 1/2 is expected to close this gap.** The one untested lever is the chimera-removal *method* itself (Section 8).

---

## 8. Step 2 rerun with `method="per-sample"` (in progress)

Rationale: `consensus` calls chimeras using a vote across all samples in the table; `per-sample` evaluates each sample independently, which can be less aggressive (retain more true positives) in high-diversity samples at the cost of being more permissive toward real chimeras. Since `consensus` was confirmed as the source of the shallow-sample read loss (Section 7d), this is the one untested lever left.

**Resource note:** the `per-sample` run initially OOM'd at 200GB memory; resubmitted with **250GB**. The default SLURM header in `02_run-chimerarem_merge.sh` only requests `--cpus-per-task=2` — worth revisiting alongside memory if runtime becomes an issue, since `per-sample` mode processes each sample's bimera detection independently and may parallelize differently than `consensus`.

*(Results pending — update this section once the run completes and compare `diff.total` for the low-yield shallow samples against the `consensus` baseline, plus confirm the well-performing deep samples don't regress.)*

---

## 9. Step 4 — ASV clustering (in progress, run in parallel with Section 8)

To further refine the final merged ASV table, `04_ASV-clustering.R` was run at 100% identity (i.e., collapsing only identical/redundant sequences at this stage, not true OTU-level clustering) as an additional QC pass on top of `02_`'s output:

```bash
Rscript scripts/preprocessing/04_ASV-clustering.R \
    data/dada2/02_nochimera_mergeruns/oneblue_16S_project/oneblue_16S_project_seqtab_final.rds \
    100 \
    data/dada2/ \
    oneblue_16S_project \
    abundance \
    0.9
```
- Clustering identity: **100** (0–100 scale) — currently used to check for redundant sequences, not for taxonomic-level OTU clustering.
- Representative method: **abundance** (most abundant ASV in a cluster becomes the representative, as opposed to the longest-sequence method).
- Min coverage: **0.9** — sequences must overlap the representative over ≥90% of positions to be clustered together.

*(Results pending.)*

---

## 10. What's next

- Compare `per-sample` (Section 8) vs. `consensus` chimera removal results once complete — check specifically whether the low-diff3 shallow samples improve without regressing the deep/high-diff3 samples.
- Review Step 4 clustering output — confirm the 283K raw ASV count (see below) is reduced sensibly without collapsing genuinely distinct organisms.
- Before final analysis, check the ASV table's **prevalence/singleton structure**:
  ```r
  prevalence <- colSums(seqtab_final > 0)
  table(cut(prevalence, breaks=c(0,1,2,5,10,50,Inf)))
  singleton_asvs <- prevalence == 1
  sum(seqtab_final[, singleton_asvs]) / sum(seqtab_final)
  ```
  A large raw ASV count (~283K, prior to Section 9 clustering) is not inherently alarming given the dataset spans sediment and seawater across seven regions and depths from 8m to 4350m — but confirm the bulk of that count is low-abundance singleton/doubleton ASVs (normal, and handled by standard downstream prevalence filtering) rather than genuine unexplained inflation.

---

## Key lessons for next time (TL;DR)

1. **Quality plots tell you where quality degrades — they don't reliably tell you the true max read length**, especially on already primer-trimmed reads. Always confirm actual read lengths (`seqkit stats` or an `awk` one-liner) before setting `truncLen`.
2. **`truncLen` too high = 100% of reads filtered**, silently, for every sample — the single most common cause of a "no reads passed the filter" error.
3. **Always calculate the overlap margin** (`truncLen(fwd) + truncLen(rev) − insert length`) and keep it comfortably above `minOverlap` (≥20bp target, not just the DADA2-default 12bp floor). Cutting this margin too tight causes global, catastrophic merge failure — not a graceful degradation.
4. **When choosing chimera/length-trim cutoffs, use read-weighted length distributions**, not raw ASV-count tables — ASV counts can make real, abundant, tightly-clustered biology look identical to rare noise, and vice versa.
5. **Isolated spikes at exactly your `truncLen` value are a red flag** for primer-dimer/non-specific artifacts — check which samples they come from.
6. **`truncLen` does not need to match across multiple runs of the same marker gene before merging.** `mergePairs()` reconstructs the true insert length regardless of the fwd/rev trim split, as long as each run's own overlap margin is sound — confirmed empirically via read-weighted length convergence across MetaLor1/2/4, not just assumed.
7. **Rename colliding sample/control names before merging seqtabs from multiple runs** — identical control names (e.g. `BLANK`, `BPCR5`) reused across runs will silently collapse into one row otherwise.
8. **Uneven per-sample read retention is not always a parameter problem.** Before re-tuning `truncLen`/`minOverlap`/length-filter width, check: (a) whether the ASV length distribution shifts for low-yield samples (composition change vs. quality issue), (b) whether loss happens at the chimera-removal step or the length-filter step specifically (quantify both, don't assume), and (c) whether metadata (depth, site, sample type) correlates with the pattern. In this dataset, depth-driven community diversity — not a fixable pipeline parameter — explained the majority of the variance, confirmed across both sediment and water-column samples.
9. **A parameter change that improves *every* sample uniformly (including already-good ones) is a sign you haven't found the real lever** — a genuine fix for a specific problem should show a disproportionate effect on the affected samples.





## 1/ Import on MARBITS
For now I created a project folder in my <emm1> account on MARBITs:
```
/mnt/smart/scratch/emm1/users/pramond/oneblue
```

Metagenomics data were processed by MACROGEN. A first batch of samples (seawater samples, 26) was received June 18th, 2026.
The transfer was performed directly from the hard drive that MACROGEN sent us to MARBITS ("only" 500GB).
We checked that the transfer was made correctly with their md5sum:
```
cd /mnt/smart/scratch/emm1/users/pramond/oneblue/MACROGEN/HN00276367/0.RawData
md5sum *.fastq.gz > marbits_md5
```

