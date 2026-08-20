# ONE-BLUE Project — DADA2 Amplicon Pipeline (16S & 18S, multi-run)

This README documents the parameter choices made while running the [dada2_guidelines](https://github.com/adriaaula/dada2_guidelines) pipeline on the `ONE-BLUE` metabarcoding datasets, and explains **why** each value was chosen. It's meant so anyone can rerun this analysis (or adapt it to a new run) without having to re-derive every decision from scratch.

Two amplicons were processed in parallel — 16S (prokaryotes) and 18S (eukaryotes) — each in its own project directory (`dada2_guidelines_16S/`, `dada2_guidelines_18S/`).

For 16S specifically, data came from **three separate sequencing runs**, each processed independently through Steps 0–2 before being combined:

| Run | Notes |
|---|---|
| **MetaLor2** | First run processed; `truncLen=225,175` |
| **MetaLor4** | `truncLen=220,190` |
| **MetaLor1** | `truncLen=220,190` |

(MetaLor3 was the corresponding 18S run and follows the 18S-specific parameters in this README; it was not part of the 16S multi-run merge.)

Sections 1–6 below describe the reasoning using the first (MetaLor2) run as the worked example. Section 7 covers how MetaLor1 and MetaLor4 were checked against that same reasoning, and how all three were combined at Step 3.

---

## 1. Primers and expected amplicon length

| | Forward primer | Reverse primer | Amplicon length (incl. primers) |
|---|---|---|---|
| **16S** | 515F-Y: `GTGYCAGCMGCCGCGGTAA` (19 nt) | 926R: `CCGYCAATTYMTTTRAGTTT` (20 nt) | ~400–420 bp |
| **18S** | TAReuk454FWD1: `CCAGCASCYGCGGTAATTCC` (20 nt) | V4RB: `ACTTTCGTTCTTGATYRR` (18 nt) | ~400–430 bp |

Reads provided to DADA2 were **already primer-trimmed** (see filenames: `*_trimmed_R1/R2.fastq.gz`), so the primer lengths above only matter for calculating the *expected biological insert length* (amplicon − primers), which is needed to check read overlap during merging:

- 16S insert (post-primer): **~361–381 bp**
- 18S insert (post-primer): **~362–392 bp**

---

## 2. Step 0 — Quality inspection (`Qscore` plots)

Before choosing `truncLen`, we inspected the per-position quality plots (`FastQC`/`Qscore`-style output) for both forward (R1) and reverse (R2) reads of each amplicon.

**What to look for:**
- Where the mean/median quality (green/orange line) drops below ~Q20–30
- The **maximum read length actually present** (don't assume the plot's x-axis max = every read's length)
- Reverse reads (R2) almost always degrade earlier/faster than forward (R1) — normal Illumina behavior

**Initial visual read of the plots** suggested:
- 16S: truncLen ≈ 240 (fwd) / 175 (rev)
- 18S: truncLen ≈ 240 (fwd) / 195 (rev)

⚠️ **This first estimate was wrong** — see Section 3.

---

## 3. Lesson learned: verify actual read lengths before trusting the quality plot

Running `01_run-dada2.sh` with `truncLen=240,195` (18S) and `truncLen=240,175` (16S) caused **100% of reads to be filtered out**, across every sample (including high-depth ones):

```
The filter removed all reads: ... not written.
Warning message: In filterAndTrim(...): No reads passed the filter.
```

**Root cause:** `filterAndTrim()` discards any read *shorter* than `truncLen` — it does not skip or pad. The quality plot's x-axis maximum is not necessarily where most (or any) reads actually reach. In our case, the primer-trimming step upstream had already shortened reads more than the plot suggested.

**Diagnosis — always check real read lengths directly**, e.g. with `seqkit stats` or:

```bash
zcat sample_R1.fastq.gz | awk 'NR%4==2{print length($0)}' | sort -n | uniq -c | tail -20
```

This revealed the true picture for 18S:

| | min_len | avg_len | max_len |
|---|---|---|---|
| R1 (fwd) | ~78–229 (outliers) | **230.0** | 235–238 |
| R2 (rev) | ~194–231 (outliers) | **232.0** | 234–238 |

i.e. reads were clustered tightly around 230–232 bp — well below the 240 bp we'd requested.

**Takeaway:** *Always cross-check the quality-plot-based truncLen guess against the actual max read length in your files (`seqkit stats` is fastest) before running `filterAndTrim`.*

---

## 4. Corrected Step 1 parameters (`01_run-dada2.sh`)

Given the real read lengths and the overlap requirement (`truncLen(fwd) + truncLen(rev) ≥ insert length + ≥20 bp overlap`):

| Amplicon | truncLen (fwd, rev) | maxEE (fwd, rev) | minOverlap | Resulting overlap margin |
|---|---|---|---|---|
| 16S | `225,175` | `2,6` | 15 | ~19–39 bp above minimum |
| 18S | `225,195` | `2,6` | 15 | comfortable (reads ~230bp, insert ~362–392bp worst case) |

```bash
# 16S
Rscript scripts/preprocessing/01_dada2-error-output.R \
    data/trimmed data/dada2/ blanes_project 225,175 2,6 15

# 18S
Rscript scripts/preprocessing/01_dada2-error-output.R \
    data/trimmed data/dada2/ blanes_project 225,195 2,6 15
```

**No pooling** was used (7th argument left blank) since samples were processed independently — this also determines the chimera-removal method in Step 2 (see Section 6).

---

## 5. Checking the error model plots

After Step 1 runs `learnErrors()`, inspect the generated `errors_<name>_fwd.pdf` / `_rev.pdf`.

**How to read them:**
- Each panel = one base transition (e.g. `A2C` = true A, observed as C)
- Black dots = observed error rate at each quality score; black line = fitted model DADA2 will use
- Red line = theoretical error rate implied by the nominal Q-score
- **Good sign:** error rate decreases smoothly as quality score increases, and the fitted black line tracks the point cloud (roughly parallel to, or a bit above, the red line)
- **Bad signs:** non-monotonic/erratic points, flat lines, or curves that don't decline at high Q — usually indicates too few reads passed filtering (i.e. go back and check your `truncLen`)

Both 16S and 18S error plots looked clean and well-behaved after the truncLen correction in Section 4 — confirming enough reads survived filtering to build a reliable model.

---

## 6. Step 2 — Chimera removal & length trimming (`02_chimerarem_merge.R`)

Step 1's output (`*_seqtab.rds`) still contains **chimeric sequences** (PCR artifacts from two different templates spliced together) and needs one more length-based cleanup pass before it's usable.

### 6a. Don't set the trim-length range from ASV counts alone

`table(nchar(colnames(seqtab)))` counts **distinct ASVs per length** — it does *not* tell you how many actual reads support each length. A length bin with few ASVs can still hold huge numbers of reads (real, abundant biology), and a bin with many ASVs can be almost entirely low-abundance noise. **Always check read-weighted length distribution before picking a cutoff:**

```r
fina <- readRDS("path/to/blanes_project_seqtab.rds")
tapply(colSums(fina), nchar(colnames(fina)), sum)
```

### 6b. What the read-weighted data showed

**16S:** reads piled up almost entirely in one tight, contiguous block (365–378 bp), with millions of reads in that range and only negligible counts (tens–hundreds) outside it.
→ **Trim range: `362,380`**

**18S:** reads showed a genuinely **multi-modal** distribution (real peaks around 327, 333, 374, 379–386, 390, 394 bp) — expected, since different eukaryotic lineages (diatoms, dinoflagellates, ciliates, etc.) have naturally different V4 lengths. A narrow cutoff here would have wrongly excluded whole taxonomic groups.
→ **Trim range: `320,401`**

One suspicious signal was flagged and excluded: a large, isolated spike of reads at exactly 225 bp (the `truncLen` floor) with no other nearby real peak — most likely a primer-dimer/non-specific product rather than true 18S signal. Confirmed by checking which sample(s) it came from:

```r
idx <- which(nchar(colnames(fina)) == 225)
rowSums(fina[, idx, drop = FALSE])
```

### 6c. Chimera removal method

**`consensus`** (the script default) — correct choice here because Step 1 was run **without** pooling. If pooling (`pool` or `pseudo`) had been used instead, the `pooled` method would be required.

### 6d. Commands run

```bash
# 16S
Rscript scripts/preprocessing/02_chimerarem_merge.R \
    data/dada2/01_errors-output/blanes_project/blanes_project_seqtab.rds \
    data/dada2/ \
    blanes_project \
    362,380 \
    consensus

# 18S
Rscript scripts/preprocessing/02_chimerarem_merge.R \
    data/dada2/01_errors-output/blanes_project/blanes_project_seqtab.rds \
    data/dada2/ \
    blanes_project \
    320,401 \
    consensus
```

---

## 7. Adding more 16S runs: MetaLor1 and MetaLor4

`03_merge` is used **only when you have multiple sequencing runs of the same marker gene** to combine into one table (as opposed to `02_`'s job, which is chimera removal + length filtering *within* a single run). For 16S, we had three runs to bring together — MetaLor2, MetaLor4, and MetaLor1 — so this step was actually needed here.

### 7a. Different `truncLen` per run is expected, not a problem

MetaLor4 and MetaLor1 were each run through Step 1 with their own `truncLen`, chosen from their own quality plots / actual read lengths, independently of MetaLor2:

```bash
Rscript scripts/preprocessing/01_dada2-error-output.R \
    data/trimmed data/dada2/ blanes_project 220,190 2,6 15
```

This is fine, and doesn't need to be "matched" to MetaLor2's `225,175`. The reason: `mergePairs()` reconstructs the **true biological insert length**, not `truncLen(fwd) + truncLen(rev)` — the overlapping region is trimmed away during merging, so:

```
merged_length = truncLen(fwd) + truncLen(rev) − overlap_length
```

As long as each run's `truncLen` sum clears the overlap requirement (`≥` insert length `+ ~20bp`, comfortably above `minOverlap`), the resulting merged sequences converge on the same length regardless of how the trimming was split between forward/reverse. Different runs/flow cells legitimately need different quality handling — that's exactly what per-run Step 1 processing is designed for.

### 7b. Confirming convergence before trusting this

Before merging, we validated that the "different truncLen is fine" reasoning actually held for our data — this is worth doing rather than assuming it, since a genuine amplicon-length or primer problem in one run would also show up as a mismatched peak. Read-weighted length distributions (Section 6a's `tapply` approach) were pulled for all three runs' post-Step-1 `seqtab.rds` files:

| Run | Dominant peak(s) |
|---|---|
| MetaLor2 | ~370–374 bp |
| MetaLor4 | 374 bp (4.26M reads at that single length) |
| MetaLor1 | 372–374 bp (1.5–2.6M reads) |

All three runs converge tightly on the **same real biological insert length window (~365–380 bp)**, despite different `truncLen` splits. This confirmed it was safe to proceed with per-run Step 2 + a Step 3 merge.

### 7c. Step 2 applied per run

Since all three runs share essentially the same true length distribution, the **same trim range** was used for all three at Step 2, for consistency (and so it's simple to describe in a methods section — "a single ASV length filter was applied uniformly across all sequencing runs prior to merging"):

```bash
# Run for MetaLor2, MetaLor4, and MetaLor1 individually, each with its own seqtab.rds path
Rscript scripts/preprocessing/02_chimerarem_merge.R \
    <run-specific path>/blanes_project_seqtab.rds \
    data/dada2/ \
    blanes_project \
    362,380 \
    consensus
```

### 7d. Step 3 — merging the three runs

Once each run has an independent, chimera-free, length-filtered `seqtab.rds`, `mergeSequenceTables()` combines them by matching identical ASV sequence strings across tables — which is exactly why the Section 7b convergence check mattered: if the runs weren't truly measuring the same insert length/region, this step would silently produce a bloated table of "unique" ASVs that are really the same organism trimmed differently.

```bash
Rscript scripts/preprocessing/03_merge-runs.R \
    <MetaLor2 seqtab path>,<MetaLor4 seqtab path>,<MetaLor1 seqtab path> \
    data/dada2/ \
    blanes_project_merged
```

(Adjust script name/arguments to whatever `03_` actually expects — check its header comments the same way we did for `01_` and `02_` before running.)

### 7e. Sanity checks after merging

- Confirm the merged table's ASV length distribution still sits in the expected 362–380 bp range (it should, since each input table was already filtered to that range).
- Compare per-sample read totals in the merged table against each run's own `*_track_analysis.tsv` — total reads per sample shouldn't change from the per-run Step 2 output, only the number of samples/columns grows.
- If downstream taxonomy assignment or diversity analysis shows an unexpectedly large number of near-identical ASVs (e.g., sequences differing by 1bp that are almost certainly the same variant), that's a sign worth double-checking that the per-run filtering ranges were applied consistently.

---

## 8. What's next

- After Step 2 (and Step 3, if merging multiple runs), check the read-tracking table (`*_track_analysis.tsv`) generated in Step 1 for each run — in particular:
  - `diff3` (merged/denoised ratio) — should be high given the overlap margins calculated per run; a big drop here means the merge step is failing more than expected for that specific run/sample
  - Compare **ASV counts vs. read counts** before/after chimera removal — a healthy run removes a modest share of *ASVs* but only a small share of *reads*, since chimeras tend to be low-abundance
- Low `diff3` values isolated to a handful of samples within an otherwise healthy run (rather than affecting the whole run) are more likely sample-specific (community composition, library quality) than a parameter problem — worth inspecting those samples' length distributions individually rather than re-tuning `truncLen` for the whole run.

---

## Key lessons for next time (TL;DR)

1. **Quality plots tell you where quality degrades — they don't reliably tell you the true max read length.** Always confirm actual read lengths (`seqkit stats` or an `awk` one-liner) before setting `truncLen`, especially on already-trimmed reads.
2. **`truncLen` too high = 100% of reads filtered**, silently, for every sample — this is the single most common cause of a "no reads passed the filter" error.
3. **Always calculate the overlap margin** (`truncLen(fwd) + truncLen(rev) − insert length`) and keep it comfortably above the `minOverlap` you're using (≥20 bp is a safe target, don't rely on the DADA2 default minimum of 12 bp).
4. **When choosing chimera/length-trim cutoffs, use read-weighted length distributions (`tapply(colSums(...), nchar(...), sum)`), not raw ASV-count tables** — ASV counts can make real, abundant, tightly-clustered biology look identical to rare noise, and vice versa.
5. **Isolated spikes at exactly your `truncLen` value are a red flag** — check which samples they come from before assuming they're real.
6. **`truncLen` does not need to match across multiple runs of the same marker gene before merging with `03_merge`.** `mergePairs()` reconstructs the true insert length regardless of the fwd/rev trim split, as long as each run's own overlap margin is sound — but always confirm this with a read-weighted length check per run before merging, rather than assuming it.







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

