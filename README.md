# oneblue
Processing metagenomics and metabarcoding data from the project ONEBLUE (https://one-blue.eu/).

# oneblue Project — DADA2 Amplicon Pipeline (16S & 18S)

This README documents the parameter choices made while running the [dada2_guidelines](https://github.com/adriaaula/dada2_guidelines) pipeline on the `oneblue` metabarcoding datasets, and explains **why** each value was chosen. It's meant so anyone can rerun this analysis (or adapt it to a new run) without having to re-derive every decision from scratch.

Two independent pipelines were run in parallel — one for 16S (prokaryotes), one for 18S (eukaryotes) — each in its own project directory (`dada2_guidelines_16S/`, `dada2_guidelines_18S/`).
The data was uploaded manually to MARBITS from the hard drives that Ulises received from MACROGEN. Raw and processed data can be found there:
```/mnt/smart/scratch/emm1/projects/one.blue/data/metagenomes/1_reads/metabarcoding_oneblue```

---

## 1. Primers and expected amplicon length

| | Forward primer | Reverse primer | Amplicon length (incl. primers) |
|---|---|---|---|
| **16S** | 515F-Y: `GTGYCAGCMGCCGCGGTAA` (19 nt) | 926R: `CCGYCAATTYMTTTRAGTTT` (20 nt) | ~400–420 bp |
| **18S** | TAReuk454FWD1: `CCAGCASCYGCGGTAATTCC` (20 nt) | V4RB: `ACTTTCGTTCTTGATYRR` (18 nt) | ~400–430 bp |

Reads provided to DADA2 were **primer-trimmed** with cutadapt so the primer lengths above only matter for calculating the *expected biological insert length* (amplicon − primers), which is needed to check read overlap during merging:

- 16S insert (post-primer): **~361–381 bp**
- 18S insert (post-primer): **~362–392 bp**

The script used is called *helper00_stats_and_cutadapt.sh*

---

## 2. Step 0 — Quality inspection (`Qscore` plots)

Before choosing `truncLen`, we inspect the per-position quality plots (`FastQC`/`Qscore`-style output) for both forward (R1) and reverse (R2) reads of each amplicon.
Running the script *00_run-qscore.sh* gives you this plot for the first 9 samples (R1 and R2).

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
    with parameters 225,175 2,6 15

# 18S
Rscript scripts/preprocessing/01_dada2-error-output.R \
    with parameters 225,195 2,6 15
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

## 7. What's next

- `03_merge` is **only needed if you have multiple sequencing runs of the same marker gene** to combine into one table. We had a single run per amplicon, so this step was skipped entirely.
- After Step 2, check the read-tracking table (`*_track_analysis.tsv`) generated in Step 1 (and any equivalent output from Step 2) — in particular:
  - `diff3` (merged/denoised ratio) — should be high given our generous overlap margins; a big drop here means the merge step is failing more than expected
  - Compare **ASV counts vs. read counts** before/after chimera removal — a healthy run removes a modest share of *ASVs* but only a small share of *reads*, since chimeras tend to be low-abundance

---

## Key lessons for next time (TL;DR)

1. **Quality plots tell you where quality degrades — they don't reliably tell you the true max read length.** Always confirm actual read lengths (`seqkit stats` or an `awk` one-liner) before setting `truncLen`, especially on already-trimmed reads.
2. **`truncLen` too high = 100% of reads filtered**, silently, for every sample — this is the single most common cause of a "no reads passed the filter" error.
3. **Always calculate the overlap margin** (`truncLen(fwd) + truncLen(rev) − insert length`) and keep it comfortably above the `minOverlap` you're using (≥20 bp is a safe target, don't rely on the DADA2 default minimum of 12 bp).
4. **When choosing chimera/length-trim cutoffs, use read-weighted length distributions (`tapply(colSums(...), nchar(...), sum)`), not raw ASV-count tables** — ASV counts can make real, abundant, tightly-clustered biology look identical to rare noise, and vice versa.
5. **Isolated spikes at exactly your `truncLen` value are a red flag** — check which samples they come from before assuming they're real.

## From DADA2 to phyloseq object

The pipeline of Adria and Aleix was run on MARBITS with default names so I renamed the files.  
From now on I worked on local (Ramiro's *oneblue* folder in MARBITS is already full)
```
blanes_project_seqtab_final.fasta => oneblue_16S_seqtab_final.fasta ; oneblue_18S_seqtab_final.fasta
blanes_project_seqtab_final.rds => oneblue_16S_seqtab_final.rds ; oneblue_18S_seqtab_final.rds
blanes_project_track_analysis_final.tsv => oneblue_16S_track_analysis_final.tsv ; oneblue_18S_track_analysis_final.tsv
```
I also *cleaned* the headers of the ASV representative sequences:
```
sed 's:;.*::' oneblue_16S_seqtab_final.fasta > oneblue_16S.fasta
sed 's:;.*::' oneblue_18S_seqtab_final.fasta > oneblue_18S.fasta
```

The phylogenetic trees were computed following recomendations from this article: https://doi.org/10.12688/f1000research.8986.2
Using tools *mafft* (align sequences) and *FastTree* (make a distance matrix and tree from the alignment).
```
On MARBITS (in tmux sessions these can take a couple of hours easily)
module load mafft
module load fasttree

mafft oneblue_16S.fasta > oneblue_16S.aln
FastTree -nt -gtr oneblue_16S.aln > oneblue_16S.tree

mafft oneblue_18S.fasta > oneblue_18S.aln
FastTree -nt -gtr oneblue_18S.aln > oneblue_18S.tree
```

Taxonomic classification was also performed on local (in tmux sessions too)
```
tmux
R

library(dada2);library(Biostrings)
# Eukaryotes
euk_db <- "~/Desktop/DATA/taxonomy/pr2_version_5.1.1_SSU_dada2.fasta.gz" # This DB was downloaded locally
euk_asv_file <- readDNAStringSet("oneblue_18S.fasta")
euk_tax <- assignTaxonomy(seqs = euk_asv_file, refFasta = euk_db, multithread = TRUE,minBoot = 80, tryRC = TRUE, taxLevels = c("Domain", "Supergroup", "Division", "Subdivision", "Class", "Order", "Family", "Genus", "Species"))
saveRDS(euk_tax, file = "petrimed_18S_taxo_MB80.RData")

library(dada2);library(Biostrings)
# Prokaryotes
prok_db <- "~/Desktop/DATA/taxonomy/silva_nr99_v138.2_toGenus_trainset.fa.gz" # This DB was downloaded locally
prok_db_sp <- "~/Desktop/DATA/taxonomy/silva_v138.2_assignSpecies.fa.gz"      # This DB was downloaded locally
prok_asv_file <- readDNAStringSet("~/Desktop/oneblue/DATA/METAB/oneblue_16S.fasta")
prok_tax <- assignTaxonomy(seqs = prok_asv_file, refFasta = prok_db, multithread = TRUE, minBoot = 80, tryRC = TRUE)
prok_tax_sp <- addSpecies(taxtab = prok_tax , refFasta = prok_db_sp) # this failed originally so we worked by chunks
prok_tax_sp_20000 <- addSpecies(taxtab = prok_tax[1:20000,] , refFasta = prok_db_sp)
prok_tax_sp_40000 <- addSpecies(taxtab = prok_tax[20001:40000,] , refFasta = prok_db_sp)
prok_tax_sp_60000 <- addSpecies(taxtab = prok_tax[40001:60000,] , refFasta = prok_db_sp)
prok_tax_sp_top <- addSpecies(taxtab = prok_tax[60001:nrow(prok_tax),] , refFasta = prok_db_sp)
prok_tax_sp<-as.data.frame(rbind(prok_tax_sp_20000, prok_tax_sp_40000, prok_tax_sp_60000, prok_tax_sp_top))
saveRDS(prok_tax_sp, file = "~/Desktop//oneblue/DATA/METAB/oneblue_16S_taxo_mb80.RData")
```

The final phyloseq object was built with the R script *phyloseq_metab_oneblue.R*.
It includes: 
1/ merging files into a phyloseq object
2/ filtering of ASVs considered as contaminations (based on ASV presence in the controls)
3/ filtering of ASVs based on taxonomy
4/ removing controls
5/ rarefying the samples to the lowest number of reads across samples

Have a look at the script for further comments and feel free to come discuss any step if something is left unclear.
Now have fun with the metabarcoding data!








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

