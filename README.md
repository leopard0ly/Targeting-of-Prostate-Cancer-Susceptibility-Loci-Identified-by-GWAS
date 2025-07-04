# Targeting-of-Prostate-Cancer-Susceptibility-Loci-Identified-by-GWAS
Full CRISPR GWAS &amp; epigenomic validation pipeline (Prostate cancer)
Plain-text README (no command snippets).
Last updated : 19 Jun 2025  Requires R 4.1 or newer.

======================================================================
0) Folder map
goldlayer_gwas_pipeline.R – main script
README.md – this file
LICENSE.md
data\ – create this folder and place the external
data files listed in Section 2
results\ – created automatically; holds CSVs and plot

======================================================================

How to run in three sentences
======================================================================
Clone or download the repository.
Create the data folder and copy every external file named in
Section 2 into it.
Open R (or RStudio) and source « goldlayer_gwas_pipeline.R » – the
script installs all required packages by itself.

======================================================================
2) Required external data
A – GTEx prostate eQTL (version 8)
• Go to the adult GTEx QTL download page:
https://gtexportal.org/home/downloads/adult-gtex/qtl
• Download the archive GTEx_Analysis_v8_eQTL.tar (about 1.5 GB).
• Extract the archive; locate the file
« Prostate.v8.signif_variant_gene_pairs.txt.gz ».
• Decompress it and copy the resulting text file into the data
folder.
B – H3K27ac ChIP-seq files (VCaP prostate cell line)
• Primary download page (use this as the authoritative source):
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE105760
• In the “Supplementary files” section download the three hg38 files:
1. Raw peaks (.bed or .bed.gz)
2. IDR-replicated peaks (.bed.gz)
3. Signal track (−log10 p) (.bigWig)
• Decompress the peak files if they are gzipped and place all three
files in the data folder.
======================================================================
3) What the script does
Checks and installs every R / Bioconductor package needed.

Retrieves significant prostate-cancer SNPs (p < 5 × 10⁻⁸) from the
EBI-GWAS catalogue.

Fetches a ±50 bp hg38 sequence window for each SNP.

Designs every possible SpCas9 spacer and scores them with Azimuth.

Flags SNPs that overlap ENCODE cCRE elements and the GTEx prostate
eQTL list.

Keeps guides with Azimuth > 0.5 that hit the intersection
cCRE ∩ eQTL – these 305 guides form the “GoldLayer 305” set.

Validates GoldLayer guides against raw peaks, IDR peaks and bigWig
signal from the H3K27ac dataset.

Writes CSV result tables, a bar-plot and an RDS object into the
results folder.

======================================================================
4) Outputs generated
results\GoldLayer_305_fullEpi.csv – every guide with all flags
results\Epigenomic_validated_guides_IDR.csv – subset inside IDR peaks
results\epigenomic_barplot_full.png – layer summary chart
high_efficacy_gRNAs.rds – Azimuth > 0.5 guide set

======================================================================

======================================================================
5-bis) Adapting the pipeline to a different tissue / cell-type
======================================================================

You must update FOUR inputs, then re-run.

┌────────┬──────────────────────────────┬──────────────────────────────────────────┐
│ Step   │ What to change               │ Where / how                              │
├────────┼──────────────────────────────┼──────────────────────────────────────────┤
│ 1      │ GWAS phenotype               │ In R script: function                    │
│        │ (EFO string *or* EFO ID)     │ fetch_prostate_assocs()                  │
│        │ Example:                     │ • change  efo_trait = "lung carcinoma"   │
│        │lung carcinoma  → EFO:0001072 │      or  efo_id   = "EFO_0001072"        │
├────────┼──────────────────────────────┼──────────────────────────────────────────┤
│ 2      │ Single-tissue GTEx eQTL file │ data\<Tissue>.v8.signif_variant_gene_pairs│
│        │ Example: Lung → download     │ .txt.gz  from GTEx page and place in     │
│        │   Lung.v8.signif_variant…    │ data\  ; update eqtl_file path in script │
├────────┼──────────────────────────────┼──────────────────────────────────────────┤
│ 3      │ Cell-specific epigenomic set │ Replace the three H3K27ac files:         │
│        │ used for validation          │ • raw peaks (.bed)                       │
│        │ (H3K27ac or other mark)      │ • IDR peaks (.bed.gz)                    │
│        │                              │ • signal track (.bigWig)                 │
│        │ Obtain them from GEO or      │ Put them in data\ and keep file names    │
│        │ ENCODE for the *new* cell.   │ consistent; edit file paths in Section 7 │
├────────┼──────────────────────────────┼──────────────────────────────────────────┤
│ 4      │ (Optional) restrict cCRE set │ By default the script uses the global    │
│        │ to cell-type (if desired)    │ ENCODE cCRE combined file. If you need   │
│        │                              │ a cell-specific regulatory set, download │
│        │                              │ the relevant BED, place it in data\ and  │
│        │                              │ change the cCRE file path in Section 5.  │
└────────┴──────────────────────────────┴──────────────────────────────────────────┘

► After updating these inputs run the script again.  
  All downstream steps (Azimuth scoring, cCRE ∩ eQTL filter, epigenomic validation)
  use the new files automatically.

Reminder – ENCODE experiment lookup:
 • GEO search:  “H3K27ac <your cell> ChIP-seq GEO”  
 • ENCODE filter:  Assay = ChIP-seq, Target = H3K27ac, Organism = human, Biosample

======================================================================
======================================================================
6) Licence
--------
Academic-only License — Free for non-commercial academic research.
Commercial use requires explicit permission. Contact [admin@leopard.ly].

======================================================================
7) Citation
--------
GTEx Consortium (2020) — Science 369
ENCODE Project Consortium (2020) — Nature 583

For reference to the matched ENCODE experiment (VCaP cell line):
https://www.encodeproject.org/experiments/ENCSR597ULV/

======================================================================
End of README

