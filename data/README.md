# Data

These are the project files linked from the Module 4 project lesson, "Counts to results", unchanged,
so your numbers match the lesson. **Synthetic, made for teaching**: no real sample or person is
behind them.

## practice_counts.txt

featureCounts output: 12 genes (PRG001 to PRG012) × 4 samples, two control and two treated.
The first line is a comment starting with `#` that records the featureCounts command. Each gene's row
starts with six columns of information about the gene (`Geneid`, `Chr`, `Start`, `End`, `Strand`,
`Length`), followed by one count column per sample, headed by the name of the alignment file it came
from (`aligned/ctrl_1.bam` and so on).

Made from simulated paired-end reads on a simulated 60 kb practice chromosome, aligned with HISAT2 and
counted with featureCounts v2.0.6. Four genes were built to change after treatment. You can make the
same table yourself from the raw reads in the free course Omics, Bioinformatics and AI, Module 5.

## practice_samples.csv

The sample sheet: columns `sample` and `condition` (`control` or `treated`), one row per sample, in
the same order as the count columns.

**Licence.** Made by MultiOmics Intelligence Ltd for LSMA teaching; free to reuse under the MIT
licence of this repository.
