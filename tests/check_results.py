"""Checks for the Module 4 project: counts to results.

Run from the repository folder:  python tests/check_results.py
Python 3 standard library only. It reads the files you commit in scripts/
and results/; it does not run R or DESeq2. Fold changes are checked with a
tolerance, so small differences between DESeq2 versions do not matter.
"""

import csv
import hashlib
import math
import re
import struct
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Fingerprint of the untouched starter script, so the check can tell
# whether you have started work on it.
STARTER = {
    "scripts/counts_to_results.R": "4d169bbfafb378fdc8f0ad18a59b1457fc8bb354ac95f18827c3c67fb01e5201",
}

LIBRARY = {"ctrl_1": 8998, "ctrl_2": 10999, "treat_1": 10000, "treat_2": 8500}
BUILT_IN = {"PRG003": "up", "PRG004": "up", "PRG007": "down", "PRG010": "up"}
TRUTH = {"PRG003": 4.0, "PRG004": 3.0, "PRG010": 2.0, "PRG007": 0.25}
# log2 fold changes from DESeq2 1.48 with design ~ condition (unshrunk); apeglm-shrunk values also pass
LFC = {"PRG003": 1.503, "PRG004": 1.238, "PRG007": -1.980, "PRG010": 0.683}
CPM = {"PRG003": 2.302, "PRG004": 1.925, "PRG010": 1.312, "PRG007": 0.204}


class Fail(Exception):
    pass


def fingerprint(path):
    text = path.read_text(encoding="utf-8", errors="replace").replace("\r\n", "\n")
    text = "\n".join(line.rstrip() for line in text.split("\n")).strip()
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def need_edited(rel):
    path = ROOT / rel
    if not path.is_file():
        raise Fail(f"{rel} is missing; keep the file and write your code in it")
    if fingerprint(path) == STARTER.get(rel):
        raise Fail(f"{rel} is still the untouched starter; write your code under the TODO comments")


def read_csv(rel, min_rows=1):
    """Rows as dicts with lower-case keys. An unnamed first column (row names
    written by write.csv) is kept as 'gene' when there is no gene column."""
    path = ROOT / rel
    if not path.is_file():
        raise Fail(f"{rel} not found; run your script and commit the file it writes")
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise Fail(f"{rel} is empty")
        names = [(h or "").strip().lower() for h in reader.fieldnames]
        rows = []
        for raw in reader:
            row = {}
            for key, value in raw.items():
                key = (key or "").strip().lower()
                if key == "":
                    if "gene" in names:
                        continue
                    key = "gene"
                row[key] = (value or "").strip()
            rows.append(row)
    if len(rows) < min_rows:
        raise Fail(f"{rel} has {len(rows)} rows; expected at least {min_rows}")
    return rows


def need_columns(rows, rel, cols):
    missing = [c for c in cols if c not in rows[0]]
    if missing:
        raise Fail(f"{rel} needs the column(s) {', '.join(missing)}; found {', '.join(rows[0])}")


def num(row, col, rel, allow_na=False):
    value = row.get(col, "")
    if allow_na and value.upper() in ("NA", "NAN", ""):
        return None
    try:
        return float(value)
    except ValueError:
        raise Fail(f"{rel}: column {col} should hold a number, found '{value}'")


def near(value, target, tol):
    return abs(value - target) <= tol


# ---------------------------------------------------------------- tasks

def task1():
    need_edited("scripts/counts_to_results.R")
    rel = "results/library_sizes.csv"
    rows = read_csv(rel, 1)
    need_columns(rows, rel, ["sample", "library_size"])
    sizes = {r["sample"]: num(r, "library_size", rel) for r in rows}
    if len(rows) != 4:
        raise Fail(f"expected 4 rows in {rel} (one per sample), found {len(rows)}")
    if set(sizes) != set(LIBRARY):
        raise Fail(f"the sample names in {rel} should match practice_samples.csv (ctrl_1 ...), not the alignment file names")
    for sample, size in LIBRARY.items():
        if sizes[sample] != size:
            raise Fail(f"the library size of {sample} does not match the checkpoint; keep only the four count columns and rename them in the order of the sample sheet")


def task2():
    need_edited("scripts/counts_to_results.R")
    rel = "results/deseq2_results.csv"
    rows = read_csv(rel, 1)
    need_columns(rows, rel, ["gene", "log2foldchange", "padj"])
    if len(rows) != 12:
        raise Fail(f"expected 12 genes in {rel}, found {len(rows)}; save the full results table, not only the significant genes")
    table = {r["gene"]: r for r in rows}
    sig = set()
    for gene, r in table.items():
        padj = num(r, "padj", rel, allow_na=True)
        if padj is not None and padj < 0.05:
            sig.add(gene)
    if sig != set(BUILT_IN):
        raise Fail(f"{len(sig)} genes have padj < 0.05 in {rel}; the checkpoint expects the four built-in changes. Check the design (~ condition) and that the sample sheet order matches the count columns")
    for gene, direction in BUILT_IN.items():
        lfc = num(table[gene], "log2foldchange", rel)
        if (lfc > 0) != (direction == "up"):
            raise Fail("the fold changes point the wrong way; is control the reference level?")
        if not near(lfc, LFC[gene], 0.15):
            raise Fail(f"the log2 fold change of {gene} is further from DESeq2's estimate than expected; use design = ~ condition on the raw counts")
    rel = "results/significant_genes.csv"
    rows = read_csv(rel, 1)
    need_columns(rows, rel, ["gene"])
    if {r["gene"] for r in rows} != set(BUILT_IN) or len(rows) != 4:
        raise Fail(f"{rel} should list exactly the genes with padj < 0.05, with NA values left out")


def task3():
    need_edited("scripts/counts_to_results.R")
    rel = "results/truth_comparison.csv"
    rows = read_csv(rel, 1)
    need_columns(rows, rel, ["gene", "truth", "deseq2", "cpm"])
    table = {r["gene"]: r for r in rows}
    if set(table) != set(TRUTH):
        raise Fail(f"{rel} should have one row for each built-in gene: {', '.join(TRUTH)}")
    for gene, truth in TRUTH.items():
        r = table[gene]
        if not near(num(r, "truth", rel), truth, 0.01):
            raise Fail(f"the truth column for {gene} should be the built-in change as a plain ratio")
        des = num(r, "deseq2", rel)
        if near(des, LFC[gene], 0.15) and not near(des, 2 ** LFC[gene], 0.3):
            raise Fail("the deseq2 column looks like log2 fold changes; convert them to plain ratios with 2^log2FoldChange")
        if des <= 0 or abs(math.log2(des) - LFC[gene]) > 0.15:
            raise Fail(f"the DESeq2 ratio for {gene} does not match your DESeq2 results; use 2^log2FoldChange")
        cpm = num(r, "cpm", rel)
        if cpm <= 0 or abs(cpm / CPM[gene] - 1) > 0.02:
            raise Fail(f"the CPM ratio for {gene} looks wrong; divide by each sample's library size, average each group, then divide treated by control")


def section(text, name):
    match = re.search(rf"^#+\s*{name}\b(.*?)(?=^#+\s|\Z)", text, re.S | re.M | re.I)
    return match.group(1) if match else ""


def task4():
    need_edited("scripts/counts_to_results.R")
    path = ROOT / "results" / "volcano.png"
    if not path.is_file():
        raise Fail("results/volcano.png not found; save the volcano plot with ggsave() and commit it")
    head = path.read_bytes()[:32]
    if head[:8] != b"\x89PNG\r\n\x1a\n":
        raise Fail("results/volcano.png is not a PNG file")
    width, height = struct.unpack(">II", head[16:24])
    if width < 600 or height < 600:
        raise Fail("results/volcano.png is very small; save it at a readable size and dpi")
    path = ROOT / "results" / "methods.md"
    if not path.is_file():
        raise Fail("results/methods.md not found; write the methods paragraph and limitations and commit the file")
    text = path.read_text(encoding="utf-8", errors="replace")
    low = text.lower()
    needs = [
        ("DESeq2 with its version number", r"deseq2\D{0,20}\d+\.\d+"),
        ("the R version", r"\br\b\D{0,20}\d+\.\d+"),
        ("the aligner (HISAT2)", r"hisat"),
        ("the counting tool (featureCounts)", r"featurecounts"),
        ("the design formula", r"~\s*condition"),
        ("the reference level", r"control"),
        ("the cut-off", r"0\.05"),
    ]
    missing = [name for name, pattern in needs if not re.search(pattern, low)]
    if missing:
        raise Fail("your methods do not yet give " + "; ".join(missing))
    lim = section(text, "Limitations").lower()
    if not lim.strip():
        raise Fail("results/methods.md needs a '## Limitations' section with two limitations")
    lim_missing = []
    if not re.search(r"(two|2)\s+(biological\s+)?(replicates?|samples?)|replicates", lim):
        lim_missing.append("how many samples each group has")
    if not re.search(r"simulat", lim):
        lim_missing.append("where the data came from")
    if lim_missing:
        raise Fail("your limitations do not yet mention " + " or ".join(lim_missing))


TASKS = [
    ("Task 1 · Counts", task1),
    ("Task 2 · Differential expression", task2),
    ("Task 3 · Compare with the truth", task3),
    ("Task 4 · A volcano plot and a methods paragraph", task4),
]


def main():
    failed = 0
    for name, func in TASKS:
        try:
            func()
            print(f"PASS  {name}")
        except Fail as err:
            failed += 1
            print(f"FAIL  {name}: {err}")
    print()
    if failed:
        print(f"{failed} of {len(TASKS)} checks failed. Read the message, fix that task, then push again.")
        sys.exit(1)
    print("All checks passed. Well done.")


if __name__ == "__main__":
    main()
