# Bioinformatics Sequence Analysis 

This repository contains two Python programs for essential bioinformatics tasks:

1. **Transversion Ratio** – computes the transversion ratio among species descending from a shared ancestral DNA sequence.
2. **Transmembrane Region Predictor** – identifies transmembrane regions within protein sequences using hydrophobicity analysis.

Both programs use basic Python syntax and file-based input/output to maintain simplicity and clarity.

---

## 1. Transversion Ratio Calculator

### Purpose

Calculates the **transversion ratio** between DNA sequences of related species.
A transversion is a substitution where a **purine (A, G)** is replaced by a **pyrimidine (C, T)** or vice versa. The ratio helps measure evolutionary divergence between species sharing a common ancestor.

### Input

A plain text file (e.g., `CS Midterm Problem 2.txt`) containing aligned DNA sequences.
Each line represents a single sequence.

**Example input:**

```
ACGTACGT
AGGTACGA
ACCTTCGT
```

### Output

Displays the **transversion ratio** on the console or writes it to an output file.

**Formula:**

```
transversion_ratio = total_transversions / total_sites_compared
```

**Example output:**

```
Transversion Ratio: 0.8333
```

### Key Features

* Handles variable-length aligned sequences.
* Skips positions with only purines or only pyrimidines.
* Uses `try-except` blocks for file reading errors.
* Uses only basic built-in Python functions.

### Usage

```
python transversion_ratio.py
```

---

## 2. Transmembrane Region Predictor

### Purpose

Identifies **transmembrane regions** in protein sequences based on hydropathy analysis.
Uses a sliding window approach to detect hydrophobic segments likely to span the cell membrane.

### Input

FASTA or plain text file containing one or more protein sequences.

**Example input:**

```
>protein1
MVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF
```

### Output

Lists predicted transmembrane regions with their start and end positions.

**Example output:**

```
Protein: protein1
Transmembrane Region: 23–41
```

### Key Features

* Calculates average hydropathy per window.
* Marks regions exceeding a defined hydrophobicity threshold.
* Adjustable parameters:

  * Window size (default: 19)
  * Hydropathy threshold (default: 1.6)

### Usage

```
python transmembrane_region.py
```

---

## Project Structure

```
bioinformatics-toolkit/
│
├── transversion_ratio.py
├── transmembrane_region.py
├── example_inputs/
│   ├── dna_sequences.txt
│   └── protein_sequences.fasta
├── outputs/
└── README.md
```

---

## Requirements

* Python 3.8 or higher
* No external libraries required

---

## Running the Programs

1. Clone the repository:

   ```
   git clone https://github.com/yourusername/bioinformatics-toolkit.git
   cd bioinformatics-toolkit
   ```
2. Place your input files in the directory.
3. Run the programs:

   ```
   python transversion_ratio.py
   python transmembrane_region.py
   ```

---

