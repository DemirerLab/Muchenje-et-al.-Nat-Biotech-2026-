# junction_seq_table

Analyses Nanopore PCR amplicon reads spanning the R2 retrotransposon 3′ junction (3′ UTR payload → 25S rDNA). Extracts and ranks junction sequences and renders a color-coded nucleotide table.

Demirer Lab, Caltech Division of Chemistry and Chemical Engineering (CCE)

---

## What it does

1. Orients reads using forward and reverse primers (≤2 mismatches)
2. Locates the junction by anchoring to a 20 bp reference sequence at the canonical nick site
3. Extracts an 11 bp window on each side of the junction
4. Groups and ranks identical junction sequences by frequency
5. Renders a color-coded nucleotide table (PNG + SVG) showing the top 12 sequences alongside the reference

## Outputs

- `{LABEL}_junction_table.png` — raster figure (200 dpi)
- `{LABEL}_junction_table.svg` — vector figure (editable in Illustrator)

## Usage

Edit the `Config` block at the top of the script with your FASTQ path, primer sequences, and reference template, then run:

```bash
python junction_seq_table.py
```

## Dependencies

```bash
pip install matplotlib
```

Python ≥ 3.8

## License

MIT — see [LICENSE](LICENSE)
