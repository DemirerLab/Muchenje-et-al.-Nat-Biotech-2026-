# fcs_pipeline

End-to-end flow cytometry analysis pipeline for *N. benthamiana* protoplast experiments. Loads FCS files, applies a sequential 4-gate hierarchy, computes statistics, and saves publication-quality figures.

Demirer Lab, Caltech Division of Chemistry and Chemical Engineering (CCE)

---

## Gate hierarchy

```
All events
    └── P1  (FSC-A vs SSC-A)              — intact protoplast scatter gate
            └── P2  (FSC-A vs FSC-Width)  — singlet gate (doublet exclusion)
                    └── P_MGOLD  (FSC-A vs 525Blue-H)  — transformation reporter / viability gate
                                └── MC  (FSC-A vs 610Yellow-H)  — target fluorescence readout gate
```

## Figures produced

- `gating_strategy.svg` — 4 samples × 4 gate steps grid with pseudocolor density scatter
- `final_fluorescence_gate.pdf` — transformation reporter and target fluorescence across samples

## Usage

1. Set `DATA_DIR` to the folder containing your `.fcs` files
2. Set `OUT_DIR` to where you want figures saved
3. Edit `GATING_SAMPLES` and `SUPP_SAMPLES` with your filenames and labels
4. Run: `python fcs_pipeline.py`

## Dependencies

```bash
pip install fcsparser numpy matplotlib scipy
```

Python ≥ 3.8

## License

MIT — see [LICENSE](LICENSE)
