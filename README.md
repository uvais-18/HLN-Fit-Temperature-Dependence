
# HLN Fit on Polycrystalline Sample

This Python script performs the Hikami-Larkin-Nagaoka (HLN) fit on magneto-conductance (MC) data of a polycrystalline sample at various temperatures.

It reads experimental CSV files of MR vs Magnetic Field for temperatures 2K to 10K, fits the data using the HLN formula, and extracts fitting parameters:

- α (alpha) — dimensionless parameter indicating WAL/WL contribution
- L_phit (Phase Coherence Length) — extracted in nanometers (nm)

Finally, it plots:

- HLN fits on experimental data
- Temperature dependence of α and L_phi
- Tells about the localization.

