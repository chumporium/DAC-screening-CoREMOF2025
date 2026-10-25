# Data and code: screening metal–organic frameworks for direct air capture with the full composition of humid air

This repository contains the data and scripts behind the article

> E. Munandar, T. Tokumasu, N. Nasruddin, *Screening metal–organic frameworks for direct air capture with the full composition of humid air*, Chemical Engineering Journal (submitted).

Please cite the article and this dataset when using the material.

## Contents

| Folder | Content |
|---|---|
| `data/` | Result tables (CSV) of all screening stages, GCMC simulations, working-capacity iterations, energy balance and levelised cost of capture; `Supplementary_Data_S1.xlsx` (the spreadsheet published with the article, English column names; it contains the data underlying all figures and tables) |
| `data/compositions/` | Gas-phase compositions (mole fractions, temperature, pressure) used for every desorption and adsorption GCMC run |
| `inputs/` | Example RASPA3 input sets (simulation, force field and component files) for one adsorption run (KENVOR, 253 K, 50% RH) and one desorption run (KENVOR, adsorption 298 K, desorption 373 K and 1 kPa, iteration 2) |
| `scripts/screening/` | Structure selection, Widom insertion input generation, harvesting and convergence classification of Henry coefficients, pore-limiting diameter, heat capacity, percolation analysis |
| `scripts/gcmc/` | Five-component GCMC input generation (`gen_gcmc.py`, `gen_pindai.py`, `bangun_kerja.py`), harvesting and convergence diagnostics; `bangun_kerja.py` and `gen_pindai.py` import `gen_gcmc.py` and `gen_widom.py` (`scripts/screening/`), which must be on the Python path |
| `scripts/working_capacity/` | Self-consistent iteration of the desorption gas composition (`iterasi_kerja.py`), the 1 kPa desorption surrogate (`surrogat_desorpsi.py`), the comparison of desorption treatments (`banding_desorpsi.py`) and the test of Eq. (8) (`uji_kriteria_eq8.py`) |
| `scripts/process_cost/` | Equilibrium energy balance (`model_energi.py`), climate grid (`peta_energi.py`), experiment-corrected sensitivity (`sensitivitas_koreksi.py`), coupling to the plant model for the levelised cost of capture (`dacsys_gcmc.py`), amine reference in the same plant model (`lcoc_amina.py`) |

## Workflow and correspondence with the article

1. **Structure set (Section 2.1).** `stage_semua.py` selects one representative structure per material (7,820). Structures are **not** included: they are available from CoRE MOF 2025, and the entries derived from the Cambridge Structural Database are subject to its licence. All materials are identified by their CoRE MOF 2025 names (CSD refcode plus suffix `_ASR`, `_FSR` or `_ion`).
2. **Stages 1–5 (Sections 2.4, 3.1; Tables S4–S5, Fig. S1).** `gen_co2.py`, `gen_gas.py`, `gen_air.py`, `gen_widom.py` generate Widom inputs; `panen.py` and `golongkan_konvergensi.py` harvest and classify them. Results: `hasil_screening_co2.csv`, `hasil_screening_n2.csv`, `hasil_screening_selektivitas.csv`, `hasil_screening_air.csv`, `label_konvergensi_co2.csv`; water Henry coefficients of further materials: `kh_h2o.csv`.
3. **Percolation labels (Section 2.4).** `peta_saja.py`, `perkolasi.py`, `uji_perkolasi.py`; results `hasil_perkolasi_semua.csv`, `hasil_screening_perkolasi.csv` (`E_halangan` is the percolation energy E_p).
4. **Five-component GCMC at 298 K (Sections 2.5, 3.1–3.3; Fig. 2–4, Tables S7–S11).** `gen_gcmc.py`, `cek_gcmcriil.py`; results `gcmcriil_final.csv` and labels; tests A–C in `ujidac_hasil.csv`, `uji_B_22.csv`, `uji_C_hasil.csv`. Henry estimate versus GCMC: `korelasi_kh_gcmc.py` -> `korelasi_KH_GCMC.csv`. Further entries of the Ag-framework family (Table S8) and the six materials of the 415 ppm check (Table S6): `panen_gcmc.py kenfam` / `panen_gcmc.py uji415` -> `gcmc_kenfam.csv`, `gcmc_uji415.csv`; water Henry coefficients of the latter: `kh_air415.csv`.
5. **Adsorption grid and working capacity (Sections 2.5, 2.7, 3.4–3.5; Table 3, Fig. 5).** `bangun_kerja.py` builds runs from `data/compositions/komp_*.csv`; `panen_kerja.py` harvests them into `kerja_hasil.csv`; `iterasi_kerja.py` computes `kerja_dq.csv` (every iteration) and `kerja_final.csv` (last iteration) and writes the composition file for the next iteration. Conditions without working capacity or without convergence are listed in `kerja_label_tanpa_kapasitas.csv`. In `kerja_hasil.csv`, the run `KENVIL_ion__ads258_rh050` is the result of an extended simulation (2.5 × 10^5 production cycles), which replaced the original run because its water loading did not meet σ ≤ 10%. The column `murni` in these files is the purity on a wet basis; the article reports the dry basis (Supplementary Data S1, `purity_dry_pct`).
6. **Surrogate (Section 2.7).** `surrogat_desorpsi.py validasi` gives `surrogat_validasi.csv`; independent GCMC checks in `uji_surrogat_prediksi.csv` and `uji_surrogat_hasil.csv`; the 5 kPa test in `surrogat_validasi_p05.csv` (`surrogat_desorpsi.py validasi p05`).
7. **Desorption atmosphere and reachability (Section 3.4; Table 2, Fig. 6).** `banding_desorpsi.py` compares the released-gas treatment (single-stage flash and differential release) with preset desorption states and with N2, O2 and Ar treated as inert: `banding_desorpsi.csv`. `uji_kriteria_eq8.py` tests Eq. (8) on all 360 conditions: `uji_kriteria_eq8.csv`.
8. **Energy balance and climate grid (Sections 2.8, 3.6–3.8; Fig. 7–8, Table 4, Table S17).** `model_energi.py` gives `energi_ekonomi.csv`, `energi_terbaik.csv`, `energi_sensitivitas.csv` (purity on a dry basis in `murni_pct`, wet basis in `murni_basah_pct`); `peta_energi.py` gives `peta_iklim.csv`; `sensitivitas_koreksi.py` gives `sensitivitas_koreksi.csv` (Table S17).
9. **Levelised cost of capture (Sections 2.9, 3.9; Fig. 9, Tables S3, S18 and S19).** `lcoc_amina.py` evaluates the amine reference with the same plant model, energy supply and 21 sites (`lcoc_amina.csv`; run it first, since it also provides the site climate statistics). `dacsys_gcmc.py` replaces the sorbent model of the plant model dacsys by the GCMC results and gives `lcoc_gcmc.csv`; all MOFs are assigned the uniform PrISMa sorbent price (30 EUR/kg) and are replaced every five years. For the system-level effect of the common approximations, `PERLAKUAN=inert` or `PERLAKUAN=umpan python peta_energi.py` builds the maps with N2, O2 and Ar treated as inert or with a preset desorption gas of feed composition (`peta_iklim_inert.csv`, `peta_iklim_umpan.csv`), and the same variable in `dacsys_gcmc.py` gives `lcoc_gcmc_inert.csv` and `lcoc_gcmc_umpan.csv`. The plant model itself is **not** part of this repository; set the environment variable `DACSYS_PATH` to the folder that contains it.
10. **Operating conditions (Sections 2.9, 3.9; Table S19).** `pareto_lcoc.py` evaluates all combinations of desorption temperature and drying target and marks the conditions that are non-dominated with respect to LCOC and CO2 purity: `pareto_lcoc.csv`.

Scripts are provided as used in this work. File names and some column names are in Indonesian; the glossary below gives the translations. Paths to external data (CoRE MOF 2025 files, the plant model) are set through environment variables or command-line arguments, as stated at the top of each script.

## Glossary of column names and labels

| Term | Meaning |
|---|---|
| `nama` | material name (CoRE MOF 2025) |
| `rh`, `RH` | relative humidity; `rh050` = 50% |
| `Ta`, `Td` | adsorption and desorption temperature (K) |
| `P` | desorption pressure: `p01` = 1 kPa, `p05` = 5 kPa, `p20` = 20.265 kPa |
| `it`, `iN` (in run names) | iteration of the desorption gas composition |
| `komp` | component (N2, O2, AR, H2O, CO2) |
| `q_molkg`, `sigma_pct` | loading (mol kg-1) and relative standard error (%), σ = Δ95/2.7764 |
| `dH_kJmol` | mixture adsorption enthalpy of the component (kJ mol-1) |
| `dq`, `sd` | CO2 working capacity (mol kg-1) and its standard error |
| `dN2` … `dCO2` | working capacity of each component |
| `lepas_pct` | fraction of adsorbed CO2 released (%) |
| `murni`, `murni_pct` | CO2 mole fraction of the released gas (%) |
| `beda`, `batas`, `stabil` | change of dq between successive iterations, convergence limit 2·sqrt(sd_n² + sd_n-1²), converged (True/False) |
| `KH_*`, `serapan_CO2` | Henry coefficient (mol kg-1 Pa-1); Henry-regime CO2 uptake at 40.5 Pa (mol kg-1) |
| `S_CO2_N2` | CO2/N2 selectivity |
| `status`: `LOLOS`, `GAGAL`/`DIKELUARKAN`, `TAK DAPAT DINILAI`, `TAK_KONVERGEN` | passed, failed/excluded, not assessable, not converged |
| `label`: `sah` | valid |
| `tanpa kapasitas kerja` | no working capacity |
| `tidak konvergen` | not converged |
| `berosilasi` | oscillating |
| `CO2 terusir air` | CO2 displaced by water |
| `air renik`, `gas renik` | trace water, trace gas (σ > 10% at loadings of about one molecule per cell) |
| `terjangkau`: `bebas`, `termal`, `tersegel` | percolation label: free, thermally activated, sealed |
| `E_halangan`, `U_min` | percolation barrier and minimum energy (kJ mol-1) |
| `skenario` | process scenario A (humid feed) or B (dried feed) |
| `rh_toleransi` | humidity tolerated by the sorbent (% RH) |
| `panas_sensibel`, `panas_desorpsi`, `pengeringan`, `listrik`, `energi_GJt` | sensible heat, desorption heat, drying energy, electricity, total (GJ t-1 CO2) |
| `energi_tanpa_kering` | total energy without drying (GJ t-1 CO2) |
| `ongkos_*` | cost terms of the lower-bound estimate (USD t-1 CO2) |
| `sumber` | source run or campaign |
| `catatan`, `keterangan` | note: `Widom didominasi satu kejadian` = Widom average dominated by a single insertion; `nilai 20.000 siklus salah; pakai nilai siklus panjang` = the 20,000-cycle value is discarded and the long-run value is used; `produksi 600k masih naik` = loading still increasing after 6 × 10^5 production cycles; `dua lari panjang tidak cocok` = two long runs disagree |

## Software

RASPA3 (versions 3.0.30 and 3.1.0), Python 3.13 with pandas 2.2, NumPy 2.2, SciPy 1.14 and Matplotlib 3.10. Force fields and simulation settings are listed in Table S1 of the article.

## Licence

Scripts (`scripts/`): MIT licence (`LICENSE`). Data (`data/`, `inputs/`): Creative Commons Attribution 4.0 International, CC BY 4.0 (`LICENSE-DATA`).

## Contact

N. Nasruddin, Research Group of Energy and Environment, Department of Mechanical Engineering, Universitas Indonesia, Depok, 16424, Indonesia (nasruddin@eng.ui.ac.id)
