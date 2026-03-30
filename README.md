# Boltz Hinge Protein Prediction Pipeline

Pipeline to evaluate how [Boltz](https://github.com/jwohlwend/boltz) models hinge proteins in their apo (ligand-free) and holo (ligand-bound) conformations. For each target, we run two Boltz predictions (apo and holo), align against true crystal structures, compute RMSD, generate contact maps, and run DynDom domain motion analysis.

## Target Proteins

| Protein | Apo PDB | Holo PDB | Ligand |
|---|---|---|---|
| Calmodulin | 1EXR | 1PRW | Calcium ions (CA) |
| DNA polymerase β (rat) | 1BPD | 2BPG | DNA template |
| Ribose binding protein | 1URP | 2DRI | Ribose (RIB) |
| Uracil DNA glycosylase | 1AKZ | 1SSP | DNA/uracil |
| Glutamate binding protein | 1WDN | 1GGG | Glutamate (GLU) |
| Maltose binding protein | 1N3X | 1NL5 | Maltose (MAL) |
| IgG | 8VEV | 8TCA | (large, run separately) |

## Requirements

**Python packages:**

```bash
pip install -r requirements.txt
```

**External tools:**

- **PyMOL** (Incentive or open-source) — used headlessly via `pymol -cq`
- **Boltz** — `pip install boltz`
- **DynDom** — https://dyndom.cmp.uea.ac.uk/dyndom/ (web server or local install)

## Pipeline Steps

| Step | Script | Description |
|---|---|---|
| 1 | `scripts/01_fetch_sequences.py` | Fetch FASTA sequences from RCSB |
| 2 | `scripts/02_check_sequence_match.py` | Align apo/holo sequences, standardize on apo |
| 3 | `scripts/03_make_yamls.py` | Generate Boltz input YAMLs (apo + holo) |
| 4 | `scripts/04_run_boltz.sh` | Run Boltz predictions |
| 5 | `scripts/05_fetch_true_structures.py` | Download true crystal structures from RCSB |
| 6 | `scripts/06_align_and_rmsd.py` | Structural alignment and RMSD (PyMOL) |
| 7 | `scripts/07_contact_maps.py` | Cα contact maps and overlap metrics |
| 8 | `scripts/08_render_images.py` | PyMOL image rendering (pLDDT, ligand views) |
| 9 | `scripts/09_run_dyndom.sh` | DynDom domain motion analysis |

## Quick Start

Run the full pipeline:

```bash
bash run_all.sh
```

Or run individual steps:

```bash
python scripts/01_fetch_sequences.py
python scripts/02_check_sequence_match.py
python scripts/03_make_yamls.py
bash scripts/04_run_boltz.sh
python scripts/05_fetch_true_structures.py
pymol -cq scripts/06_align_and_rmsd.py
python scripts/07_contact_maps.py
pymol -cq scripts/08_render_images.py
bash scripts/09_run_dyndom.sh
```

## Output Structure

```
results/
├── metrics.csv                          # RMSD + contact map metrics
├── sequence_match_report.csv            # Apo/holo sequence identity
├── images/                              # PyMOL rendered PNGs
├── contact_maps/                        # Contact map PNGs
├── aligned_pdbs/                        # Post-alignment PDB files
└── dyndom/                              # DynDom output per pair
```

## Notes

- **IgG** (8VEV/8TCA) is flagged `skip_by_default: true` in `config.yaml` due to its large size. Set to `false` when ready to run it with sufficient GPU memory.
- **Boltz output format** must be PDB (`--output_format pdb`) for DynDom compatibility.
- **Sequence mismatch** between apo/holo structures is handled by step 2, which standardizes on the apo sequence.
- **pLDDT** is stored in the B-factor column of Boltz output; regions with pLDDT < 70 should be interpreted cautiously.
