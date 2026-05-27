# Changelog

All notable changes to Protein Design Studio.

## [0.4.0] — 2026-05-22

### Added
- ESM Validation Panel: ProTherm benchmark dataset validation (s669_mini, 14 verified mutations)
- Statistical metrics: Spearman ρ, Pearson r, directional accuracy, confusion matrix
- Residue type breakdown and ESM score binning
- Canvas scatter plot visualization (no third-party chart library dependency)
- Validation result caching for instant reload
- `scripts/download_benchmarks.py` — ProTherm/S669 dataset downloader with multi-source fallback
- `scripts/download_via_windows.py` — Windows proxy-aware downloader

### Changed
- Validator uses chain-aware caching (`pdb_id, chain, position` triple)
- PDB wild-type consistency check automatically skips mismatched entries

## [0.3.0] — 2026-05-21

### Added
- Knowledge Engine: AI-generated literature review with file-based cache
- Docker deployment: Dockerfile + docker-compose.yml for one-command startup
- PDB Search: search RCSB by name
- Design Report Generator: one-click report download
- 3D superposition visualization for multi-PDB comparison
- Collapsible protocol cards in UI
- Batch CSV export
- Ramachandran plot visualization
- AI combinatorial mutation suggestions

### Changed
- Frontend calibration panel with tab-based experimental data input

## [0.2.0] — 2026-05-21

### Added
- DSSP phi/psi angle calculation and export
- Experimental Protocol Generator (5 templates: point mutation, saturation, insertion, truncation, combinatorial)
- Private Calibration: linear regression + prediction intervals with persistence
- Multi-PDB Comparison: sequence alignment, RMSD, SASA, secondary structure, conserved regions

## [0.1.0] — 2026-05-21

### Added
- PDB parser with SEQRES sequence and three-layer residue mapping
- ESM-2 650M single-residue scoring (single forward pass for all 20 amino acids)
- 3Dmol.js 3D viewer with click-to-select residues
- AI Tutor: DeepSeek v4-pro streaming SSE chat with protein context injection
- Batch surface residue scan with CSS grid heatmap rendering
- 7-point startup diagnostics (Python/CUDA/VRAM/API/port/disk/ESM)
- PDB shortcut buttons (GFP, lysozyme, RNase, hemoglobin)
- Primer designer with PCR conditions and Tm calculation
- `download_esm.py` for pre-downloading ESM-2 model weights
