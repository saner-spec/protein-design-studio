# Protein Design Studio

AI-guided protein mutation design tool. Load a PDB → 3D browse → click a residue → ESM-2 scores all 20 amino acids → AI explains → export experimental protocols.

![Version](https://img.shields.io/badge/version-0.4.0-blue)
![License](https://img.shields.io/badge/license-MIT-green)
![Python](https://img.shields.io/badge/python-3.9%2B-blue)
![PyTorch](https://img.shields.io/badge/pytorch-2.5.1%2Bcu121-orange)

## Features

### Core
- **PDB Parser** — Full SEQRES sequence with three-layer residue mapping (PDB index → ESM index → residue type), multi-chain support, SASA, hydrogen bonds, phi/psi angles
- **ESM-2 650M Scoring** — Single forward pass scores all 20 amino acids at clicked residue (CUDA-accelerated)
- **3D Viewer** — 3Dmol.js with click-to-select residues, highlighting, and detail panel
- **AI Tutor** — DeepSeek v4-pro streaming SSE chat with protein context injection

### Advanced
- **Batch Scan + Heatmap** — Surface residue batch scoring rendered as CSS grid heatmap
- **ESM Validation Panel** — ProTherm benchmark dataset validation with Spearman/Pearson/directional accuracy/scatter plot
- **DSSP + Ramachandran** — phi/psi angle calculation and Ramachandran plot visualization
- **Multi-PDB Comparison** — Sequence alignment, RMSD, SASA, secondary structure, conserved regions
- **Experimental Protocol Generator** — 5 templates (point mutation, saturation, insertion, truncation, combinatorial)
- **Primer Designer** — PCR conditions + Tm calculation
- **Private Calibration** — Linear regression + prediction intervals with user-submitted experimental data
- **Knowledge Engine** — AI-generated literature review with caching
- **PDB Search** — Search RCSB by name
- **Report Generator** — One-click design report download

### Deployment
- **Docker** — One-command startup via docker-compose
- **Diagnostics** — 7-point startup check (Python/CUDA/VRAM/API/port/disk/ESM)
- **tmux Persistence** — Long-running server via tmux script

## Quick Start

### Docker (Recommended)

```bash
docker compose up -d
# Open http://localhost:8899
```

No Python setup needed. First run downloads ESM-2 model (~2.5GB, cached in Docker volume).

### Manual (WSL2 / Linux)

```bash
# 1. Clone and install dependencies
git clone https://github.com/YOUR_USER/protein-design-studio.git
cd protein-design-studio
pip install -r requirements.txt

# 2. Install PyTorch with CUDA support
pip install torch --index-url https://download.pytorch.org/whl/cu121
pip install fair-esm

# 3. (Optional) Pre-download ESM-2 650M model (~2.5GB)
python download_esm.py

# 4. Set API key for AI tutor
export DEEPSEEK_API_KEY="sk-..."

# 5. Start
python run.py
# Open http://localhost:8899
```

> **Note:** Windows native is not supported for fair-esm. Use WSL2 or Docker.
> If you encounter WSL2 localhost forwarding issues, use the WSL2 VM IP directly: `curl http://$(hostname -I | awk '{print $1}'):8899/api/health`

## API Endpoints

| Method | Path | Description |
|--------|------|-------------|
| GET | `/api/health` | Health check + 7-point diagnostics |
| GET | `/api/structure/{pdb_id}` | Load PDB + return residue list |
| POST | `/api/score/{pdb_id}/{chain}/{pos}` | ESM-2 score all 20 amino acids |
| POST | `/api/scan/{pdb_id}/{chain}` | Batch scan surface residues |
| POST | `/api/chat` | AI tutor streaming chat |
| POST | `/api/explain/{pdb_id}/{chain}/{pos}` | AI explain a residue |
| GET | `/api/dssp/{pdb_id}` | DSSP phi/psi data |
| GET | `/api/ramachandran/{pdb_id}` | Ramachandran plot data |
| POST | `/api/protocol/{pdb_id}` | Generate experimental protocol |
| POST | `/api/calibrate/create` | Create private calibration set |
| POST | `/api/calibrate/{name}/add` | Add calibration data point |
| POST | `/api/calibrate/{name}/predict` | Get calibrated predictions |
| POST | `/api/compare` | Multi-PDB comparison |
| GET | `/api/knowledge/{pdb_id}` | Literature knowledge review |
| GET | `/api/search` | PDB search |
| GET | `/api/export/csv/{pdb_id}/{chain}` | Batch CSV export |
| POST | `/api/report/{pdb_id}` | Design report download |
| POST | `/api/design-primers/...` | Primer design |
| POST | `/api/design-workflow/...` | Complete design workflow |
| GET | `/api/validation/benchmarks` | List available benchmarks |
| GET | `/api/validation/{name}` | Run ESM validation |
| GET | `/api/validation/{name}/scatter` | Scatter plot data (lightweight) |
| POST | `/api/validation/{name}/rerun` | Force re-run validation |

## ESM Validation Panel

The validation panel answers the meta-question: *"How trustworthy are these ESM scores?"*

Built on a curated subset of ProTherm experimental data (14 single-point mutations from T4 lysozyme, Barnase, and Staph Nuclease, verified for PDB wild-type consistency). Statistical metrics include:

- **Spearman ρ** — rank correlation
- **Pearson r** — linear correlation  
- **Directional Accuracy** — fraction of correct sign predictions
- **Confusion Matrix** — stabilizing/destabilizing classification
- **Residue Type Breakdown** — accuracy by hydrophobic/polar/charged/special

> **Known limitation:** ESM-2 shows systematic bias toward predicting all mutations as stabilizing. Directional accuracy on s669_mini is ~14%. Use validation results to stay honest about ESM's limitations, not to blindly trust scores.

For larger validation sets, run `scripts/download_benchmarks.py` to fetch the full ProTherm dataset (~30,000 entries, requires manual download from ProThermDB).

## Tech Stack

```
Backend:  Python 3.12 + FastAPI + fair-esm + PyTorch 2.5.1+cu121
Frontend: Vanilla HTML/CSS/JS + 3Dmol.js (CDN)
AI:       DeepSeek v4-pro (api.deepseek.com/v1)
Model:    ESM-2 650M (~2.5GB VRAM)
```

## Requirements

- Python 3.9+
- NVIDIA GPU with 6GB+ VRAM (or CPU, significantly slower)
- PyTorch 2.x with CUDA support
- fair-esm 2.x

## Project Structure

```
protein_designer/
├── run.py                     # One-click start (diagnostics → ESM → uvicorn)
├── download_esm.py            # Pre-download ESM-2 650M (~2.5GB)
├── Dockerfile                 # Docker build
├── docker-compose.yml         # One-command Docker deployment
├── start.bat                  # Windows batch launcher
├── start_tmux.sh              # tmux persistent server
├── requirements.txt
│
├── backend/
│   ├── main.py                # FastAPI app (~31KB)
│   ├── config.py              # Global configuration
│   ├── diagnostics.py         # 7-point startup diagnostics
│   ├── pdb_parser.py          # PDB parser (SEQRES + SASA + H-bonds + phi/psi)
│   ├── esm_scorer.py          # ESM-2 scoring pipeline
│   ├── residue_map.py         # Non-standard → standard amino acid mapping
│   ├── ai_client.py           # DeepSeek streaming SSE client
│   ├── validator.py           # ESM validation engine
│   ├── calibrator.py          # Private calibration (linear regression)
│   ├── dssp_exporter.py       # DSSP phi/psi + secondary structure
│   ├── knowledge_engine.py    # AI literature review + cache
│   ├── pdb_search.py          # RCSB PDB search
│   ├── primer_designer.py     # Primer design + PCR conditions
│   ├── protocol_generator.py  # Experimental protocol (5 templates)
│   ├── report_generator.py    # One-click report download
│   ├── structure_comparator.py # Multi-PDB comparison
│   └── tests/
│       └── test_esm_scorer.py
│
├── frontend/
│   └── index.html             # SPA (~85KB, vanilla HTML/CSS/JS + 3Dmol.js)
│
├── data/
│   └── benchmarks/
│       └── s669_mini.json     # Validation benchmark (14 verified mutations)
│
└── scripts/
    ├── download_benchmarks.py  # ProTherm/S669 dataset downloader
    └── download_via_windows.py # Windows proxy downloader
```

## Known Limitations

- **CUDA only on WSL2/Linux** — fair-esm does not support Windows natively
- **ESM-2 directional accuracy ~14%** — model systematically predicts all mutations as stabilizing; scores are useful for *ranking* but not for predicting sign of ΔΔG
- **Single-user, single-process** — no concurrency protection
- **Validation dataset is small** (14 entries) — constrained by PDB wild-type residue consistency with ProTherm
- **DSSP may label some residues as coil** if backbone atoms are missing

## Environment Variables

| Variable | Required | Description |
|----------|----------|-------------|
| `DEEPSEEK_API_KEY` | Yes | DeepSeek API key for AI tutor |
| `PORT` | No | Server port (default: 8899) |
| `ESM_MODEL_DIR` | No | Custom ESM model directory |

## Development History

- **2026-05-20** — Concept discussion
- **2026-05-21** — v0.1.0: Backend skeleton, PDB parsing, ESM pipeline, frontend panel
- **2026-05-21** — v0.2.0: DSSP, protocol generator, calibration, multi-PDB comparison
- **2026-05-21** — v0.3.0: Knowledge engine, Docker, PDB search, report download, Ramachandran
- **2026-05-22** — v0.4.0: ESM validation panel with ProTherm benchmarks

## License

MIT — see [LICENSE](LICENSE) for details.
