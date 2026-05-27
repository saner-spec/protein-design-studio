#!/bin/bash
cd /mnt/e/AI_Agents/protein_designer
.venv/bin/python << 'PYEOF'
# Check for circular imports
import sys
sys.path.insert(0, ".")
from backend import main
print("All imports OK, no circular imports")
print("Routes:", len(main.app.routes))
# Check specific module imports
from backend import pdb_parser, esm_scorer, ai_client, residue_map
from backend import dssp_exporter, protocol_generator, calibrator
from backend import structure_comparator, knowledge_engine
from backend import pdb_search, report_generator, primer_designer
from backend import config, diagnostics
print("All modules load cleanly")
PYEOF
