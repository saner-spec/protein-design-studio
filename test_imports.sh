#!/bin/bash
cd /mnt/e/AI_Agents/protein_designer
.venv/bin/python -c 'from backend.dssp_exporter import export_dssp; print("dssp ok")'
.venv/bin/python -c 'from backend.protocol_generator import generate_protocol; print("protocol ok")'
.venv/bin/python -c 'from backend.calibrator import calibration_mgr; print("calibrator ok")'
.venv/bin/python -c 'from backend.primer_designer import design_mutagenesis_primers; print("primer ok")'
echo "---all imports passed---"
