#!/bin/bash
cd /mnt/e/AI_Agents/protein_designer
.venv/bin/python -c 'from backend.main import app; print("imports OK, routes:", len(app.routes))'
