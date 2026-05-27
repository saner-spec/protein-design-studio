#!/bin/bash
echo "=========================================="
echo "1. DSSP 二级结构导出"
echo "=========================================="
curl -s "http://localhost:8899/api/dssp/1ema" | python3 -m json.tool | head -30
echo ""
echo "=========================================="
echo "2. 实验方案生成 (SDM GLU222→ALA)"
echo "=========================================="
curl -s -X POST "http://localhost:8899/api/protocol/1ema" \
  -H "Content-Type: application/json" \
  -d '{"protocol_type":"sdm","mutation_desc":"GLU222→ALA"}' | python3 -m json.tool | head -30
echo ""
echo "=========================================="
echo "3. 校准集创建"
echo "=========================================="
curl -s -X POST "http://localhost:8899/api/calibrate/create" \
  -H "Content-Type: application/json" \
  -d '{"name":"gfp_test","entries":[
    {"position":222,"wt":"E","mut":"A","exp_ddg":-0.3,"esm_score":-3.72},
    {"position":222,"wt":"E","mut":"D","exp_ddg":-1.1,"esm_score":-2.15},
    {"position":222,"wt":"E","mut":"W","exp_ddg":-0.8,"esm_score":-5.02}
  ]}' | python3 -m json.tool
echo ""
echo "=========================================="
echo "4. 校准预测"
echo "=========================================="
curl -s -X POST "http://localhost:8899/api/calibrate/gfp_test/predict?esm_score=-4.0" | python3 -m json.tool
echo ""
echo "=========================================="
echo "All tests complete"
echo "=========================================="
