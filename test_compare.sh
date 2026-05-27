#!/bin/bash
curl -s -X POST "http://localhost:8899/api/compare" \
  -H "Content-Type: application/json" \
  -d '{"pdb_ids":["1ema","1lz1"],"chain":"A"}' | \
  python3 -c "
import sys, json
d = json.load(sys.stdin)
print('Identity:', d.get('identity'), '%')
print('RMSD:', d.get('rmsd'), 'Å')
print('1ema SS: H', d['ss_comparison']['1ema']['H']['pct'], '% E', d['ss_comparison']['1ema']['E']['pct'], '% C', d['ss_comparison']['1ema']['C']['pct'], '%')
print('1lz1 SS: H', d['ss_comparison']['1lz1']['H']['pct'], '% E', d['ss_comparison']['1lz1']['E']['pct'], '% C', d['ss_comparison']['1lz1']['C']['pct'], '%')
print('Conserved regions:', len(d.get('conserved_regions', [])))
print('SASA pairs:', len(d.get('sasa_comparison', [])))
for r in d.get('conserved_regions', [])[:3]:
    print(f'  Region {r[\"start\"]}-{r[\"end\"]} ({r[\"length\"]}aa, {r[\"identity\"]}%)')
"
