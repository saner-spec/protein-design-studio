#!/usr/bin/env python3
"""Download ProTherm data via Windows curl.exe.
Strategy: submit ProThermDB download form, parse response for download link.
"""
import subprocess
import sys
import json
import re
import urllib.parse
from pathlib import Path

PROJECT = Path("/mnt/e/AI_Agents/protein_designer")
OUTFILE = PROJECT / "data" / "s669_downloaded.csv"
BENCHMARK_OUT = PROJECT / "data" / "benchmarks" / "protherm_full.json"
CURL_EXE = "/mnt/c/Windows/System32/curl.exe"

def curl(url, timeout=15, method="GET", data=None):
    """Run Windows curl.exe, return (http_code, body, effective_url)"""
    cmd = [CURL_EXE, "-sL", "-w", "EFFECTIVE:%{url_effective}\nHTTP:%{http_code}",
           "--max-time", str(timeout),
           "-H", "User-Agent: Mozilla/5.0"]
    if method == "POST":
        cmd.extend(["-X", "POST"])
    if data:
        cmd.extend(["-d", data])
    cmd.append(url)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout+5)
        body = result.stdout
        effective = ""
        code = "000"
        for line in body.split("\n"):
            if line.startswith("EFFECTIVE:"):
                effective = line[10:]
            elif line.startswith("HTTP:"):
                code = line[5:]
        # Remove our marker lines from body
        body = re.sub(r'^(EFFECTIVE:|HTTP:).*\n?', '', body, flags=re.MULTILINE)
        return code.strip(), body.strip(), effective.strip()
    except Exception as e:
        return "ERR", str(e), ""

def main():
    # Step 1: Get the download page to find the form
    print("Step 1: Accessing ProThermDB Downloads page...")
    code, body, _ = curl("https://web.iitm.ac.in/bioinfo2/prothermdb/Downloads.html")
    print(f"  Response: {code}, {len(body)} bytes")

    # Find form action
    form_action = re.search(r'action\s*=\s*"([^"]*)"', body)
    ncform = re.search(r'name="__ncforminfo"\s+value="([^"]*)"', body)
    print(f"  Form action: {form_action.group(1) if form_action else 'NOT FOUND'}")
    print(f"  ncforminfo: {ncform.group(1)[:40] if ncform else 'NOT FOUND'}")

    if form_action and ncform:
        action = form_action.group(1)
        # Build full URL
        if action.startswith(".."):
            action = "https://web.iitm.ac.in/bioinfo2/prothermdb/" + action.lstrip("./")
        elif action.startswith("/"):
            action = "https://web.iitm.ac.in" + action

        # Step 2: Submit form
        print(f"\nStep 2: Submitting form to {action}...")
        form_data = urllib.parse.urlencode({
            "name": "Researcher",
            "email": "protein.design@research.org",
            "role": "Scientist",
            "insti": "University",
            "message": "Academic research on protein stability prediction",
            "__ncforminfo": ncform.group(1),
        })
        code, body, effective = curl(action, timeout=30, method="POST", data=form_data)
        print(f"  Response: {code}, {len(body)} bytes")
        print(f"  Effective URL: {effective}")

        # Look for download links in response
        links = re.findall(r'href="([^"]*\.(?:csv|tsv|txt|zip|xls)[^"]*)"', body, re.I)
        print(f"  Download links: {links}")

        # If the response IS the data (CSV/TSV)
        if len(body) > 500 and ("pdb" in body[:500].lower() or "ddg" in body[:500].lower()):
            print(f"\n  *** Response appears to be the data itself! ***")
            OUTFILE.write_text(body, encoding="utf-8")
            print(f"  Saved {len(body.splitlines())} lines to {OUTFILE}")
            convert_to_json(body)
            return

        # Follow redirect to get the file
        if effective and effective != action:
            print(f"\nStep 3: Following redirect to fetch data...")
            code, body, _ = curl(effective, timeout=30)
            print(f"  Response: {code}, {len(body)} bytes")
            if len(body) > 500:
                OUTFILE.write_text(body, encoding="utf-8")
                print(f"  Saved to {OUTFILE}")
                convert_to_json(body)
                return

    # Step 3: Try direct access to known ProThermDB data paths
    print("\nStep 3: Trying direct data paths on ProThermDB...")
    direct_paths = [
        "https://web.iitm.ac.in/bioinfo2/prothermdb/data/prothermdb_all.txt",
        "https://web.iitm.ac.in/bioinfo2/prothermdb/data/download.txt",
        "https://web.iitm.ac.in/bioinfo2/prothermdb/data/ProThermDB.tsv",
    ]
    for url in direct_paths:
        code, body, _ = curl(url)
        print(f"  {code}: {url}")
        if code == "200" and len(body) > 500:
            OUTFILE.write_text(body, encoding="utf-8")
            convert_to_json(body)
            return

    # Step 4: GitHub fallback - try direct repo clones with more patterns
    print("\nStep 4: GitHub fallback with expanded patterns...")
    github_repos = [
        ("KULL-Centre/_Data", ["main", "master"], ["s669.csv", "data/s669.csv", "S669.csv"]),
        ("allydunham/thesis", ["main"], ["s669.csv", "data/s669.csv"]),
        ("allydunham/DDGun", ["master", "main"], ["data/s669.csv", "data/S669.csv"]),
        ("haskellinger/pub-data", ["main"], ["s669.csv"]),
    ]
    for repo, branches, paths in github_repos:
        for branch in branches:
            for path in paths:
                url = f"https://raw.githubusercontent.com/{repo}/refs/heads/{branch}/{path}"
                code, body, _ = curl(url)
                if code == "200" and len(body) > 100:
                    print(f"  FOUND: {url}")
                    OUTFILE.write_text(body, encoding="utf-8")
                    convert_to_json(body)
                    return

    print("\n" + "=" * 60)
    print("ALL AUTOMATIC SOURCES EXHAUSTED")
    print("=" * 60)
    print("Manual download required:")
    print("1. Open in browser: https://web.iitm.ac.in/bioinfo2/prothermdb/Downloads.html")
    print("2. Fill the form and download the TSV file")
    print(f"3. Save to: {PROJECT / 'data' / 'protherm_raw.tsv'}")
    print("4. Run: python scripts/download_benchmarks.py --source prothermdb")
    sys.exit(1)


AA3TO1 = {
    "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H",
    "ILE":"I","LYS":"K","LEU":"L","MET":"M","ASN":"N","PRO":"P","GLN":"Q",
    "ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y","MSE":"M",
}

def convert_to_json(csv_text):
    lines = csv_text.strip().split("\n")
    header = lines[0].strip().lower().replace('"', '').split(",")

    col = {}
    for i, h in enumerate(header):
        h = h.strip()
        if "pdb" in h: col["pdb"] = i
        elif "chain" in h: col["chain"] = i
        elif h in ("position", "pos", "resi", "resnum"): col["pos"] = i
        elif h in ("wt", "wild_type", "wild"): col["wt"] = i
        elif h in ("mut", "mutation", "mutant"): col["mut"] = i
        elif h in ("ddg", "exp_ddg", "ddg_exp"): col["ddg"] = i
        elif "ph" in h: col["ph"] = i
        elif h in ("temp", "temperature", "t"): col["temp"] = i

    print(f"  Columns detected: {list(col.keys())}")

    def norm_aa(aa):
        aa = aa.strip().upper().strip('"')
        if len(aa) == 3: return AA3TO1.get(aa, "X")
        if len(aa) == 1 and aa in "ACDEFGHIKLMNPQRSTVWY": return aa
        return "X"

    entries = []
    for line in lines[1:]:
        if not line.strip(): continue
        parts = [p.strip().strip('"') for p in line.split(",")]
        try:
            entry = {"chain": "A", "ph": 7.0, "temp": 25.0}
            if "pdb" not in col: continue
            pdb = parts[col["pdb"]].lower()
            if not pdb or len(pdb) != 4: continue
            entry["pdb_id"] = pdb
            if "chain" in col and col["chain"] < len(parts): entry["chain"] = parts[col["chain"]] or "A"
            entry["position"] = int(float(parts[col["pos"]])) if "pos" in col else 0
            entry["wt"] = norm_aa(parts[col["wt"]]) if "wt" in col else "X"
            entry["mut"] = norm_aa(parts[col["mut"]]) if "mut" in col else "X"
            entry["exp_ddg"] = float(parts[col["ddg"]]) if "ddg" in col else 0.0
            if "ph" in col and col["ph"] < len(parts):
                try: entry["ph"] = float(parts[col["ph"]])
                except: pass
            if "temp" in col and col["temp"] < len(parts):
                try: entry["temp"] = float(parts[col["temp"]])
                except: pass

            if entry["wt"] == "X" or entry["mut"] == "X" or entry["position"] <= 0: continue
            if entry["wt"] == entry["mut"]: continue
            if abs(entry["exp_ddg"]) > 15: continue
            entries.append(entry)
        except (ValueError, IndexError):
            continue

    seen = set()
    unique = []
    for e in entries:
        key = (e["pdb_id"], e["chain"], e["position"], e["mut"])
        if key not in seen:
            seen.add(key)
            unique.append(e)

    print(f"  Parsed {len(unique)} unique entries")
    data = {"name": "protherm_full", "description": f"ProTherm 完整数据集，{len(unique)} 条实验突变。", "entries": unique}
    BENCHMARK_OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(BENCHMARK_OUT, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    print(f"  Saved to {BENCHMARK_OUT}")

if __name__ == "__main__":
    main()
