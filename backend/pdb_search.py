"""Protein Design Studio — PDB 搜索

通过 RCSB Search API 按蛋白名称/关键词搜索PDB条目。
"""
import json
import httpx

RCSB_SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"


async def search_pdb(query: str, limit: int = 10) -> list[dict]:
    """按名称/关键词搜索PDB结构

    返回: [{pdb_id, title, organism, resolution, residue_count, method}, ...]
    """
    # RCSB Search API: 全文搜索
    search_json = {
        "query": {
            "type": "terminal",
            "service": "full_text",
            "parameters": {
                "value": query
            }
        },
        "return_type": "entry",
        "request_options": {
            "paginate": {
                "start": 0,
                "rows": limit
            },
            "results_content_type": ["experimental"],
            "sort": [{"sort_by": "score", "direction": "desc"}],
        }
    }

    async with httpx.AsyncClient(timeout=15.0) as client:
        try:
            resp = await client.post(
                RCSB_SEARCH_URL,
                json=search_json,
                headers={"Content-Type": "application/json"},
            )
            if resp.status_code != 200:
                return []

            data = resp.json()
            results = []

            for entry in data.get("result_set", []):
                pdb_id = entry.get("identifier", "")
                title = entry.get("title", "")
                # 提取结构信息
                struct_info = _extract_entry_info(entry)

                results.append({
                    "pdb_id": pdb_id,
                    "title": title[:120],
                    "organism": struct_info.get("organism", ""),
                    "resolution": struct_info.get("resolution"),
                    "residue_count": struct_info.get("residue_count"),
                    "method": struct_info.get("method", ""),
                })

            return results

        except (httpx.TimeoutException, httpx.ConnectError):
            return []


def _extract_entry_info(entry: dict) -> dict:
    """从RCSB条目中提取结构信息"""
    info = {}

    # 来源生物
    try:
        for entity in entry.get("rcsb_entry_info", {}).get("nonpolymer_bound_components", []):
            pass
        src = entry.get("rcsb_entry_info", {}).get("source_organism", [])
        if src:
            info["organism"] = src[0].get("scientific_name", "")
    except Exception:
        info["organism"] = ""

    # 分辨率
    try:
        info["resolution"] = entry.get("rcsb_entry_info", {}).get("resolution_combined", [None])[0]
    except Exception:
        info["resolution"] = None

    # 残基数
    try:
        info["residue_count"] = entry.get("rcsb_entry_info", {}).get("deposited_polymer_entity_instance_count")
    except Exception:
        info["residue_count"] = None

    # 实验方法
    try:
        methods = entry.get("rcsb_entry_info", {}).get("experimental_method", [])
        info["method"] = ", ".join(methods) if methods else ""
    except Exception:
        info["method"] = ""

    return info


async def quick_search(query: str) -> list[dict]:
    """快速搜索 — 只返回PDB ID和标题，用于前端自动补全"""
    results = await search_pdb(query, limit=8)
    return [
        {"pdb_id": r["pdb_id"], "title": r["title"][:80]}
        for r in results
    ]
