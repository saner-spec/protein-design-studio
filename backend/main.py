"""Protein Design Studio — FastAPI 主应用"""
import json
from pathlib import Path
from fastapi import FastAPI, HTTPException, Query
from fastapi.staticfiles import StaticFiles
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, StreamingResponse

from .config import HOST, PORT, DS_API_KEY, DS_API_BASE, DS_MODEL, DEFAULT_SASA_THRESHOLD, ESM_FULL_SCAN_MAX_LEN
from .diagnostics import full_diagnostics, format_diagnostics
from .pdb_parser import (
    parse_pdb, fetch_pdb, get_residue_info, is_surface_residue,
)
from .esm_scorer import score_single_position, score_batch_positions, load_model, check_vram
from .residue_map import STANDARD_AAS
from .ai_client import chat_stream, explain_residue
from .dssp_exporter import export_dssp, get_ss_summary
from .protocol_generator import generate_protocol
from .calibrator import calibration_mgr, CalibrationDataset
from .structure_comparator import compare_structures
from .knowledge_engine import generate_knowledge, get_cached_knowledge, generate_knowledge_brief
from .pdb_search import search_pdb
from .report_generator import generate_report
from .primer_designer import (
    design_mutagenesis_primers,
    design_sequencing_primer,
    assess_primer,
    calc_tm_basic,
    gc_content,
    aa_to_dna,
)
from .validator import (
    list_benchmarks,
    get_benchmark,
    run_validation,
    BenchmarkDataset,
    ValidationRunner,
)

from pydantic import BaseModel

app = FastAPI(title="Protein Design Studio", version="0.4.0")

# CORS — 允许本地开发
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

# 静态文件
FRONTEND = Path(__file__).resolve().parent.parent / "frontend"
FRONTEND.mkdir(exist_ok=True)
app.mount("/static", StaticFiles(directory=str(FRONTEND)), name="static")

# 运行时状态
_app_state = {"structures": {}, "diagnostics": None}


# ──────────────────────── API 路由 ────────────────────────

@app.get("/api/health")
async def health():
    """健康检查 + 环境诊断"""
    diag = full_diagnostics()
    _app_state["diagnostics"] = diag
    return {
        "status": "ok" if diag["all_ok"] else "degraded",
        "diagnostics": diag,
    }


@app.get("/api/structure/{pdb_id}")
async def get_structure(pdb_id: str):
    """加载PDB结构并返回摘要"""
    pdb_id = pdb_id.lower()

    # 检查缓存
    if pdb_id in _app_state["structures"]:
        s = _app_state["structures"][pdb_id]
        return {
            "pdb_id": s["pdb_id"],
            "seqres_sequence": s["seqres_sequence"],
            "residue_count": s["residue_count"],
            "seqres_length": s.get("seqres_length", len(s["seqres_sequence"])),
            "chains": s["chains"],
            "residue_list": [
                {
                    "resnum": r["resnum"],
                    "aa": r["aa"],
                    "resname": r["resname"],
                    "chain": r["chain"],
                    "rel_sasa": round(r["rel_sasa"], 3),
                    "ss": r["secondary_structure"],
                    "surface": is_surface_residue(r),
                }
                for r in s["residues"]
            ],
            "ligands": [
                {"resnum": l["resnum"], "resname": l["resname"]}
                for l in s["ligands"]
            ],
        }

    # 加载
    try:
        pdb_path = await fetch_pdb(pdb_id)
    except Exception as e:
        raise HTTPException(404, '找不到蛋白结构 "' + pdb_id.upper() + '"。请检查PDB ID是否正确。PDB ID通常为4位字母数字组合（如 1EMA、4HHB、1LZ1）。')

    try:
        structure = parse_pdb(pdb_path)
    except Exception as e:
        raise HTTPException(500, "蛋白结构解析失败，PDB文件格式可能已损坏，请尝试其他PDB ID。")

    _app_state["structures"][pdb_id] = structure

    return {
        "pdb_id": structure["pdb_id"],
        "seqres_sequence": structure["seqres_sequence"],
        "residue_count": structure["residue_count"],
        "seqres_length": structure.get("seqres_length", len(structure["seqres_sequence"])),
        "chains": structure["chains"],
        "residue_list": [
            {
                "resnum": r["resnum"],
                "aa": r["aa"],
                "resname": r["resname"],
                "chain": r["chain"],
                "rel_sasa": round(r["rel_sasa"], 3),
                "ss": r["secondary_structure"],
                "surface": is_surface_residue(r),
            }
            for r in structure["residues"]
        ],
        "ligands": [
            {"resnum": l["resnum"], "resname": l["resname"]}
            for l in structure["ligands"]
        ],
    }


@app.get("/api/structure/{pdb_id}/pdb")
async def get_pdb_file(pdb_id: str):
    """返回原始PDB文件内容（给3Dmol.js渲染）"""
    try:
        pdb_path = await fetch_pdb(pdb_id)
    except Exception as e:
        raise HTTPException(404, '找不到蛋白结构文件 "' + pdb_id.upper() + '"')
    return FileResponse(pdb_path, media_type="chemical/x-pdb")


@app.get("/api/residue/{pdb_id}/{chain}/{resnum}")
async def get_residue(pdb_id: str, chain: str, resnum: int):
    """获取单个残基的详细信息"""
    pdb_id = pdb_id.lower()
    if pdb_id not in _app_state["structures"]:
        # 自动加载
        await get_structure(pdb_id)

    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "结构未加载")

    residue = get_residue_info(structure, chain, resnum)
    if not residue:
        raise HTTPException(404, f"残基 {chain}/{resnum} 未找到")

    return {
        "resnum": residue["resnum"],
        "chain": residue["chain"],
        "resname": residue["resname"],
        "aa": residue["aa"],
        "sasa": round(residue["sasa"], 1),
        "rel_sasa": round(residue["rel_sasa"], 3),
        "secondary_structure": residue["secondary_structure"],
        "phi": round(residue.get("phi"), 1) if residue.get("phi") is not None else None,
        "psi": round(residue.get("psi"), 1) if residue.get("psi") is not None else None,
        "surface": is_surface_residue(residue),
        "neighbors": residue["neighbors"],
        "hbonds": residue["hbonds"][:10],
        "seqres_idx": residue.get("seqres_idx", -1),
    }


class ChatRequest(BaseModel):
    message: str
    pdb_id: str = ""


class ScanRequest(BaseModel):
    sasa_threshold: float = 0.25


@app.post("/api/score/{pdb_id}/{chain}/{resnum}")
async def score_residue(
    pdb_id: str,
    chain: str,
    resnum: int,
):
    """对残基位置的所有20种氨基酸打分"""
    pdb_id = pdb_id.lower()
    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "请先加载结构: GET /api/structure/{pdb_id}")

    residue = get_residue_info(structure, chain, resnum)
    if not residue:
        raise HTTPException(404, f"残基 {chain}/{resnum} 未找到")

    seq_idx = residue.get("seqres_idx", -1)
    if seq_idx < 0:
        raise HTTPException(400, "该残基不在SEQRES序列中，无法打分")

    sequence = structure["seqres_sequence"]
    if not sequence or seq_idx >= len(sequence):
        raise HTTPException(400, "序列索引错误")

    # 确保ESM已加载
    load_model()

    wt_aa = residue["aa"]
    scores = score_single_position(sequence, seq_idx)

    # 构建返回: 包含分数+化学性质
    aa_properties = {
        "A": "疏水小", "R": "碱性大", "N": "极性", "D": "酸性",
        "C": "含硫", "Q": "极性", "E": "酸性", "G": "最小",
        "H": "芳香碱性", "I": "疏水β支", "L": "疏水", "K": "碱性长",
        "M": "含硫疏水", "F": "芳香疏水", "P": "环状刚性",
        "S": "极性小", "T": "极性β支", "W": "芳香大", "Y": "芳香极性", "V": "疏水β支",
    }

    result = []
    for aa, score in sorted(scores.items(), key=lambda x: x[1]):
        result.append({
            "aa": aa,
            "score": round(score, 4),
            "property": aa_properties.get(aa, ""),
            # 颜色编码
            "color": "good" if score < -0.5 else "neutral" if abs(score) <= 0.5 else "bad",
        })

    return {
        "pdb_id": pdb_id,
        "chain": chain,
        "resnum": resnum,
        "wt_aa": wt_aa,
        "scores": result,
        "note": "ESM-2 650M 预测，整体相关性约0.4-0.6，仅供筛选参考",
    }


@app.post("/api/scan/{pdb_id}/{chain}")
async def scan_surface(pdb_id: str, chain: str, sasa_threshold: float = 0.25):
    """批量扫描表面残基"""
    pdb_id = pdb_id.lower()
    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "请先加载结构")

    # 收集表面残基
    surface_residues = [
        r for r in structure["residues"]
        if r["chain"] == chain and is_surface_residue(r, sasa_threshold)
    ]

    if not surface_residues:
        return {"pdb_id": pdb_id, "chain": chain, "surface_count": 0, "heatmap": []}

    # 检查长度限制
    if len(structure["seqres_sequence"]) > ESM_FULL_SCAN_MAX_LEN:
        raise HTTPException(400, f"序列过长({len(structure['seqres_sequence'])}aa)，不支持全扫描")

    sequence = structure["seqres_sequence"]
    positions = [r["seqres_idx"] for r in surface_residues if r.get("seqres_idx", -1) >= 0]

    if not positions:
        return {"pdb_id": pdb_id, "chain": chain, "surface_count": 0, "heatmap": []}

    # 确保ESM已加载
    load_model()

    # 批量打分
    all_scores = score_batch_positions(sequence, positions)

    # 构建热图数据
    heatmap = []
    for r in surface_residues:
        pos = r.get("seqres_idx", -1)
        if pos < 0 or pos not in all_scores:
            continue
        row = {
            "resnum": r["resnum"],
            "resname": r["resname"],
            "aa": r["aa"],
            "rel_sasa": round(r["rel_sasa"], 3),
            "scores": [
                {"aa": aa, "score": round(sc, 4)}
                for aa, sc in sorted(all_scores[pos].items(), key=lambda x: x[1])
            ],
        }
        heatmap.append(row)

    return {
        "pdb_id": pdb_id,
        "chain": chain,
        "surface_count": len(heatmap),
        "threshold": sasa_threshold,
        "heatmap": heatmap,
    }


@app.post("/api/chat")
async def chat(req: ChatRequest):
    """AI导师对话 (流式SSE)"""
    context_parts = []
    if req.pdb_id and req.pdb_id in _app_state.get("structures", {}):
        s = _app_state["structures"][req.pdb_id]
        context_parts.append(f"当前蛋白质: PDB {req.pdb_id}, {s['residue_count']}残基")

    context = "\n".join(context_parts) if context_parts else ""

    async def generate():
        async for token in chat_stream(req.message, context=context):
            yield f"data: {token}\n\n"
        yield "data: [DONE]\n\n"

    return StreamingResponse(generate(), media_type="text/event-stream; charset=utf-8")


@app.post("/api/explain/{pdb_id}/{chain}/{resnum}")
async def explain_residue_endpoint(pdb_id: str, chain: str, resnum: int):
    """让AI解释一个残基（非流式）"""
    pdb_id = pdb_id.lower()
    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "结构未加载")
    
    residue = get_residue_info(structure, chain, resnum)
    if not residue:
        raise HTTPException(404, "残基未找到")

    explanation = await explain_residue({
        "resnum": residue["resnum"],
        "resname": residue["resname"],
        "aa": residue["aa"],
        "chain": residue["chain"],
        "sasa": round(residue["sasa"], 1),
        "secondary_structure": residue["secondary_structure"],
        "phi": round(residue.get("phi"), 1) if residue.get("phi") is not None else None,
        "psi": round(residue.get("psi"), 1) if residue.get("psi") is not None else None,
        "surface": is_surface_residue(residue),
    }, pdb_id=pdb_id)

    return {"explanation": explanation or "⚠️ AI解释暂时不可用（API未配置或超时）"}


# ──────────────────────── DSSP API ──────────────────────────────

@app.get("/api/dssp/{pdb_id}")
async def get_dssp(pdb_id: str):
    """导出DSSP格式二级结构分配"""
    pdb_id = pdb_id.lower()
    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        await get_structure(pdb_id)
        structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "结构未加载")
    return {"dssp": export_dssp(structure), "summary": get_ss_summary(structure)}


@app.get("/api/ramachandran/{pdb_id}")
async def get_ramachandran(pdb_id: str):
    """获取所有残基的phi/psi值用于拉氏图"""
    pdb_id = pdb_id.lower()
    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        await get_structure(pdb_id)
        structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "结构未加载")

    points = []
    for r in structure["residues"]:
        phi = r.get("phi")
        psi = r.get("psi")
        if phi is None or psi is None:
            continue
        points.append({
            "resnum": r["resnum"],
            "aa": r["aa"],
            "phi": round(phi, 1),
            "psi": round(psi, 1),
            "ss": r.get("secondary_structure", "C"),
        })
    return {"pdb_id": pdb_id, "points": points, "count": len(points)}


# ──────────────────────── 实验方案 API ──────────────────────────

class ProtocolRequest(BaseModel):
    protocol_type: str = "full"  # sdm|expression|purification|characterization|full
    mutation_desc: str = ""
    primers: str = ""             # 预设计的引物序列(可选, 留空则提示用户使用primer_designer API)
    strain_info: str = "BL21(DE3)"
    culture_volume: str = "500"
    iptg_conc: float = 0.5
    induction_temp: int = 25
    induction_time: int = 16


@app.post("/api/protocol/{pdb_id}")
async def get_protocol(pdb_id: str, req: ProtocolRequest):
    """生成实验方案"""
    protocol = generate_protocol(
        protocol_type=req.protocol_type,
        pdb_id=pdb_id,
        mutation_desc=req.mutation_desc,
        primers=req.primers or "(请使用 primer_designer API 获取精确引物序列)",
        strain_info=req.strain_info,
        culture_volume=req.culture_volume,
        iptg_conc=req.iptg_conc,
        induction_temp=req.induction_temp,
        induction_time=req.induction_time,
    )
    return {"protocol": protocol, "format": "markdown"}


# ──────────────────────── 私有校准 API ──────────────────────────

class CalibrationAddRequest(BaseModel):
    name: str
    entries: list[dict] = []  # [{position, wt, mut, exp_ddg, esm_score}]


@app.post("/api/calibrate/create")
async def create_calibration(req: CalibrationAddRequest):
    """创建校准数据集并添加实验数据"""
    if calibration_mgr.get(req.name):
        raise HTTPException(409, f"校准集 '{req.name}' 已存在，请先删除或使用其他名称")
    ds = calibration_mgr.create(req.name)
    for entry in req.entries:
        ds.add(
            position=entry.get("position", 0),
            wt=entry.get("wt", ""),
            mut=entry.get("mut", ""),
            exp_ddg=entry.get("exp_ddg", 0.0),
            esm_score=entry.get("esm_score", 0.0),
        )
    result = ds.fit()
    calibration_mgr.save(req.name)
    return result


@app.post("/api/calibrate/{name}/add")
async def add_calibration_entry(name: str, position: int, wt: str, mut: str,
                                 exp_ddg: float, esm_score: float):
    """向已有校准集添加单个数据点"""
    ds = calibration_mgr.get(name)
    if not ds:
        raise HTTPException(404, f"校准集 '{name}' 不存在")
    ds.add(position=position, wt=wt, mut=mut, exp_ddg=exp_ddg, esm_score=esm_score)
    result = ds.fit()
    calibration_mgr.save(name)
    return result


@app.post("/api/calibrate/{name}/predict")
async def calibrate_prediction(name: str, esm_score: float):
    """用校准模型预测校准后的ΔΔG"""
    ds = calibration_mgr.get(name)
    if not ds:
        raise HTTPException(404, f"校准集 '{name}' 不存在")
    return ds.predict(esm_score)


@app.get("/api/calibrations")
async def list_calibrations():
    """列出所有校准集"""
    return {"calibrations": calibration_mgr.list_all()}


@app.get("/api/calibrate/{name}")
async def get_calibration(name: str):
    """获取校准集详情"""
    ds = calibration_mgr.get(name)
    if not ds:
        raise HTTPException(404, f"校准集 '{name}' 不存在")
    return ds.summary()


@app.delete("/api/calibrate/{name}")
async def delete_calibration(name: str):
    """删除校准集"""
    if calibration_mgr.delete(name):
        return {"ok": True}
    raise HTTPException(404, f"校准集 '{name}' 不存在")


# ──────────────────────── 多PDB对比 API ──────────────────────────

class CompareRequest(BaseModel):
    pdb_ids: list[str]
    chain: str = "A"


@app.post("/api/compare")
async def compare_pdbs(req: CompareRequest):
    """对比多个PDB结构"""
    if len(req.pdb_ids) < 2:
        raise HTTPException(400, "需要至少2个PDB ID")
    if len(req.pdb_ids) > 2:
        raise HTTPException(400, "当前版本仅支持2个PDB对比")

    structures = []
    for pid in req.pdb_ids:
        pid_lower = pid.lower()
        if pid_lower not in _app_state["structures"]:
            await get_structure(pid_lower)
        s = _app_state["structures"].get(pid_lower)
        if not s:
            raise HTTPException(404, f"PDB '{pid}' 加载失败")
        structures.append(s)

    result = compare_structures(structures, req.pdb_ids, req.chain)
    return result


# ──────────────────────── 知识引擎 API ──────────────────────────

@app.get("/api/knowledge/{pdb_id}")
async def get_knowledge(pdb_id: str, refresh: bool = False):
    """获取蛋白质知识页面 (缓存优先, refresh=true强制重新生成)"""
    pdb_id = pdb_id.lower()
    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        await get_structure(pdb_id)
        structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "结构未加载")

    result = await generate_knowledge(pdb_id, structure, refresh=refresh)
    return result


@app.get("/api/knowledge/{pdb_id}/brief")
async def get_knowledge_brief(pdb_id: str):
    """获取蛋白质简短摘要 (用于AI上下文注入)"""
    pdb_id = pdb_id.lower()
    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        await get_structure(pdb_id)
        structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "结构未加载")

    brief = await generate_knowledge_brief(pdb_id, structure)
    return {"brief": brief}


# ──────────────────────── PDB搜索 API ───────────────────────────

@app.get("/api/search")
async def search_pdb_endpoint(q: str, limit: int = 10):
    """按名称/关键词搜索PDB"""
    if not q or len(q) < 2:
        return {"results": []}
    results = await search_pdb(q, limit)
    return {"query": q, "results": results}


# ──────────────────────── CSV导出 API ──────────────────────────

@app.get("/api/export/csv/{pdb_id}/{chain}")
async def export_csv(pdb_id: str, chain: str, sasa_threshold: float = 0.25):
    """导出一个链所有表面残基的ESM打分CSV"""
    pdb_id = pdb_id.lower()
    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        await get_structure(pdb_id)
        structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "结构未加载")

    surface = [r for r in structure["residues"]
               if r["chain"] == chain and is_surface_residue(r, sasa_threshold)]
    if not surface:
        raise HTTPException(404, "无表面残基")

    if len(structure["seqres_sequence"]) > 500:
        raise HTTPException(400, "序列过长")

    sequence = structure["seqres_sequence"]
    positions = [r["seqres_idx"] for r in surface if r.get("seqres_idx", -1) >= 0]
    load_model()
    all_scores = score_batch_positions(sequence, positions)

    csv = "Position,WT,"
    csv += ",".join(sorted(STANDARD_AAS)) + "\n"

    for r in surface:
        pos = r.get("seqres_idx", -1)
        if pos < 0 or pos not in all_scores:
            continue
        csv += f"{r['resnum']},{r['aa']}"
        for aa in sorted(STANDARD_AAS):
            csv += f",{all_scores[pos].get(aa, ''):.4f}" if aa in all_scores[pos] else ","
        csv += "\n"

    return {"csv": csv, "filename": f"{pdb_id}_{chain}_esm_scores.csv"}


# ──────────────────────── 设计报告 API ──────────────────────────

class ReportRequest(BaseModel):
    include_protocol: bool = True
    include_knowledge: bool = True
    calibration_name: str = ""


@app.post("/api/report/{pdb_id}")
async def get_report(pdb_id: str, req: ReportRequest = ReportRequest()):
    """生成综合设计报告"""
    pdb_id = pdb_id.lower()
    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        await get_structure(pdb_id)
        structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "结构未加载")

    # 收集可选数据
    protocol = None
    calibration = None
    scores = None

    if req.calibration_name:
        ds = calibration_mgr.get(req.calibration_name)
        if ds:
            calibration = ds.summary()

    report = generate_report(
        structure=structure,
        scores=scores,
        calibration=calibration,
        protocol=protocol,
    )
    return {"report": report, "format": "markdown"}


# ──────────────────────── PCR引物设计 API ────────────────────────

class PrimerDesignRequest(BaseModel):
    target_aa: str
    primer_length: int = 30
    use_preferred_codons: bool = True


@app.post("/api/design-primers/{pdb_id}/{chain}/{resnum}")
async def design_primers(
    pdb_id: str,
    chain: str,
    resnum: int,
    req: PrimerDesignRequest,
):
    """设计定点突变PCR引物"""
    pdb_id = pdb_id.lower()
    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "请先加载结构")

    residue = get_residue_info(structure, chain, resnum)
    if not residue:
        raise HTTPException(404, f"残基 {chain}/{resnum} 未找到")

    wt_aa = residue["aa"]
    target_aa = req.target_aa.upper()
    if target_aa not in "ACDEFGHIKLMNPQRSTVWY":
        raise HTTPException(400, f"无效的目标氨基酸: {target_aa}")

    # 定位seqres_idx
    seq_idx = residue.get("seqres_idx", -1)
    if seq_idx < 0:
        raise HTTPException(400, "该残基不在SEQRES序列中")

    sequence = structure["seqres_sequence"]
    if seq_idx >= len(sequence):
        raise HTTPException(400, "序列索引错误")
    if sequence[seq_idx] != wt_aa:
        raise HTTPException(400, f"序列不匹配: 期望{wt_aa} 但序列中是{sequence[seq_idx]}")

    # 设计引物
    primers = design_mutagenesis_primers(
        wt_sequence=sequence,
        mutation_position=seq_idx,
        target_aa=target_aa,
        wt_aa=wt_aa,
        primer_length=req.primer_length,
        use_preferred_codons=req.use_preferred_codons,
        pdb_resnum=resnum,
    )

    return {
        "pdb_id": pdb_id,
        "chain": chain,
        "resnum": resnum,
        "resname": residue["resname"],
        "mutation": primers["mutation"],
        "wt_aa": wt_aa,
        "target_aa": target_aa,
        "wt_codon": primers["wt_codon"],
        "target_codon": primers["target_codon"],
        "codon_change": primers["codon_change"],
        "forward_primer": primers["forward_primer"],
        "reverse_primer": primers["reverse_primer"],
        "primer_tm_diff": primers["primer_tm_diff"],
        "recommended_annealing_temp": primers["recommended_annealing_temp"],
    }


@app.post("/api/design-workflow/{pdb_id}/{chain}/{resnum}/{target_aa}")
async def design_workflow(
    pdb_id: str,
    chain: str,
    resnum: int,
    target_aa: str,
    primer_length: int = 30,
):
    """完整设计工作流: ESM打分 + 引物设计 + 建议"""
    pdb_id = pdb_id.lower()
    structure = _app_state["structures"].get(pdb_id)
    if not structure:
        raise HTTPException(404, "请先加载结构")

    residue = get_residue_info(structure, chain, resnum)
    if not residue:
        raise HTTPException(404, f"残基 {chain}/{resnum} 未找到")

    wt_aa = residue["aa"]
    target_aa = target_aa.upper()
    if target_aa not in "ACDEFGHIKLMNPQRSTVWY":
        raise HTTPException(400, f"无效的目标氨基酸: {target_aa}")

    seq_idx = residue.get("seqres_idx", -1)
    if seq_idx < 0:
        raise HTTPException(400, "该残基不在SEQRES序列中")

    sequence = structure["seqres_sequence"]

    # 1. ESM打分
    load_model()
    esm_scores = score_single_position(sequence, seq_idx)
    esm_mutation_score = esm_scores.get(target_aa, None)

    # 2. 引物设计
    primers = design_mutagenesis_primers(
        wt_sequence=sequence,
        mutation_position=seq_idx,
        target_aa=target_aa,
        wt_aa=wt_aa,
        primer_length=primer_length,
        pdb_resnum=resnum,
    )

    # 3. 综合建议
    suggestions = []
    if esm_mutation_score is not None:
        if esm_mutation_score < -0.5:
            suggestions.append("✅ ESM预测该突变可能稳定蛋白质结构")
        elif esm_mutation_score < 0:
            suggestions.append("⚠️ ESM预测该突变影响较小，可能中性")
        else:
            suggestions.append("❌ ESM预测该突变可能去稳定化，建议谨慎")

    fwd = primers["forward_primer"]
    rev = primers["reverse_primer"]
    if fwd["quality"] == "good" and rev["quality"] == "good":
        suggestions.append("✅ 引物质量良好，可直接合成")
    elif fwd["quality"] == "fair" or rev["quality"] == "fair":
        suggestions.append("⚠️ 引物质量一般，建议优化长度或GC含量")
    else:
        suggestions.append("❌ 引物质量差，建议调整设计参数")

    if primers.get("primer_tm_diff") and primers["primer_tm_diff"] > 5:
        suggestions.append("⚠️ 正反向引物Tm差异较大(>5°C)，建议重新设计")

    return {
        "pdb_id": pdb_id,
        "chain": chain,
        "resnum": resnum,
        "resname": residue["resname"],
        "wt_aa": wt_aa,
        "target_aa": target_aa,
        "mutation": f"{wt_aa}{resnum}{target_aa}",
        "secondary_structure": residue.get("secondary_structure", "C"),
        "rel_sasa": round(residue.get("rel_sasa", 0), 3),
        "surface": is_surface_residue(residue),
        "esm_score": round(esm_mutation_score, 4) if esm_mutation_score is not None else None,
        "codon_change": primers["codon_change"],
        "forward_primer": primers["forward_primer"],
        "reverse_primer": primers["reverse_primer"],
        "primer_tm_diff": primers["primer_tm_diff"],
        "recommended_annealing_temp": primers["recommended_annealing_temp"],
        "suggestions": suggestions,
    }


# ──────────────────────── ESM 验证 API ───────────────────────────


class ValidationRequest(BaseModel):
    force_rerun: bool = False


@app.get("/api/validation/benchmarks")
async def get_benchmarks():
    """列出所有可用的基准验证数据集"""
    return {"benchmarks": list_benchmarks()}


@app.get("/api/validation/{benchmark_name}")
async def validate(benchmark_name: str, force_rerun: bool = False):
    """对基准数据集运行 ESM 验证，返回完整指标和散点图数据"""
    result = await run_validation(benchmark_name, force_rerun=force_rerun)
    if "error" in result:
        raise HTTPException(status_code=404, detail=result["error"])
    return result


@app.get("/api/validation/{benchmark_name}/scatter")
async def validation_scatter(benchmark_name: str):
    """仅返回散点图数据（用于轻量级刷新）"""
    ds = get_benchmark(benchmark_name)
    if ds is None:
        raise HTTPException(status_code=404, detail=f"数据集 '{benchmark_name}' 不存在")
    runner = ValidationRunner(ds)
    if not runner.load_cache():
        raise HTTPException(
            status_code=409,
            detail="ESM 缓存未生成，请先调用 GET /api/validation/{benchmark_name}"
        )
    runner._build_results()
    metrics = runner.compute_metrics()
    if "error" in metrics:
        raise HTTPException(status_code=500, detail=metrics["error"])
    return {"scatter": metrics.get("scatter", [])}


@app.post("/api/validation/{benchmark_name}/rerun")
async def rerun_validation(benchmark_name: str):
    """强制重新运行 ESM 验证（清除缓存）"""
    result = await run_validation(benchmark_name, force_rerun=True)
    if "error" in result:
        raise HTTPException(status_code=404, detail=result["error"])
    return result


# ──────────────────────── 前端入口 ────────────────────────

@app.get("/")
async def index():
    return FileResponse(FRONTEND / "index.html")


# ──────────────────────── 启动 ───────────────────────────

def start():
    """程序入口"""
    import uvicorn
    print(format_diagnostics(full_diagnostics()))
    uvicorn.run(app, host=HOST, port=PORT, log_level="info")


if __name__ == "__main__":
    start()
