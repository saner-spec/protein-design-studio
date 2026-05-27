"""ESM打分正确性测试

CC审核要求: 对已知突变效应的蛋白做回归测试, 防止位置偏移等silent bug。
"""
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from backend.esm_scorer import score_single_position, score_batch_positions
from backend.pdb_parser import parse_pdb, get_residue_info
from backend.residue_map import STANDARD_AAS


def test_esm_output_sanity():
    """ESM输出合理性检查"""
    # 短序列: 10个A
    seq = "AAAAAAAAAA"
    scores = score_single_position(seq, 4)

    # 野生型分数为0
    assert "A" in scores
    assert scores["A"] == 0.0

    # 所有20种氨基酸都有分数
    assert len(scores) == 20

    # 分数应该是有限的
    for aa, s in scores.items():
        assert abs(s) < 100, f"{aa} score {s} out of range"


def test_glycine_preference():
    """甘氨酸应该偏好在loop区域, 在疏水核心不受欢迎"""
    # 纯疏水序列
    seq_hydrophobic = "LLLLLLLLLL"
    scores_G_in_hydrophobic = score_single_position(seq_hydrophobic, 4).get("G", 999)

    # 纯甘氨酸序列
    seq_gly = "GGGGGGGGGG"
    scores_L_in_gly = score_single_position(seq_gly, 4).get("L", 999)

    # G在疏水环境中应该得分较低 (更负 = 更不被认可)
    # L在Gly环境中也不应该受欢迎
    assert scores_G_in_hydrophobic < 0, "G should be disfavored in hydrophobic core"
    assert scores_L_in_gly < 0, "L should be disfavored in Gly-rich region"


def test_consensus_position():
    """保守位置: 重复的氨基酸应该强烈偏好野生型"""
    seq = "A" * 20
    scores = score_single_position(seq, 10)

    # 所有其他氨基酸的分数应 <= 0 (不如野生型)
    for aa, s in scores.items():
        if aa != "A":
            assert s <= 0, f"{aa} should not be preferred over wild-type A in poly-A, got {s}"


def test_position_index_bounds():
    """位置边界检查"""
    seq = "ACDEFGHIKL"

    # 合法位置
    scores = score_single_position(seq, 0)  # N端
    assert len(scores) == 20

    scores = score_single_position(seq, 9)  # C端
    assert len(scores) == 20

    # 非法位置
    with pytest.raises(ValueError):
        score_single_position(seq, 10)  # 超出范围

    with pytest.raises(ValueError):
        score_single_position(seq, -1)


def test_batch_scoring_consistency():
    """批量打分与单残基打分一致性"""
    seq = "MKFLILFNILV"
    pos = 3

    single = score_single_position(seq, pos)
    batch = score_batch_positions(seq, [pos])

    assert pos in batch
    for aa in STANDARD_AAS:
        assert abs(single[aa] - batch[pos][aa]) < 0.001, f"Batch/single mismatch for {aa}"


def test_pdb_index_offset():
    """PDB残基编号 vs ESM序列索引的一致性

    用1EMA的解析结果验证 seqres_idx 映射正确。
    """
    pdb_path = Path(__file__).resolve().parent.parent.parent / "static" / "pdbs" / "1EMA.pdb"
    if not pdb_path.exists():
        pytest.skip("1EMA.pdb not found")

    s = parse_pdb(pdb_path)
    r = get_residue_info(s, "A", 222)

    assert r is not None, "GLU222 not found"
    assert r["aa"] == "E"
    assert r["seqres_idx"] >= 0, "SEQRES index not assigned"
    assert r["seqres_idx"] < len(s["seqres_sequence"]), "SEQRES index out of range"

    # SEQRES序列中该位置应该是E
    assert s["seqres_sequence"][r["seqres_idx"]] == "E", (
        f"SEQRES[{r['seqres_idx']}] = {s['seqres_sequence'][r['seqres_idx']]}, expected E"
    )
