"""
PHB_gtdb — 工具函数模块
"""

import os
import sys
import gzip
import subprocess
import logging
import shutil
from pathlib import Path
from datetime import datetime
from typing import List, Dict, Tuple, Optional

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# ==============================================================================
# 日志配置
# ==============================================================================

def setup_logging(log_dir: Path, name: str = "phb_gtdb") -> logging.Logger:
    """配置双输出日志（文件 + 终端）。"""
    log_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    log_file = log_dir / f"{name}_{timestamp}.log"

    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)

    # 文件 handler
    fh = logging.FileHandler(log_file, encoding="utf-8")
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(funcName)s | %(message)s"
    ))

    # 终端 handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setLevel(logging.INFO)
    ch.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(message)s",
        datefmt="%H:%M:%S"
    ))

    logger.addHandler(fh)
    logger.addHandler(ch)
    logger.info(f"日志文件: {log_file}")
    return logger


# ==============================================================================
# FASTA 读写
# ==============================================================================

def read_fasta(path: Path) -> List[SeqRecord]:
    """读取 FASTA 文件（支持 .gz 压缩）。"""
    if not path.exists():
        raise FileNotFoundError(f"文件不存在: {path}")

    if path.suffix == ".gz":
        with gzip.open(path, "rt") as f:
            return list(SeqIO.parse(f, "fasta"))
    else:
        return list(SeqIO.parse(path, "fasta"))


def write_fasta(records: List[SeqRecord], path: Path, compress: bool = False):
    """写入 FASTA 文件。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    if compress:
        with gzip.open(path.with_suffix(".fasta.gz"), "wt") as f:
            SeqIO.write(records, f, "fasta")
    else:
        with open(path, "w") as f:
            SeqIO.write(records, f, "fasta")


def count_sequences(path: Path) -> int:
    """快速计数 FASTA 中的序列数。"""
    if path.suffix == ".gz":
        with gzip.open(path, "rt") as f:
            return sum(1 for line in f if line.startswith(">"))
    else:
        with open(path) as f:
            return sum(1 for line in f if line.startswith(">"))


# ==============================================================================
# 系统命令
# ==============================================================================

def run_cmd(cmd: str, desc: str = "", timeout: int = None,
            check: bool = True, logger: logging.Logger = None) -> subprocess.CompletedProcess:
    """运行 shell 命令，记录日志。

    Returns:
        subprocess.CompletedProcess 对象
    """
    if logger:
        logger.info(f"执行: {desc or cmd[:120]}")

    bash = shutil.which("bash") or "/bin/bash"
    result = subprocess.run(
        cmd, shell=True, capture_output=True, text=True,
        timeout=timeout, executable=bash
    )

    if result.returncode != 0 and check:
        err_msg = result.stderr[:500] if result.stderr else "no stderr"
        if logger:
            logger.error(f"命令失败 (code={result.returncode}): {err_msg}")
        raise RuntimeError(f"命令失败:\nCMD: {cmd[:200]}\n{err_msg}")

    return result


# ==============================================================================
# GTDB 分类学
# ==============================================================================

def parse_gtdb_taxonomy(tax_string: str) -> Dict[str, str]:
    """解析 GTDB 分类学字符串。

    GTDB 格式: d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;...

    Returns:
        {rank: taxon_name} 字典
    """
    rank_map = {"d__": "domain", "p__": "phylum", "c__": "class",
                 "o__": "order", "f__": "family", "g__": "genus", "s__": "species"}
    taxonomy = {}
    for part in tax_string.split(";"):
        part = part.strip()
        for prefix, rank in rank_map.items():
            if part.startswith(prefix):
                taxonomy[rank] = part[len(prefix):]
                break
    return taxonomy


def load_gtdb_taxonomy(tax_file: Path = None) -> pd.DataFrame:
    """加载 GTDB 分类学文件。

    Returns:
        DataFrame: columns = [accession, taxonomy_string]
    """
    if tax_file is None:
        from config import GTDB_TAXONOMY
        tax_file = GTDB_TAXONOMY

    df = pd.read_csv(tax_file, sep="\t", header=None,
                     names=["accession", "taxonomy"])
    return df


def load_gtdb_metadata(meta_file: Path = None) -> pd.DataFrame:
    """加载 GTDB 元数据文件。"""
    if meta_file is None:
        from config import GTDB_METADATA
        meta_file = GTDB_METADATA

    return pd.read_csv(meta_file, sep="\t", low_memory=False)


# ==============================================================================
# 序列工具
# ==============================================================================

def sanitize_id(seq_id: str) -> str:
    """清洗序列 ID。"""
    import re
    return re.sub(r"[^\w\-\.]", "_", seq_id)


def filter_by_length(records: List[SeqRecord],
                     min_len: int = 50,
                     max_len: int = 100000) -> List[SeqRecord]:
    """按序列长度筛选。"""
    before = len(records)
    filtered = [r for r in records if min_len <= len(r.seq) <= max_len]
    if before > len(filtered):
        logging.getLogger("phb_gtdb").info(
            f"长度筛选: {before} → {len(filtered)} (范围 {min_len}-{max_len} aa)"
        )
    return filtered
