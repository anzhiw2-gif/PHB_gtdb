"""
PHB_gtdb — 全局配置文件
所有路径、参数集中管理，便于复现分析
"""

import os
from pathlib import Path

# ==============================================================================
# 项目路径
# ==============================================================================
PROJECT_ROOT = Path("/home/data/haoyu/PHB_gtdb")
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
PROCESSED_DIR = DATA_DIR / "processed"
EXTERNAL_DIR = DATA_DIR / "external"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = RESULTS_DIR / "figures"
TABLES_DIR = RESULTS_DIR / "tables"
LOGS_DIR = RESULTS_DIR / "logs"
SCRIPTS_DIR = PROJECT_ROOT / "scripts"

# GTDB 数据库路径 (R232, 服务器 10.16.1.141)
GTDB_ROOT = Path("/home/data/haoyu/GTDB")
GTDB_GENOMES = GTDB_ROOT / "gtdb_genomes_reps_r232" / "database"
GTDB_METADATA = GTDB_ROOT / "metadata" / "bac120_metadata_r232.tsv"
GTDB_TAXONOMY = GTDB_ROOT / "taxonomy" / "bac120_taxonomy_r232.tsv"
GTDB_TAXONOMY_AR53 = GTDB_ROOT / "taxonomy" / "ar53_taxonomy_r232.tsv"
GTDB_TREE = GTDB_ROOT / "GTDB_tree" / "bac120_r232.tree"

# Conda 环境路径 (无需 sudo)
CONDA_ENV = Path("/home/data/haoyu/miniconda3/envs/phb_gtdb")
CONDA_BIN = Path("/home/data/haoyu/miniconda3/bin")

# ==============================================================================
# PHB 相关基因配置
# ==============================================================================

# PHB 降解与合成相关基因
PHB_GENES = {
    "phaZ": {
        "name": "PHB depolymerase",
        "pfam": "PF14556",
        "description": "PHA depolymerase — 分解 PHB 塑料",
    },
    "phaC": {
        "name": "PHA synthase (Class I)",
        "pfam": "PF07167",
        "description": "PHA synthase N-terminal domain — 合成 PHB",
    },
    "phaC_class2": {
        "name": "PHA synthase (Class II)",
        "pfam": "PF06863",
        "description": "PHA synthase — 合成中链 PHA",
    },
    "phaA": {
        "name": "Acetyl-CoA acetyltransferase",
        "kegg": "K00626",
        "description": "Beta-ketothiolase — PHB 合成第一步",
    },
    "phaB": {
        "name": "Acetoacetyl-CoA reductase",
        "kegg": "K00023",
        "description": "Acetoacetyl-CoA reductase — PHB 合成第二步",
    },
    "phaP": {
        "name": "Phasin",
        "pfam": "PF09650",
        "description": "PHA granule-associated protein",
    },
    "phaR": {
        "name": "PHA synthesis repressor",
        "pfam": "PF08695",
        "description": "PHA synthesis regulatory protein",
    },
}

# 额外通用水解酶域 (可能与PHB降解有关)
ADDITIONAL_PFAMS = [
    "PF00561",  # Alpha/beta hydrolase fold (PhaZ属于此家族)
    "PF07859",  # Alpha/beta hydrolase family
    "PF12697",  # Alpha/beta hydrolase family
]

# ==============================================================================
# HMM 搜索参数
# ==============================================================================
HMM_EVALUE = 1e-5
HMM_DOM_EVALUE = 1e-3

# ==============================================================================
# BLAST 参数
# ==============================================================================
BLAST_EVALUE = 1e-10
BLAST_MAX_TARGET_SEQS = 500

# ==============================================================================
# 序列比对参数
# ==============================================================================
ALIGNMENT_TOOL = "mafft"  # mafft / muscle
ALIGNMENT_THREADS = 30

# ==============================================================================
# 系统发育参数
# ==============================================================================
PHYLO_BOOTSTRAP = 1000
PHYLO_MODEL = "MFP"  # ModelFinder Plus
PHYLO_THREADS = 30

# ==============================================================================
# 计算资源
# ==============================================================================
MAX_THREADS = 30
MEMORY_GB = 100

# ==============================================================================
# 可视化参数
# ==============================================================================
FIGURE_DPI = 300
FIGURE_FORMAT = "pdf"

# ==============================================================================
# 确保所有必要目录存在
# ==============================================================================
def ensure_project_dirs():
    """创建所有必要的项目目录。"""
    for d in [
        RAW_DIR, PROCESSED_DIR, EXTERNAL_DIR,
        FIGURES_DIR, TABLES_DIR, LOGS_DIR,
    ]:
        d.mkdir(parents=True, exist_ok=True)
