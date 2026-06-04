#!/bin/bash
# =============================================================================
# PHB_gtdb — 主分析流程
#
# 完整运行 PHB 降解相关基因的系统发育与功能分析
# 用法: bash main_pipeline.sh [--step N] [--threads 30]
# =============================================================================

set -euo pipefail

# === 默认参数 ===
THREADS=30
START_STEP=1
END_STEP=6
SCRIPTS_DIR="$(cd "$(dirname "$0")/scripts" && pwd)"
LOG_FILE="$(cd "$(dirname "$0")" && pwd)/results/logs/pipeline_$(date +%Y%m%d_%H%M%S).log"

# === 颜色输出 ===
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log() { echo -e "${BLUE}[$(date +%H:%M:%S)]${NC} $1" | tee -a "$LOG_FILE"; }
success() { echo -e "${GREEN}[DONE]${NC} $1" | tee -a "$LOG_FILE"; }
error() { echo -e "${RED}[ERROR]${NC} $1" | tee -a "$LOG_FILE"; }
warn() { echo -e "${YELLOW}[WARN]${NC} $1" | tee -a "$LOG_FILE"; }

# === 参数解析 ===
while [[ $# -gt 0 ]]; do
    case "$1" in
        --step) START_STEP="$2"; END_STEP="$2"; shift 2 ;;
        --threads) THREADS="$2"; shift 2 ;;
        --help|-h)
            echo "PHB_gtdb 主分析流程"
            echo "用法: bash main_pipeline.sh [选项]"
            echo "  --step N      仅运行第 N 步 (1-6)"
            echo "  --threads N   线程数 (默认: 30)"
            echo "  --help        显示帮助"
            exit 0
            ;;
        *) error "未知参数: $1"; exit 1 ;;
    esac
done

mkdir -p results/logs

log "============================================"
log "PHB_gtdb 分析流程开始"
log "============================================"
log "线程数: $THREADS"
log "步骤范围: $START_STEP - $END_STEP"
log "日志文件: $LOG_FILE"
log ""

# === 环境检查 ===
log "检查环境..."
python3 --version | tee -a "$LOG_FILE"
for tool in mafft fasttree hmmscan prodigal blastp; do
    if command -v "$tool" &>/dev/null; then
        log "  $tool: $(command -v "$tool")"
    else
        warn "  $tool: NOT FOUND"
    fi
done
log ""

# ==============================================================================
# Step 1: HMM 搜索 PHB 基因
# ==============================================================================
step1() {
    log "============================================"
    log "Step 1/6: HMM 搜索 PHB 相关基因"
    log "============================================"

    cd "$SCRIPTS_DIR"
    python3 01_phb_search.py --threads "$THREADS" --resume 2>&1 | tee -a "$LOG_FILE"

    if [ -f "../data/processed/phb_search_results.tsv" ]; then
        success "Step 1 完成"
        n_genomes=$(wc -l < ../data/processed/phb_search_results.tsv)
        log "含 PHB 基因的基因组数: $((n_genomes - 1))"
    else
        error "Step 1 失败: 输出文件不存在"
        exit 1
    fi
}

# ==============================================================================
# Step 2: 提取序列
# ==============================================================================
step2() {
    log "============================================"
    log "Step 2/6: 提取 PHB 基因序列"
    log "============================================"

    cd "$SCRIPTS_DIR"
    python3 02_extract_sequences.py --cdhit 0.95 2>&1 | tee -a "$LOG_FILE"

    if [ -f "../data/processed/phb_proteins_dedup.fasta" ] || \
       [ -f "../data/processed/phb_proteins_annotated.fasta" ]; then
        success "Step 2 完成"
    else
        error "Step 2 失败"
        exit 1
    fi
}

# ==============================================================================
# Step 3: 多序列比对
# ==============================================================================
step3() {
    log "============================================"
    log "Step 3/6: 多序列比对 (MAFFT)"
    log "============================================"

    cd "$SCRIPTS_DIR"
    python3 03_msa.py --gene all --threads "$THREADS" 2>&1 | tee -a "$LOG_FILE"

    # 检查是否有对齐文件生成
    n_aligned=$(find ../data/processed/ -name "*_aligned.fasta" 2>/dev/null | wc -l)
    if [ "$n_aligned" -gt 0 ]; then
        success "Step 3 完成 ($n_aligned 个比对文件)"
    else
        warn "Step 3: 未生成比对文件"
    fi
}

# ==============================================================================
# Step 4: 系统发育分析
# ==============================================================================
step4() {
    log "============================================"
    log "Step 4/6: 系统发育分析 (IQ-TREE)"
    log "============================================"

    cd "$SCRIPTS_DIR"
    python3 04_phylogeny.py --gene all --threads "$THREADS" 2>&1 | tee -a "$LOG_FILE"

    n_trees=$(find ../data/processed/ -name "*.treefile" 2>/dev/null | wc -l)
    if [ "$n_trees" -gt 0 ]; then
        success "Step 4 完成 ($n_trees 个系统发育树)"
    else
        warn "Step 4: 未生成系统发育树"
    fi
}

# ==============================================================================
# Step 5: 功能注释
# ==============================================================================
step5() {
    log "============================================"
    log "Step 5/6: 功能注释"
    log "============================================"

    cd "$SCRIPTS_DIR"
    python3 05_annotation.py --threads "$THREADS" 2>&1 | tee -a "$LOG_FILE"

    if [ -f "../data/processed/annotations/eggnog.emapper.annotations" ] || \
       [ -f "../results/tables/phb_phylum_summary.tsv" ]; then
        success "Step 5 完成"
    else
        warn "Step 5: 注释可能不完整 (eggNOG-mapper 可能未安装)"
    fi
}

# ==============================================================================
# Step 6: 可视化
# ==============================================================================
step6() {
    log "============================================"
    log "Step 6/6: 可视化"
    log "============================================"

    cd "$SCRIPTS_DIR"
    python3 06_visualization.py --all 2>&1 | tee -a "$LOG_FILE"

    n_figs=$(find ../results/figures/ -name "*.pdf" 2>/dev/null | wc -l)
    if [ "$n_figs" -gt 0 ]; then
        success "Step 6 完成 ($n_figs 个图表)"
        log "图表位置: results/figures/"
    else
        warn "Step 6: 未生成图表"
    fi
}

# ==============================================================================
# 执行流程
# ==============================================================================

for step_num in $(seq "$START_STEP" "$END_STEP"); do
    case "$step_num" in
        1) step1 ;;
        2) step2 ;;
        3) step3 ;;
        4) step4 ;;
        5) step5 ;;
        6) step6 ;;
        *) error "无效步骤: $step_num"; exit 1 ;;
    esac
    log ""
done

log "============================================"
log "PHB_gtdb 分析流程全部完成!"
log "============================================"
log "结果目录: $(cd "$(dirname "$0")" && pwd)/results/"
log ""
log "生成文件:"
find "$(cd "$(dirname "$0")" && pwd)/results/" -type f 2>/dev/null | while read f; do
    log "  $f"
done
