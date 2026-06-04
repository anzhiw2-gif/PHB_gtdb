#!/bin/bash
# =============================================================================
# PHB_gtdb 实时监控脚本
# 用法:
#   bash monitor.sh              # 一次性查看
#   bash monitor.sh --watch 10   # 每10秒刷新
#   bash monitor.sh --watch 5    # 每5秒刷新
#   bash monitor.sh --log 20     # 查看最近20行日志
# =============================================================================

PROJECT_DIR="/home/data/haoyu/PHB_gtdb"
RESULT_FILE="$PROJECT_DIR/data/processed/phb_search_results.tsv"
TMP_DIR="$PROJECT_DIR/data/processed/tmp"
LOG_DIR="$PROJECT_DIR/results/logs"

# === 颜色 ===
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

print_header() {
    echo -e "${BOLD}${CYAN}╔══════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BOLD}${CYAN}║${NC}     ${BOLD}PHB_gtdb 实时监控 — $(date '+%Y-%m-%d %H:%M:%S')${NC}     ${BOLD}${CYAN}║${NC}"
    echo -e "${BOLD}${CYAN}╚══════════════════════════════════════════════════════════════╝${NC}"
    echo ""
}

# --- 1. 运行进程 ---
show_processes() {
    echo -e "${BOLD}${BLUE}[1] 运行中的进程${NC}"
    echo -e "${BLUE}────────────────────────────────────────────────────────────${NC}"

    MAIN_PID=$(ps aux | grep -v grep | grep '01_phb_search.py' | awk '{print $2}' | head -1)

    if [ -z "$MAIN_PID" ]; then
        echo -e "  ${YELLOW}⚠ 01_phb_search.py 未在运行${NC}"
    else
        MAIN_CPU=$(ps aux | grep -v grep | grep '01_phb_search.py' | awk '{print $3}' | head -1)
        MAIN_MEM=$(ps aux | grep -v grep | grep '01_phb_search.py' | awk '{print $4}' | head -1)
        MAIN_ELAPSED=$(ps -o etime= -p "$MAIN_PID" 2>/dev/null | xargs)

        echo -e "  ${GREEN}●${NC} 主进程 PID: ${BOLD}$MAIN_PID${NC}  CPU: ${MAIN_CPU}%  MEM: ${MAIN_MEM}%  运行时间: ${MAIN_ELAPSED}"

        # 统计子进程
        N_WORKERS=$(ps aux | grep -v grep | grep '01_phb_search.py' | grep -v "$MAIN_PID" | wc -l)
        echo -e "  ${GREEN}●${NC} 工作进程数: ${BOLD}$N_WORKERS${NC}"

        # 统计工具进程
        N_PRODIGAL=$(ps aux | grep -v grep | grep 'prodigal' | wc -l)
        N_BLAST=$(ps aux | grep -v grep | grep 'blastp' | wc -l)
        echo -e "  ${GREEN}●${NC} Prodigal 运行中: ${YELLOW}$N_PRODIGAL${NC}  |  BLASTP 运行中: ${YELLOW}$N_BLAST${NC}"
    fi
    echo ""
}

# --- 2. 基因组处理进度 ---
show_progress() {
    echo -e "${BOLD}${BLUE}[2] 基因组处理进度${NC}"
    echo -e "${BLUE}────────────────────────────────────────────────────────────${NC}"

    # 结果文件中的基因组数
    if [ -f "$RESULT_FILE" ]; then
        N_RESULTS=$(tail -n +2 "$RESULT_FILE" 2>/dev/null | wc -l)
        if [ "$N_RESULTS" -gt 0 ]; then
            N_PHAZ=$(tail -n +2 "$RESULT_FILE" 2>/dev/null | awk -F'\t' '{sum+=$2} END {print sum}')
            echo -e "  结果文件: ${BOLD}$N_RESULTS${NC} 个基因组 (含 ${GREEN}$N_PHAZ${NC} 个 PhaZ 基因)"
        else
            echo -e "  结果文件: 仅有表头，尚无命中"
        fi
    else
        echo -e "  结果文件: ${YELLOW}尚未创建${NC}"
    fi

    # 临时目录 (正在处理或刚完成的)
    if [ -d "$TMP_DIR" ]; then
        N_TMP=$(ls -d "$TMP_DIR"/*/ 2>/dev/null | wc -l)
        if [ "$N_TMP" -gt 0 ]; then
            echo -e "  活跃临时目录: ${BOLD}$N_TMP${NC} 个 (正在处理/待清理)"

            # 显示最新的3个
            echo -e "  ${CYAN}最近活跃:${NC}"
            ls -lt "$TMP_DIR" 2>/dev/null | head -4 | tail -3 | while read line; do
                DIRNAME=$(echo "$line" | awk '{print $NF}')
                # 检查目录中有哪些文件
                HAS_FNA=""
                HAS_FAA=""
                HAS_BLAST=""
                [ -f "$TMP_DIR/$DIRNAME"/*.fna ] 2>/dev/null && HAS_FNA="DNA" || HAS_FNA=""
                [ -f "$TMP_DIR/$DIRNAME"/*.faa ] 2>/dev/null && HAS_FAA="Prodigal✓" || HAS_FAA=""
                [ -f "$TMP_DIR/$DIRNAME"/*.blast ] 2>/dev/null && HAS_BLAST="BLAST✓" || HAS_BLAST=""
                STATUS="$HAS_FNA $HAS_FAA $HAS_BLAST"
                echo -e "    ${CYAN}$DIRNAME${NC} → $STATUS"
            done
        else
            echo -e "  活跃临时目录: 0 (全部已清理)"
        fi
    fi
    echo ""
}

# --- 3. 结果统计 ---
show_results() {
    echo -e "${BOLD}${BLUE}[3] PhaZ 搜索结果统计${NC}"
    echo -e "${BLUE}────────────────────────────────────────────────────────────${NC}"

    if [ -f "$RESULT_FILE" ] && [ "$(tail -n +2 "$RESULT_FILE" 2>/dev/null | wc -l)" -gt 0 ]; then
        echo -e "  ${BOLD}Top 20 PhaZ 命中基因组:${NC}"
        echo -e "  ${CYAN}$(printf '%-30s %s' 'Genome_ID' 'PhaZ_Count')${NC}"
        echo -e "  ${CYAN}$(printf '%-30s %s' '------------------------------' '----------')${NC}"
        tail -n +2 "$RESULT_FILE" 2>/dev/null | sort -t$'\t' -k2 -rn | head -20 | while IFS=$'\t' read -r gid count; do
            if [ "$count" -gt 0 ]; then
                printf "  ${GREEN}%-30s${NC} %s\n" "$gid" "$count"
            fi
        done

        TOTAL_GENOMES=$(tail -n +2 "$RESULT_FILE" 2>/dev/null | wc -l)
        TOTAL_PHAZ=$(tail -n +2 "$RESULT_FILE" 2>/dev/null | awk -F'\t' '{sum+=$2} END {print sum}')
        WITH_PHAZ=$(tail -n +2 "$RESULT_FILE" 2>/dev/null | awk -F'\t' '{if($2>0) print}' | wc -l)
        echo ""
        echo -e "  ${BOLD}汇总:${NC} 共 ${BOLD}$TOTAL_GENOMES${NC} 基因组, ${GREEN}$WITH_PHAZ${NC} 含 PhaZ, 共 ${GREEN}$TOTAL_PHAZ${NC} 个 PhaZ 基因"
    else
        echo -e "  ${YELLOW}暂无结果 (搜索进行中...)${NC}"
    fi
    echo ""
}

# --- 4. 最新日志 ---
show_log() {
    local N=${1:-10}
    echo -e "${BOLD}${BLUE}[4] 最新日志 (最近 $N 行)${NC}"
    echo -e "${BLUE}────────────────────────────────────────────────────────────${NC}"

    # 找最新日志文件
    LATEST_LOG=$(ls -t "$LOG_DIR"/01_phb_search_*.log 2>/dev/null | head -1)

    if [ -z "$LATEST_LOG" ]; then
        echo -e "  ${YELLOW}无日志文件${NC}"
    else
        echo -e "  ${CYAN}日志: $(basename "$LATEST_LOG")${NC}"
        echo ""
        tail -"$N" "$LATEST_LOG" | while IFS= read -r line; do
            if echo "$line" | grep -q "ERROR"; then
                echo -e "  ${RED}$line${NC}"
            elif echo "$line" | grep -q "WARNING\|WARN"; then
                echo -e "  ${YELLOW}$line${NC}"
            elif echo "$line" | grep -q "完成\|成功\|SUCCESS"; then
                echo -e "  ${GREEN}$line${NC}"
            else
                echo -e "  $line"
            fi
        done
    fi
    echo ""
}

# --- 5. 磁盘使用 ---
show_disk() {
    echo -e "${BOLD}${BLUE}[5] 磁盘使用${NC}"
    echo -e "${BLUE}────────────────────────────────────────────────────────────${NC}"

    PROJECT_SIZE=$(du -sh "$PROJECT_DIR" 2>/dev/null | cut -f1)
    TMP_SIZE=$(du -sh "$TMP_DIR" 2>/dev/null | cut -f1)
    RESULTS_SIZE=$(du -sh "$PROJECT_DIR/results" 2>/dev/null | cut -f1)

    echo -e "  项目总大小: ${BOLD}$PROJECT_SIZE${NC}"
    echo -e "  tmp/ 目录:  ${YELLOW}$TMP_SIZE${NC}"
    echo -e "  results/:   $RESULTS_SIZE"
    echo ""
}

# ==============================================================================
# 主逻辑
# ==============================================================================

WATCH_MODE=false
INTERVAL=10
LOG_LINES=10

while [[ $# -gt 0 ]]; do
    case "$1" in
        --watch|-w)
            WATCH_MODE=true
            INTERVAL="${2:-10}"
            shift 2
            ;;
        --log|-l)
            LOG_LINES="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

if $WATCH_MODE; then
    # 持续刷新模式
    while true; do
        clear
        print_header
        show_processes
        show_progress
        show_results
        show_log "$LOG_LINES"
        show_disk
        echo -e "${CYAN}每 ${INTERVAL} 秒刷新 | Ctrl+C 退出${NC}"
        sleep "$INTERVAL"
    done
else
    # 单次查看
    print_header
    show_processes
    show_progress
    show_results
    show_log "$LOG_LINES"
    show_disk
fi
