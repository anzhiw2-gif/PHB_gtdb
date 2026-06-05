#!/bin/bash
# =============================================================================
# PHB_gtdb — 搜索实时监控脚本
# 用法:
#   bash monitor_search.sh                  # 一次性查看 (自动检测搜索类型)
#   bash monitor_search.sh --watch 15       # 每15秒刷新
#   bash monitor_search.sh --watch 30 --mode archaea  # 指定搜索模式
# =============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="$SCRIPT_DIR/data"
RESULTS_DIR="$SCRIPT_DIR/results"
LOGS_DIR="$RESULTS_DIR/logs"
PROCESSED_DIR="$DATA_DIR/processed"

# 颜色
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
MAGENTA='\033[0;35m'
BOLD='\033[1m'
DIM='\033[2m'
NC='\033[0m'

# ==============================================================================
# 工具函数
# ==============================================================================
format_time() {
    local seconds=$1
    if [ "$seconds" -gt 86400 ]; then
        local d=$((seconds / 86400))
        local h=$(((seconds % 86400) / 3600))
        echo "${d}d ${h}h"
    elif [ "$seconds" -gt 3600 ]; then
        local h=$((seconds / 3600))
        local m=$(((seconds % 3600) / 60))
        echo "${h}h ${m}m"
    elif [ "$seconds" -gt 60 ]; then
        local m=$((seconds / 60))
        local s=$((seconds % 60))
        echo "${m}m ${s}s"
    else
        echo "${seconds}s"
    fi
}

format_size_bytes() {
    local bytes=$1
    if [ "$bytes" -gt 1073741824 ]; then
        awk "BEGIN {printf \"%.1f GB\", $bytes/1073741824}"
    elif [ "$bytes" -gt 1048576 ]; then
        awk "BEGIN {printf \"%.1f MB\", $bytes/1048576}"
    else
        awk "BEGIN {printf \"%.1f KB\", $bytes/1024}"
    fi
}

draw_bar() {
    local percent=$1
    local width=${2:-50}
    local bars=$(awk "BEGIN {printf \"%.0f\", $percent * $width / 100}")
    [ "$bars" -gt "$width" ] && bars=$width

    local bar=""
    for ((i=0; i<bars; i++)); do bar="${bar}█"; done
    for ((i=bars; i<width; i++)); do bar="${bar}░"; done

    # 颜色
    if awk "BEGIN {exit !($percent >= 90)}"; then
        echo -e "${GREEN}${bar}${NC}"
    elif awk "BEGIN {exit !($percent >= 50)}"; then
        echo -e "${YELLOW}${bar}${NC}"
    else
        echo -e "${BLUE}${bar}${NC}"
    fi
}

# ==============================================================================
# 检测搜索模式
# ==============================================================================
detect_mode() {
    # 检查是否有指定模式
    if [ -n "$SEARCH_MODE" ]; then
        echo "$SEARCH_MODE"
        return
    fi

    # 自动检测
    local has_bacteria=false
    local has_archaea=false

    # 检查进程
    if ps aux | grep -v grep | grep -q "01_phb_search"; then
        has_bacteria=true
    fi
    if ps aux | grep -v grep | grep -q "01b_archaea_search"; then
        has_archaea=true
    fi

    # 检查结果文件
    if [ -f "$PROCESSED_DIR/phb_search_results.tsv" ]; then
        has_bacteria=true
    fi
    if [ -f "$PROCESSED_DIR/archaea_phb_search_results.tsv" ]; then
        has_archaea=true
    fi

    if $has_bacteria && ! $has_archaea; then
        echo "bacteria"
    elif $has_archaea && ! $has_bacteria; then
        echo "archaea"
    elif $has_bacteria && $has_archaea; then
        echo "both"
    else
        echo "none"
    fi
}

# ==============================================================================
# 各模块
# ==============================================================================
show_header() {
    clear
    local mode=$(detect_mode)
    local mode_label=""
    case "$mode" in
        bacteria) mode_label="${GREEN}细菌 Bacteria${NC}" ;;
        archaea)  mode_label="${MAGENTA}古菌 Archaea${NC}" ;;
        both)     mode_label="${GREEN}细菌${NC}+${MAGENTA}古菌${NC}" ;;
        *)        mode_label="${DIM}未知${NC}" ;;
    esac

    echo -e "${BOLD}${CYAN}╔══════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BOLD}${CYAN}║${NC}     ${BOLD}PHB_gtdb PhaZ 搜索监控${NC} — $(date '+%Y-%m-%d %H:%M:%S')     ${BOLD}${CYAN}║${NC}"
    echo -e "${BOLD}${CYAN}║${NC}     目标: ${mode_label}  |  工具: Pyrodigal + DIAMOND          ${BOLD}${CYAN}║${NC}"
    echo -e "${BOLD}${CYAN}╚══════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
}

show_bacteria_progress() {
    echo -e "${BOLD}${GREEN}[B] 细菌搜索 — PhaZ 降解基因${NC}"
    echo -e "${GREEN}──────────────────────────────────────────────────────────────────────${NC}"

    local PID=$(ps aux | grep -v grep | grep "01_phb_search" | awk '{print $2}' | head -1)
    local RESULT_FILE="$PROCESSED_DIR/phb_search_results.tsv"
    local TMP_DIR="$PROCESSED_DIR/tmp"

    # 进程状态
    if [ -n "$PID" ]; then
        local CPU=$(ps -p "$PID" -o %cpu --no-headers 2>/dev/null | xargs)
        local MEM=$(ps -p "$PID" -o %mem --no-headers 2>/dev/null | xargs)
        local ELAPSED=$(ps -o etime= -p "$PID" 2>/dev/null | xargs)
        local THREADS=$(ps -p "$PID" -o nlwp --no-headers 2>/dev/null | xargs)

        echo -e "  状态: ${GREEN}● 运行中${NC}  PID: ${BOLD}$PID${NC}  "
        echo -e "  CPU: ${CPU:-?}%  |  内存: ${MEM:-?}%  |  线程: ${THREADS:-?}  |  已运行: ${ELAPSED:-?}"
    else
        echo -e "  状态: ${YELLOW}● 未运行${NC}"
    fi

    # 进度统计
    local COMPLETED=0
    local HITS=0
    local TOTAL_HITS=0

    if [ -f "$RESULT_FILE" ] && [ -s "$RESULT_FILE" ]; then
        COMPLETED=$(tail -n +2 "$RESULT_FILE" 2>/dev/null | wc -l)
        HITS=$(tail -n +2 "$RESULT_FILE" 2>/dev/null | awk -F'\t' '{if($3>0) print}' | wc -l)
        TOTAL_HITS=$(tail -n +2 "$RESULT_FILE" 2>/dev/null | awk -F'\t' '{s+=$3}END{print int(s)}')
    fi

    # 从 tmp 目录估算已处理数
    local TMP_COUNT=0
    if [ -d "$TMP_DIR" ]; then
        TMP_COUNT=$(find "$TMP_DIR" -maxdepth 1 -type d 2>/dev/null | wc -l)
    fi
    local PROCESSED=$((COMPLETED > TMP_COUNT ? COMPLETED : TMP_COUNT))
    local TOTAL_EST=165468  # 细菌基因组总数

    if [ "$PROCESSED" -gt 0 ]; then
        local PERCENT=$(awk "BEGIN {printf \"%.1f\", $PROCESSED * 100 / $TOTAL_EST}")
        [ "$(echo "$PERCENT > 100" | bc 2>/dev/null || echo 0)" -eq 1 ] && PERCENT="100.0"

        echo -e "  进度: ${BOLD}${PROCESSED}${NC} / ${TOTAL_EST} 基因组  (${BOLD}${PERCENT}%${NC})"
        draw_bar "$PERCENT" 50
        echo ""

        # 速度 & ETA
        if [ -n "$PID" ]; then
            local ELAPSED_RAW=$(ps -o etime= -p "$PID" 2>/dev/null | xargs)
            if [ -n "$ELAPSED_RAW" ]; then
                local TOTAL_SEC=$(echo "$ELAPSED_RAW" | awk -F'[-:]' '{
                    if (NF==4) { print $1*86400 + $2*3600 + $3*60 + $4 }
                    else if (NF==3) { print $1*3600 + $2*60 + $3 }
                    else { print 0 }
                }')
                if [ "$TOTAL_SEC" -gt 0 ] 2>/dev/null; then
                    local RATE=$(awk "BEGIN {printf \"%.1f\", $PROCESSED / $TOTAL_SEC}")
                    local REMAINING=$((TOTAL_EST - PROCESSED))
                    local ETA_SEC=$(awk "BEGIN {printf \"%.0f\", $REMAINING / ($RATE + 0.001)}")
                    echo -e "  速率: ${BOLD}${RATE}${NC} 基因组/秒  |  "
                    echo -e "  预计剩余: ${BOLD}$(format_time $ETA_SEC)${NC}  (${REMAINING} 基因组)"
                fi
            fi
        fi
    else
        echo -e "  进度: 尚未开始"
    fi

    echo -e "  命中: ${RED}${BOLD}${HITS}${NC} 个基因组含 PhaZ  |  共 ${RED}${BOLD}${TOTAL_HITS}${NC} 个 PhaZ 基因"
    echo ""
}

show_archaea_progress() {
    echo -e "${BOLD}${MAGENTA}[A] 古菌搜索 — 验证 phaZ 不存在假设${NC}"
    echo -e "${MAGENTA}──────────────────────────────────────────────────────────────────────${NC}"

    local PID=$(ps aux | grep -v grep | grep "01b_archaea_search" | awk '{print $2}' | head -1)
    local RESULT_FILE="$PROCESSED_DIR/archaea_phb_search_results.tsv"
    local TMP_DIR="$PROCESSED_DIR/tmp_archaea"

    # 进程状态
    if [ -n "$PID" ]; then
        local CPU=$(ps -p "$PID" -o %cpu --no-headers 2>/dev/null | xargs)
        local MEM=$(ps -p "$PID" -o %mem --no-headers 2>/dev/null | xargs)
        local ELAPSED=$(ps -o etime= -p "$PID" 2>/dev/null | xargs)
        local THREADS=$(ps -p "$PID" -o nlwp --no-headers 2>/dev/null | xargs)

        echo -e "  状态: ${GREEN}● 运行中${NC}  PID: ${BOLD}$PID${NC}"
        echo -e "  CPU: ${CPU:-?}%  |  内存: ${MEM:-?}%  |  线程: ${THREADS:-?}  |  已运行: ${ELAPSED:-?}"
    else
        echo -e "  状态: ${YELLOW}● 未运行${NC}"
    fi

    # 进度统计
    local COMPLETED=0
    local HITS=0
    local TOTAL_HITS=0

    if [ -f "$RESULT_FILE" ] && [ -s "$RESULT_FILE" ]; then
        COMPLETED=$(tail -n +2 "$RESULT_FILE" 2>/dev/null | wc -l)
        HITS=$(tail -n +2 "$RESULT_FILE" 2>/dev/null | awk -F'\t' '{if($3>0) print}' | wc -l)
        TOTAL_HITS=$(tail -n +2 "$RESULT_FILE" 2>/dev/null | awk -F'\t' '{s+=$3}END{print int(s)}')
    fi

    local TMP_COUNT=0
    if [ -d "$TMP_DIR" ]; then
        TMP_COUNT=$(find "$TMP_DIR" -maxdepth 1 -type d 2>/dev/null | wc -l)
    fi
    local PROCESSED=$((COMPLETED > TMP_COUNT ? COMPLETED : TMP_COUNT))
    local TOTAL_EST=34455  # 古菌基因组估算

    if [ "$PROCESSED" -gt 0 ]; then
        local PERCENT=$(awk "BEGIN {printf \"%.1f\", $PROCESSED * 100 / $TOTAL_EST}")
        [ "$(echo "$PERCENT > 100" | bc 2>/dev/null || echo 0)" -eq 1 ] && PERCENT="100.0"

        echo -e "  进度: ${BOLD}${PROCESSED}${NC} / ~${TOTAL_EST} 基因组  (${BOLD}${PERCENT}%${NC})"
        draw_bar "$PERCENT" 50
        echo ""

        if [ -n "$PID" ]; then
            local ELAPSED_RAW=$(ps -o etime= -p "$PID" 2>/dev/null | xargs)
            if [ -n "$ELAPSED_RAW" ]; then
                local TOTAL_SEC=$(echo "$ELAPSED_RAW" | awk -F'[-:]' '{
                    if (NF==4) { print $1*86400 + $2*3600 + $3*60 + $4 }
                    else if (NF==3) { print $1*3600 + $2*60 + $3 }
                    else { print 0 }
                }')
                if [ "$TOTAL_SEC" -gt 0 ] 2>/dev/null; then
                    local RATE=$(awk "BEGIN {printf \"%.1f\", $PROCESSED / $TOTAL_SEC}")
                    local REMAINING=$((TOTAL_EST - PROCESSED))
                    local ETA_SEC=$(awk "BEGIN {printf \"%.0f\", $REMAINING / ($RATE + 0.001)}")
                    echo -e "  速率: ${BOLD}${RATE}${NC} 基因组/秒  |  "
                    echo -e "  预计剩余: ${BOLD}$(format_time $ETA_SEC)${NC}  (${REMAINING} 基因组)"
                fi
            fi
        fi
    else
        echo -e "  进度: 尚未开始"
    fi

    echo -e "  命中: ${RED}${BOLD}${HITS}${NC} 个古菌基因组含 PhaZ  |  共 ${RED}${BOLD}${TOTAL_HITS}${NC} 个 PhaZ 同源序列"
    echo ""
}

show_recent_hits() {
    local N=${1:-10}
    echo -e "${BOLD}${BLUE}[最近命中] 最新发现 PhaZ 的基因组${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────────────────────${NC}"

    local BACTERIA_FILE="$PROCESSED_DIR/phb_search_results.tsv"
    local ARCHAEA_FILE="$PROCESSED_DIR/archaea_phb_search_results.tsv"

    local shown=0
    for f in "$ARCHAEA_FILE" "$BACTERIA_FILE"; do
        if [ -f "$f" ] && [ -s "$f" ]; then
            local label=$(basename "$f")
            echo -e "  ${DIM}来源: $label${NC}"

            # 取 phaZ_count > 0 的最后 N 条
            tail -n +2 "$f" 2>/dev/null | tac | awk -F'\t' '$3>0' | head -"$N" | while IFS=$'\t' read -r gid domain count evalue pident refs; do
                local eval_fmt=""
                if [ -n "$evalue" ] && [ "$evalue" != "best_evalue" ]; then
                    eval_fmt="e=${evalue}"
                fi
                local pid_fmt=""
                if [ -n "$pident" ] && [ "$pident" != "best_pident" ]; then
                    pid_fmt="id=${pident}%"
                fi
                printf "  ${RED}●${NC} %-30s  PhaZ: %2s  %s  %s\n" "$gid" "$count" "$pid_fmt" "$eval_fmt"
            done
            shown=1
            break
        fi
    done

    if [ "$shown" -eq 0 ]; then
        echo -e "  ${DIM}暂无命中 (或搜索尚未开始)${NC}"
    fi
    echo ""
}

show_log_tail() {
    local N=${1:-5}
    echo -e "${BOLD}${BLUE}[最新日志]${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────────────────────${NC}"

    local LOG=""
    if [ -f "$LOGS_DIR/01b_archaea_search.log" ]; then
        LOG="$LOGS_DIR/01b_archaea_search.log"
    elif [ -f "$LOGS_DIR/01_phb_search.log" ]; then
        LOG="$LOGS_DIR/01_phb_search.log"
    fi

    if [ -n "$LOG" ]; then
        echo -e "  ${DIM}$(basename $LOG)${NC}"
        tail -"$N" "$LOG" | while IFS= read -r line; do
            echo -e "  ${DIM}$line${NC}"
        done
    else
        echo -e "  ${DIM}暂无日志${NC}"
    fi
    echo ""
}

show_disk() {
    echo -e "${BOLD}${BLUE}[系统资源]${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────────────────────${NC}"

    local DISK_TOTAL=$(df -h /home/data 2>/dev/null | awk 'NR==2 {print $2}')
    local DISK_USED=$(df -h /home/data 2>/dev/null | awk 'NR==2 {print $3}')
    local DISK_AVAIL=$(df -h /home/data 2>/dev/null | awk 'NR==2 {print $4}')
    local DISK_PCT=$(df -h /home/data 2>/dev/null | awk 'NR==2 {print $5}')

    echo -e "  磁盘: 总量 ${BOLD}${DISK_TOTAL:-?}${NC}  已用 ${YELLOW}${DISK_USED:-?}${NC} (${DISK_PCT:-?})  可用 ${GREEN}${DISK_AVAIL:-?}${NC}"

    # 数据目录大小
    if [ -d "$DATA_DIR" ]; then
        local DATA_SIZE=$(du -sh "$DATA_DIR" 2>/dev/null | cut -f1)
        echo -e "  数据目录: ${BOLD}${DATA_SIZE:-?}${NC}"
    fi

    # CPU 负载
    local LOAD=$(uptime 2>/dev/null | awk -F'load average:' '{print $2}' | xargs)
    echo -e "  系统负载: ${LOAD:-?}"

    echo ""
}

show_summary() {
    echo -e "${BOLD}${CYAN}[汇总]${NC}"
    echo -e "${CYAN}──────────────────────────────────────────────────────────────────────${NC}"

    local B_FILE="$PROCESSED_DIR/phb_search_results.tsv"
    local A_FILE="$PROCESSED_DIR/archaea_phb_search_results.tsv"

    if [ -f "$B_FILE" ] && [ -s "$B_FILE" ]; then
        local B_TOTAL=$(tail -n +2 "$B_FILE" 2>/dev/null | wc -l)
        local B_HITS=$(tail -n +2 "$B_FILE" 2>/dev/null | awk -F'\t' '$3>0' | wc -l)
        local B_PHAZ=$(tail -n +2 "$B_FILE" 2>/dev/null | awk -F'\t' '{s+=$3}END{print int(s)}')
        echo -e "  ${GREEN}细菌搜索${NC}: 已处理 ${BOLD}${B_TOTAL}${NC} 基因组  |  命中 ${RED}${B_HITS}${NC} 含 PhaZ (${B_PHAZ} genes)"
    fi

    if [ -f "$A_FILE" ] && [ -s "$A_FILE" ]; then
        local A_TOTAL=$(tail -n +2 "$A_FILE" 2>/dev/null | wc -l)
        local A_HITS=$(tail -n +2 "$A_FILE" 2>/dev/null | awk -F'\t' '$3>0' | wc -l)
        local A_PHAZ=$(tail -n +2 "$A_FILE" 2>/dev/null | awk -F'\t' '{s+=$3}END{print int(s)}')
        echo -e "  ${MAGENTA}古菌搜索${NC}: 已处理 ${BOLD}${A_TOTAL}${NC} 基因组  |  命中 ${RED}${A_HITS}${NC} 含 PhaZ (${A_PHAZ} genes)"

        if [ "$A_HITS" -eq 0 ] && [ "$A_TOTAL" -gt 0 ]; then
            echo -e "  ${GREEN}✓ 古菌中未发现 phaZ — 与预期一致${NC}"
        elif [ "$A_HITS" -gt 0 ]; then
            echo -e "  ${RED}⚠ 古菌中发现 phaZ 同源序列 — 需进一步验证${NC}"
        fi
    fi

    echo ""
}

# ==============================================================================
# 主逻辑
# ==============================================================================
WATCH_MODE=false
INTERVAL=15
SEARCH_MODE=""
N_RECENT=10
N_LOG=5

while [[ $# -gt 0 ]]; do
    case "$1" in
        --watch|-w)
            WATCH_MODE=true
            INTERVAL="${2:-15}"
            shift 2
            ;;
        --mode|-m)
            SEARCH_MODE="$2"
            shift 2
            ;;
        --recent|-r)
            N_RECENT="$2"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

if $WATCH_MODE; then
    while true; do
        show_header
        local mode=$(detect_mode)
        case "$mode" in
            bacteria|both)
                show_bacteria_progress
                ;;
        esac
        case "$mode" in
            archaea|both)
                show_archaea_progress
                ;;
        esac
        if [ "$mode" = "none" ]; then
            echo -e "  ${YELLOW}未检测到运行中的搜索进程${NC}"
            echo -e "  启动方式:"
            echo -e "    ${GREEN}细菌搜索${NC}: nohup conda run -n phb_gtdb python 01_phb_search.py --threads 30 &"
            echo -e "    ${MAGENTA}古菌搜索${NC}: nohup conda run -n phb_gtdb python 01b_archaea_search.py --threads 30 &"
            echo ""
        fi
        show_recent_hits "$N_RECENT"
        show_log_tail "$N_LOG"
        show_disk
        show_summary
        echo -e "${CYAN}每 ${INTERVAL} 秒刷新 | Ctrl+C 退出 | $(date '+%H:%M:%S')${NC}"
        sleep "$INTERVAL"
    done
else
    show_header
    mode=$(detect_mode)
    case "$mode" in
        bacteria|both)
            show_bacteria_progress
            ;;
    esac
    case "$mode" in
        archaea|both)
            show_archaea_progress
            ;;
    esac
    if [ "$mode" = "none" ]; then
        echo -e "${YELLOW}未检测到运行中的搜索进程${NC}"
        echo ""
    fi
    show_recent_hits "$N_RECENT"
    show_log_tail "$N_LOG"
    show_disk
    show_summary
fi
