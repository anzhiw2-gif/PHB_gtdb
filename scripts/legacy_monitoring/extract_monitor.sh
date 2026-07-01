#!/bin/bash
# =============================================================================
# GTDB 解压监控脚本
# 用法:
#   bash extract_monitor.sh              # 一次性查看
#   bash extract_monitor.sh --watch 30   # 每30秒自动刷新
# =============================================================================

EXTRACT_DIR="/home/data/haoyu/GTDB/gtdb_genomes_reps_r232"
ARCHIVE="/home/data/haoyu/GTDB/gtdb_genomes_reps_r232.tar.gz"
ARCHIVE_SIZE=192164247525  # 179 GB
LOGFILE="/home/data/haoyu/GTDB/extract.log"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

format_size() {
    local bytes=$1
    if [ "$bytes" -gt 1073741824 ]; then
        awk "BEGIN {printf \"%.2f GB\", $bytes/1073741824}"
    elif [ "$bytes" -gt 1048576 ]; then
        awk "BEGIN {printf \"%.2f MB\", $bytes/1048576}"
    else
        awk "BEGIN {printf \"%.2f KB\", $bytes/1024}"
    fi
}

print_header() {
    clear
    echo -e "${BOLD}${CYAN}╔══════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BOLD}${CYAN}║${NC}      ${BOLD}GTDB R232 解压监控 — $(date '+%Y-%m-%d %H:%M:%S')${NC}      ${BOLD}${CYAN}║${NC}"
    echo -e "${BOLD}${CYAN}╚══════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
}

show_status() {
    echo -e "${BOLD}${BLUE}[1] 解压状态${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────────────────────${NC}"

    TAR_PID=$(ps aux | grep -v grep | grep 'tar.*gtdb_genomes_reps' | awk '{print $2}' | head -1)

    if [ -n "$TAR_PID" ]; then
        CPU=$(ps -p "$TAR_PID" -o %cpu --no-headers 2>/dev/null | xargs)
        ELAPSED=$(ps -o etime= -p "$TAR_PID" 2>/dev/null | xargs)
        echo -e "  ${GREEN}●${NC} tar 解压中  PID: ${BOLD}$TAR_PID${NC}  CPU: ${CPU}%  运行: ${ELAPSED}"
    else
        echo -e "  ${GREEN}●${NC} tar 已结束${NC}"
    fi
    echo ""
}

show_progress() {
    echo -e "${BOLD}${BLUE}[2] 解压进度${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────────────────────${NC}"

    # 统计基因组数
    if [ -d "$EXTRACT_DIR" ]; then
        GCA_COUNT=$(find "$EXTRACT_DIR/database/GCA" -name '*_genomic.fna.gz' 2>/dev/null | wc -l)
        GCF_COUNT=$(find "$EXTRACT_DIR/database/GCF" -name '*_genomic.fna.gz' 2>/dev/null | wc -l)
        TOTAL_EXT=$((GCA_COUNT + GCF_COUNT))
    else
        TOTAL_EXT=0
    fi

    # R232 代表基因组约 158,000+ 个
    EST_TOTAL=158500

    if [ "$TOTAL_EXT" -gt 0 ]; then
        PERCENT=$(awk "BEGIN {printf \"%.1f\", $TOTAL_EXT * 100 / $EST_TOTAL}")
        BARS=$(awk "BEGIN {printf \"%.0f\", $PERCENT * 50 / 100}")
        [ "$BARS" -gt 50 ] && BARS=50

        BAR=""
        for ((i=0; i<BARS; i++)); do BAR="${BAR}█"; done
        for ((i=BARS; i<50; i++)); do BAR="${BAR}░"; done

        echo -e "  ${BOLD}已解压基因组:${NC} ${TOTAL_EXT} (GCA: ${GCA_COUNT}, GCF: ${GCF_COUNT})"
        echo ""
        echo -e "  ${GREEN}${BAR}${NC}  ${BOLD}${PERCENT}%${NC}"
        echo ""

        # 速率
        TAR_PID=$(ps aux | grep -v grep | grep 'tar.*gtdb_genomes_reps' | awk '{print $2}' | head -1)
        if [ -n "$TAR_PID" ]; then
            ELAPSED_RAW=$(ps -o etime= -p "$TAR_PID" 2>/dev/null | xargs)
            if [ -n "$ELAPSED_RAW" ]; then
                TOTAL_SEC=$(echo "$ELAPSED_RAW" | awk -F'[-:]' '{
                    if (NF==4) { print $1*86400 + $2*3600 + $3*60 + $4 }
                    else if (NF==3) { print $1*3600 + $2*60 + $3 }
                    else { print 0 }
                }')
                if [ "$TOTAL_SEC" -gt 0 ] 2>/dev/null; then
                    RATE=$(awk "BEGIN {printf \"%.0f\", $TOTAL_EXT / $TOTAL_SEC}")
                    REMAINING=$((EST_TOTAL - TOTAL_EXT))
                    ETA_SEC=$(awk "BEGIN {printf \"%.0f\", $REMAINING / ($RATE + 0.01)}")
                    ETA_H=$((ETA_SEC / 3600))
                    ETA_M=$(((ETA_SEC % 3600) / 60))
                    echo -e "  ${BOLD}速率:${NC}   ${RATE} 基因组/秒"
                    echo -e "  ${BOLD}预计剩余:${NC} ${ETA_H}时${ETA_M}分 (${REMAINING} 基因组)"
                fi
            fi
        fi
    else
        echo -e "  ${YELLOW}尚未开始解压${NC}"
    fi
    echo ""
}

show_disk() {
    echo -e "${BOLD}${BLUE}[3] 磁盘空间${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────────────────────${NC}"

    DISK_TOTAL=$(df -h / | awk 'NR==2 {print $2}')
    DISK_USED=$(df -h / | awk 'NR==2 {print $3}')
    DISK_AVAIL=$(df -h / | awk 'NR==2 {print $4}')
    DISK_PCT=$(df -h / | awk 'NR==2 {print $5}')

    echo -e "  总容量: ${BOLD}${DISK_TOTAL}${NC}"
    echo -e "  已使用: ${YELLOW}${DISK_USED}${NC} (${DISK_PCT})"
    echo -e "  可用:   ${GREEN}${DISK_AVAIL}${NC}"

    if [ -d "$EXTRACT_DIR" ]; then
        DIR_SIZE=$(du -sh "$EXTRACT_DIR" 2>/dev/null | cut -f1)
        echo -e "  解压目录: ${BOLD}${DIR_SIZE}${NC}"
    fi

    echo ""
}

show_latest() {
    echo -e "${BOLD}${BLUE}[4] 最近解压的基因组${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────────────────────${NC}"

    if [ -d "$EXTRACT_DIR/database" ]; then
        find "$EXTRACT_DIR/database" -name '*_genomic.fna.gz' -newer "$EXTRACT_DIR" 2>/dev/null | tail -10 | while read f; do
            echo -e "  ${CYAN}$(basename $f)${NC}"
        done
    fi
    echo ""
}

# ==============================================================================
WATCH_MODE=false
INTERVAL=30

while [[ $# -gt 0 ]]; do
    case "$1" in
        --watch|-w)
            WATCH_MODE=true
            INTERVAL="${2:-30}"
            shift 2
            ;;
        *)
            shift
            ;;
    esac
done

if $WATCH_MODE; then
    while true; do
        print_header
        show_status
        show_progress
        show_disk
        show_latest
        echo -e "${CYAN}每 ${INTERVAL} 秒刷新 | Ctrl+C 退出${NC}"
        sleep "$INTERVAL"
    done
else
    print_header
    show_status
    show_progress
    show_disk
    show_latest
fi
