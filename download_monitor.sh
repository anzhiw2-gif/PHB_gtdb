#!/bin/bash
# =============================================================================
# GTDB 下载监控脚本
# 用法:
#   bash download_monitor.sh              # 一次性查看
#   bash download_monitor.sh --watch 10   # 每10秒自动刷新
# =============================================================================

DOWNLOAD_DIR="/home/data/haoyu/GTDB"
TARGET_FILE="$DOWNLOAD_DIR/gtdb_genomes_reps_r232.tar.gz"
LOG_FILE="$DOWNLOAD_DIR/download.log"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

TOTAL_SIZE=192164247525  # 179 GB

# --- 格式化字节 ---
format_size() {
    local bytes=$1
    if [ "$bytes" -gt 1073741824 ]; then
        echo "$(awk "BEGIN {printf \"%.2f\", $bytes/1073741824}") GB"
    elif [ "$bytes" -gt 1048576 ]; then
        echo "$(awk "BEGIN {printf \"%.2f\", $bytes/1048576}") MB"
    else
        echo "$(awk "BEGIN {printf \"%.2f\", $bytes/1024}") KB"
    fi
}

print_header() {
    echo -e "${BOLD}${CYAN}╔══════════════════════════════════════════════════════════════════════╗${NC}"
    echo -e "${BOLD}${CYAN}║${NC}       ${BOLD}GTDB R232 下载监控 — $(date '+%Y-%m-%d %H:%M:%S')${NC}       ${BOLD}${CYAN}║${NC}"
    echo -e "${BOLD}${CYAN}╚══════════════════════════════════════════════════════════════════════╝${NC}"
    echo ""
}

# --- 下载进度 ---
show_progress() {
    echo -e "${BOLD}${BLUE}[1] 下载状态${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────────────────────${NC}"

    # 检查下载进程 (aria2c 或 wget)
    DOWNLOAD_PID=$(ps aux | grep -v grep | grep -E "aria2c.*gtdb_genomes_reps|wget.*gtdb_genomes_reps" | awk '{print $2}' | head -1)
    DOWNLOAD_CMD=$(ps aux | grep -v grep | grep -E "aria2c.*gtdb_genomes_reps|wget.*gtdb_genomes_reps" | awk '{for(i=11;i<=NF;i++) printf "%s ", $i}' | head -1)

    if [ -n "$DOWNLOAD_PID" ]; then
        DOWNLOAD_CPU=$(ps -p "$DOWNLOAD_PID" -o %cpu --no-headers 2>/dev/null | xargs)
        DOWNLOAD_ELAPSED=$(ps -o etime= -p "$DOWNLOAD_PID" 2>/dev/null | xargs)
        TOOL_NAME=""
        if echo "$DOWNLOAD_CMD" | grep -q "aria2c"; then
            TOOL_NAME="aria2c"
            CONNS=$(echo "$DOWNLOAD_CMD" | grep -oP '\-x\s+\d+' | awk '{print $2}')
        elif echo "$DOWNLOAD_CMD" | grep -q "wget"; then
            TOOL_NAME="wget"
            CONNS="1"
        fi
        echo -e "  ${GREEN}●${NC} ${TOOL_NAME} 运行中  连接数: ${CONNS:-?}  PID: ${BOLD}$DOWNLOAD_PID${NC}  已运行: ${DOWNLOAD_ELAPSED}"
    else
        if [ -f "$TARGET_FILE" ]; then
            CURRENT_SIZE=$(stat --printf="%s" "$TARGET_FILE" 2>/dev/null || echo 0)
            if [ "$CURRENT_SIZE" -ge "$TOTAL_SIZE" ]; then
                echo -e "  ${GREEN}●${NC} 下载完成! ✓"
            else
                echo -e "  ${YELLOW}●${NC} wget 未运行 (可能中断，续传需手动重启)"
            fi
        else
            echo -e "  ${YELLOW}●${NC} 尚未开始下载"
        fi
    fi

    echo ""
}

# --- 文件大小和进度条 ---
show_filesize() {
    echo -e "${BOLD}${BLUE}[2] 下载进度${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────────────────────${NC}"

    # aria2c 日志解析 (更准确的进度)
    ARIA_LOG="$DOWNLOAD_DIR/download_aria2.log"
    if [ -f "$ARIA_LOG" ]; then
        # 提取最后一行进度信息 (以 [# 开头的行)
        LAST_LINE=$(grep '^\[' "$ARIA_LOG" 2>/dev/null | grep 'GiB/' | tail -1)
        if [ -z "$LAST_LINE" ]; then
            LAST_LINE=$(grep '^\[' "$ARIA_LOG" 2>/dev/null | tail -1)
        fi
        CURRENT_GB=$(echo "$LAST_LINE" | grep -oP '\[\S+\s+\K[\d.]+(?=GiB/)' 2>/dev/null)
        TOTAL_GB_TXT=$(echo "$LAST_LINE" | grep -oP '/\K[\d.]+(?=GiB\()' 2>/dev/null)
        PERCENT=$(echo "$LAST_LINE" | grep -oP '\(\K[\d.]+(?=%\))' 2>/dev/null)
        DL_SPEED=$(echo "$LAST_LINE" | grep -oP 'DL:\K\S+' 2>/dev/null)
        ETA_TIME=$(echo "$LAST_LINE" | grep -oP 'ETA:\K\S+' 2>/dev/null)

        if [ -n "$CURRENT_GB" ] && [ -n "$TOTAL_GB_TXT" ]; then
            CURRENT_SIZE=$(awk "BEGIN {printf \"%.0f\", $CURRENT_GB * 1073741824}")
            TOTAL_ARIA=$(awk "BEGIN {printf \"%.0f\", $TOTAL_GB_TXT * 1073741824}")
            PERCENT_NUM=$(printf "%.2f" "${PERCENT:-0}")

            BARS=$(awk "BEGIN {printf \"%.0f\", ${PERCENT:-0} * 50 / 100}")
            BAR=""
            for ((i=0; i<BARS; i++)); do BAR="${BAR}█"; done
            for ((i=BARS; i<50; i++)); do BAR="${BAR}░"; done

            echo -e "  ${BOLD}目标:${NC}   ${TOTAL_GB_TXT} GB"
            echo -e "  ${BOLD}已下载:${NC} ${CURRENT_GB} GB"
            echo -e "  ${BOLD}速度:${NC}   ${DL_SPEED:-?}"
            echo -e "  ${BOLD}预计剩余:${NC} ${ETA_TIME:-?}"
            echo ""
            echo -e "  ${GREEN}${BAR}${NC}  ${BOLD}${PERCENT_NUM}%${NC}"
            echo ""
            return
        fi
    fi

    # 回退到文件大小检测
    if [ -f "$TARGET_FILE" ]; then
        CURRENT_SIZE=$(stat --printf="%s" "$TARGET_FILE" 2>/dev/null || echo 0)
    else
        CURRENT_SIZE=0
    fi

    if [ "$CURRENT_SIZE" -gt 0 ]; then
        PERCENT=$(awk "BEGIN {printf \"%.2f\", $CURRENT_SIZE * 100 / $TOTAL_SIZE}")
        REMAINING=$((TOTAL_SIZE - CURRENT_SIZE))

        # 进度条 (50 字符)
        BARS=$(awk "BEGIN {printf \"%.0f\", $PERCENT * 50 / 100}")
        BAR=""
        for ((i=0; i<BARS; i++)); do BAR="${BAR}█"; done
        for ((i=BARS; i<50; i++)); do BAR="${BAR}░"; done

        echo -e "  ${BOLD}目标:${NC}  $(format_size $TOTAL_SIZE)"
        echo -e "  ${BOLD}已下载:${NC} $(format_size $CURRENT_SIZE)"
        echo -e "  ${BOLD}剩余:${NC}   $(format_size $REMAINING)"
        echo ""
        echo -e "  ${GREEN}${BAR}${NC}  ${BOLD}${PERCENT}%${NC}"
        echo ""

        # 预估速度 (支持 aria2c / wget)
        DOWNLOAD_PID=$(ps aux | grep -v grep | grep -E "aria2c.*gtdb|wget.*gtdb" | awk '{print $2}' | head -1)
        if [ -n "$DOWNLOAD_PID" ] && [ "$CURRENT_SIZE" -gt 1048576 ]; then
            DL_ELAPSED_RAW=$(ps -o etime= -p "$DOWNLOAD_PID" 2>/dev/null | xargs)
            if [ -n "$DL_ELAPSED_RAW" ]; then
                # 使用 awk 解析已用时间
                TOTAL_SEC=$(echo "$DL_ELAPSED_RAW" | awk -F'[-:]' '{
                    if (NF==4) { print $1*86400 + $2*3600 + $3*60 + $4 }
                    else if (NF==3) { print $1*3600 + $2*60 + $3 }
                    else { print 0 }
                }')

                if [ "$TOTAL_SEC" -gt 0 ] 2>/dev/null; then
                    SPEED_BPS=$(awk "BEGIN {printf \"%.0f\", $CURRENT_SIZE / $TOTAL_SEC}")
                    SPEED=$(format_size $SPEED_BPS)
                    ETA_SEC=$(awk "BEGIN {printf \"%.0f\", $REMAINING / $SPEED_BPS}")
                    ETA_H=$((ETA_SEC / 3600))
                    ETA_M=$(((ETA_SEC % 3600) / 60))
                    echo -e "  ${BOLD}速度:${NC}   ${SPEED}/s"
                    echo -e "  ${BOLD}预计剩余:${NC} ${ETA_H}时${ETA_M}分"
                fi
            fi
        fi
    else
        echo -e "  ${YELLOW}文件尚未开始下载"
    fi

    echo ""
}

# --- 磁盘使用 ---
show_disk() {
    echo -e "${BOLD}${BLUE}[3] 磁盘空间${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────────────────────${NC}"

    DISK_TOTAL=$(df -h / | awk 'NR==2 {print $2}')
    DISK_USED=$(df -h / | awk 'NR==2 {print $3}')
    DISK_AVAIL=$(df -h / | awk 'NR==2 {print $4}')
    DISK_PCT=$(df -h / | awk 'NR==2 {print $5}')

    echo -e "  总容量: ${BOLD}$DISK_TOTAL${NC}"
    echo -e "  已使用: ${YELLOW}$DISK_USED${NC} (${DISK_PCT})"
    echo -e "  可用:   ${GREEN}$DISK_AVAIL${NC}"

    # 下载目录
    if [ -d "$DOWNLOAD_DIR" ]; then
        DL_SIZE=$(du -sh "$DOWNLOAD_DIR" 2>/dev/null | cut -f1)
        echo -e "  GTDB/:   ${BOLD}$DL_SIZE${NC}"
    fi

    echo ""
}

# --- 最新日志 ---
show_log() {
    local N=${1:-5}
    echo -e "${BOLD}${BLUE}[4] 下载日志 (最近 $N 行)${NC}"
    echo -e "${BLUE}──────────────────────────────────────────────────────────────────────${NC}"

    # 检查 aria2c 或 wget 日志
    if [ -f "$DOWNLOAD_DIR/download_aria2.log" ]; then
        LOG="$DOWNLOAD_DIR/download_aria2.log"
    elif [ -f "$DOWNLOAD_DIR/download.log" ]; then
        LOG="$DOWNLOAD_DIR/download.log"
    else
        LOG=""
    fi

    if [ -n "$LOG" ]; then
        echo -e "  ${CYAN}日志: $(basename "$LOG")${NC}"
        echo ""
        tail -"$N" "$LOG" | while IFS= read -r line; do
            echo -e "  $line"
        done
    else
        echo -e "  ${YELLOW}暂无日志${NC}"
    fi
    echo ""
}

# ==============================================================================
# 主逻辑
# ==============================================================================

WATCH_MODE=false
INTERVAL=10
LOG_LINES=5

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
    while true; do
        clear
        print_header
        show_progress
        show_filesize
        show_disk
        show_log "$LOG_LINES"
        echo -e "${CYAN}每 ${INTERVAL} 秒刷新 | 按 Ctrl+C 退出${NC}"
        sleep "$INTERVAL"
    done
else
    print_header
    show_progress
    show_filesize
    show_disk
    show_log "$LOG_LINES"
fi
