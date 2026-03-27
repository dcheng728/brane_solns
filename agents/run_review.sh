#!/bin/bash
# ─── Overnight review agents for 12d_review_v3 ──────────────────────────────
#
# Each agent is defined in agents/tasks/<name>.md with YAML frontmatter:
#   model:  opus | sonnet | haiku
#   tools:  comma-separated allowed tools
#
# Usage:
#   ./agents/run_review.sh                  # run all tasks
#   ./agents/run_review.sh directions       # run one task by name
#   ./agents/run_review.sh task1 task2      # run specific tasks
#
# Results land in agents/results/<name>_<timestamp>.md
# Retries on failure (rate limit, usage cap) every 30 min.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
TASKS_DIR="$REPO_ROOT/agents/tasks"
RESULTS="$REPO_ROOT/agents/results"
TIMESTAMP=$(date +%Y%m%d_%H%M%S)

mkdir -p "$RESULTS"

# ─── Parse frontmatter from a task file ──────────────────────────────────────

parse_field() {
    # $1 = file, $2 = field name
    sed -n '/^---$/,/^---$/p' "$1" | grep "^$2:" | sed "s/^$2:[[:space:]]*//"
}

# ─── Run a single agent ─────────────────────────────────────────────────────

run_agent() {
    local task_file="$1"
    local name
    name="$(basename "$task_file" .md)"

    local model tools prompt
    model="$(parse_field "$task_file" model)"
    tools="$(parse_field "$task_file" tools)"

    # Prompt = everything after the closing ---
    prompt="$(sed '1,/^---$/{ /^---$/!d; /^---$/{ x; s/.*//; x; d; }; }; 1,/^---$/d' "$task_file")"

    # Map model names to claude flags
    local model_flag=""
    case "${model:-sonnet}" in
        opus)   model_flag="--model claude-opus-4-6" ;;
        sonnet) model_flag="--model claude-sonnet-4-6" ;;
        haiku)  model_flag="--model claude-haiku-4-5" ;;
        *)      model_flag="--model $model" ;;
    esac

    local outfile="$RESULTS/${name}_${TIMESTAMP}.md"
    local logfile="$RESULTS/${name}_${TIMESTAMP}.log"

    echo "[$(date '+%H:%M:%S')] Starting agent: $name (model: ${model:-sonnet}, tools: ${tools:-all})"

    while true; do
        claude -p "$prompt" \
            $model_flag \
            ${tools:+--allowedTools "$tools"} \
            --output-file "$outfile" \
            2>"$logfile"

        if [ $? -eq 0 ]; then
            echo "[$(date '+%H:%M:%S')] Agent '$name' finished. Output: $outfile"
            break
        else
            echo "[$(date '+%H:%M:%S')] Agent '$name' failed. Retrying in 30 min..."
            sleep 1800
        fi
    done
}

# ─── Main ────────────────────────────────────────────────────────────────────

# Collect task files
task_files=()
if [ $# -gt 0 ]; then
    for name in "$@"; do
        f="$TASKS_DIR/${name}.md"
        if [ ! -f "$f" ]; then
            echo "Error: task '$name' not found at $f" >&2
            exit 1
        fi
        task_files+=("$f")
    done
else
    for f in "$TASKS_DIR"/*.md; do
        [ -f "$f" ] && task_files+=("$f")
    done
fi

if [ ${#task_files[@]} -eq 0 ]; then
    echo "No task files found in $TASKS_DIR/" >&2
    exit 1
fi

# Launch all agents in parallel
for f in "${task_files[@]}"; do
    run_agent "$f" &
done

echo ""
echo "Launched ${#task_files[@]} agents. Waiting for completion..."
echo "Results will appear in: $RESULTS/"
echo ""
echo "To monitor: tail -f $RESULTS/*_${TIMESTAMP}.log"
echo "To stop all: kill 0"
echo ""

wait
echo ""
echo "═══════════════════════════════════════════════"
echo " All agents finished at $(date '+%H:%M:%S')"
echo " Results in: $RESULTS/"
echo "═══════════════════════════════════════════════"
ls -lh "$RESULTS"/*_${TIMESTAMP}.md 2>/dev/null
