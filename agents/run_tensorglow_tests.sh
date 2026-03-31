#!/usr/bin/env bash
#
# Stress-test tensorGlow with two Claude Code agents.
#
# Each cycle:
#   1. Physicist designs a test (may refer to psets/ for inspiration)
#   2. Engineer implements it with tensorGlow (no hard-coding)
#   3. Physicist checks the result and writes a report
#
# Output:
#   agents/tensorglow_report.md  — structured findings
#   agents/tensorglow_session.log — full transcript
#
# Usage:
#   bash agents/run_tensorglow_tests.sh [N_CYCLES]

set -euo pipefail

N_CYCLES=${1:-20}
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
REPORT="$REPO_ROOT/agents/tensorglow_report.md"
LOG="$REPO_ROOT/agents/tensorglow_session.log"

PHYSICIST_MD="$REPO_ROOT/agents/tasks/tensorglow_physicist.md"
ENGINEER_MD="$REPO_ROOT/agents/tasks/tensorglow_engineer.md"

cd "$REPO_ROOT"

# Extract system prompt from .md (strip YAML frontmatter)
get_system_prompt() {
    sed '1,/^---$/d; 1,/^---$/d' "$1"
}

PHYSICIST_SYSTEM=$(get_system_prompt "$PHYSICIST_MD")
ENGINEER_SYSTEM=$(get_system_prompt "$ENGINEER_MD")

# ─── Initialize ──────────────────────────────────────────────────

: > "$LOG"
PASS_ITEMS=""
BUG_ITEMS=""
LIMIT_ITEMS=""
PASS=0; BUG=0; CODE_ERROR=0; LIMITATION=0

log() { echo "$1" | tee -a "$LOG"; }

log "================================================================"
log "  tensorGlow stress-test — $(date)"
log "  $N_CYCLES cycles"
log "================================================================"

# ─── Main loop ───────────────────────────────────────────────────

HISTORY=""

for i in $(seq 1 "$N_CYCLES"); do
    log ""
    log "──── CYCLE $i / $N_CYCLES ────"

    # ── Phase 1: Physicist designs a test ──

    POSE_PROMPT="Design test problem $i. You may read the GR problem sheets in src/tensorGlow/psets/ and the examples in src/tensorGlow/examples/ for inspiration. Follow the format in your instructions (PROBLEM / EXPECTED RESULT / VERIFICATION)."
    if [ -n "$HISTORY" ]; then
        POSE_PROMPT="$POSE_PROMPT

Previous results: $HISTORY

Pick something DIFFERENT from what has already been tested. If bugs were found, probe that area deeper."
    fi

    log "[PHYSICIST] designing..."
    PROBLEM=$(claude -p "$POSE_PROMPT" \
        --system-prompt "$PHYSICIST_SYSTEM" \
        --model sonnet \
        --allowedTools "Read,Glob,Grep" \
        2>>"$LOG") || { log "  ✗ physicist failed"; continue; }
    log "$PROBLEM"

    # ── Phase 2: Engineer implements and runs ──

    log ""
    log "[ENGINEER] implementing..."
    ENGINEER_OUTPUT=$(claude -p "Implement this test. Write the script, run it, report the full output.

$PROBLEM" \
        --system-prompt "$ENGINEER_SYSTEM" \
        --model sonnet \
        --allowedTools "Bash,Write,Read,Glob,Grep" \
        2>>"$LOG") || { log "  ✗ engineer failed"; continue; }
    log "$ENGINEER_OUTPUT"

    # ── Phase 3: Physicist reviews and reports ──

    log ""
    log "[PHYSICIST] reviewing..."
    REVIEW=$(claude -p "Review this tensorGlow test. Follow your Phase 2 instructions.

PROBLEM:
$PROBLEM

ENGINEER REPORT:
$ENGINEER_OUTPUT" \
        --system-prompt "$PHYSICIST_SYSTEM" \
        --model sonnet \
        --allowedTools "Read,Glob,Grep" \
        2>>"$LOG") || { log "  ✗ review failed"; continue; }
    log "$REVIEW"

    # ── Parse result ──

    CLASS=$(echo "$REVIEW" | grep -oi "CLASSIFICATION: *[A-Z_]*" | head -1 | sed 's/CLASSIFICATION:[[:space:]]*//' || echo "UNKNOWN")
    SEVERITY=$(echo "$REVIEW" | grep -oi "SEVERITY: *[a-z]*" | head -1 | sed 's/SEVERITY:[[:space:]]*//' || echo "n/a")
    ONELINE=$(echo "$PROBLEM" | grep -i "^PROBLEM:" | head -1 | sed 's/PROBLEM:[[:space:]]*//' | cut -c1-80 || echo "see log")

    case "$CLASS" in
        PASS)       ((PASS++)) || true
                    PASS_ITEMS="$PASS_ITEMS
- **#$i**: $ONELINE" ;;
        BUG)        ((BUG++)) || true
                    BUG_ITEMS="$BUG_ITEMS
- **#$i** ($SEVERITY): $ONELINE" ;;
        LIMITATION) ((LIMITATION++)) || true
                    LIMIT_ITEMS="$LIMIT_ITEMS
- **#$i**: $ONELINE" ;;
        CODE_ERROR) ((CODE_ERROR++)) || true ;;
    esac

    log "  → $CLASS"
    HISTORY="$HISTORY #$i($CLASS)"

    rm -f src/tensorGlow/examples/_test_scratch.py
done

# ─── Write report ────────────────────────────────────────────────

cat > "$REPORT" << EOF
# tensorGlow Test Report

Auto-generated: $(date)
Cycles: $N_CYCLES

## Confirmed working (no bugs)
$PASS_ITEMS

## Potential bugs
$BUG_ITEMS

## Missing features / limitations
$LIMIT_ITEMS

## Summary

| Result | Count |
|--------|-------|
| PASS | $PASS |
| BUG | $BUG |
| CODE_ERROR | $CODE_ERROR |
| LIMITATION | $LIMITATION |

Full session log: \`agents/tensorglow_session.log\`
EOF

log ""
log "================================================================"
log "  DONE — PASS:$PASS  BUG:$BUG  CODE_ERROR:$CODE_ERROR  LIMITATION:$LIMITATION"
log "  Report: $REPORT"
log "================================================================"
