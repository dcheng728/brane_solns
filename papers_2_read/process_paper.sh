#!/bin/bash
# ============================================================================
# process_paper.sh — Two-agent adversarial tutorial generator
#
# Usage:
#   bash papers_2_read/process_paper.sh <paper_dir_name> [num_iterations]
#
# Example:
#   bash papers_2_read/process_paper.sh <paper_dir_name>
#   bash papers_2_read/process_paper.sh <paper_dir_name> 5
#
# This script runs two Claude Code agents in a loop:
#   1. Clarity Agent  — writes/revises the tutorial for readability
#   2. Fidelity Agent — audits against the source paper for accuracy
# They alternate, each refining the tutorial, for N iterations.
#
# Rate limit handling:
#   If claude hits a rate limit, the script waits and retries automatically.
#   All progress is saved — each version is its own file (tutorial_v0.md, etc).
#   If you kill the script, just re-run it: it resumes from the last completed version.
# ============================================================================

set -uo pipefail  # no -e: we handle errors manually for retry logic

# ── Configuration ──────────────────────────────────────────────────────────
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PAPERS_DIR="${REPO_ROOT}/papers_2_read"
AGENTS_DIR="${PAPERS_DIR}/agents"
KNOWLEDGE_DIR="${REPO_ROOT}/writeups"

CLARITY_PROMPT="${AGENTS_DIR}/clarity_agent.md"
FIDELITY_PROMPT="${AGENTS_DIR}/fidelity_agent.md"

NUM_ITERATIONS="${2:-3}"  # default 3 rounds of fidelity->clarity
MAX_RETRIES=50            # max retries per agent call (for rate limits)
RETRY_WAIT=300            # seconds to wait between retries (5 min)

# ── Validate input ─────────────────────────────────────────────────────────
if [ $# -lt 1 ]; then
    echo "Usage: $0 <paper_dir_name> [num_iterations]"
    echo "  paper_dir_name: directory name under papers_2_read/ (e.g. arXiv-XXXX.XXXXXvN)"
    echo "  num_iterations: number of fidelity<->clarity rounds (default: 3)"
    exit 1
fi

PAPER_NAME="$1"
PAPER_DIR="${PAPERS_DIR}/${PAPER_NAME}"
SOURCE_FILE="${PAPER_DIR}/main.tex"
LOG_FILE="${PAPER_DIR}/agent_log.txt"

if [ ! -d "$PAPER_DIR" ]; then
    echo "Error: Paper directory not found: ${PAPER_DIR}"
    exit 1
fi

if [ ! -f "$SOURCE_FILE" ]; then
    echo "Error: Source file not found: ${SOURCE_FILE}"
    echo "Expected main.tex in ${PAPER_DIR}/"
    exit 1
fi

# ── Find resume point ─────────────────────────────────────────────────────
# Check which versions already exist so we can skip completed steps
find_latest_version() {
    local latest=-1
    for f in "${PAPER_DIR}"/tutorial_v*.md; do
        [ -f "$f" ] || continue
        local v="${f##*tutorial_v}"
        v="${v%.md}"
        if [ "$v" -gt "$latest" ] 2>/dev/null; then
            latest="$v"
        fi
    done
    echo "$latest"
}

RESUME_FROM=$(find_latest_version)

# ── Setup ──────────────────────────────────────────────────────────────────
if [ "$RESUME_FROM" -ge 0 ] 2>/dev/null; then
    echo "============================================" | tee -a "$LOG_FILE"
    echo "RESUMING: ${PAPER_NAME} from v${RESUME_FROM}" | tee -a "$LOG_FILE"
    echo "Resumed:  $(date)" | tee -a "$LOG_FILE"
    echo "============================================" | tee -a "$LOG_FILE"
else
    echo "============================================" | tee "$LOG_FILE"
    echo "Processing: ${PAPER_NAME}" | tee -a "$LOG_FILE"
    echo "Source:     ${SOURCE_FILE}" | tee -a "$LOG_FILE"
    echo "Output:     ${PAPER_DIR}/tutorial.md" | tee -a "$LOG_FILE"
    echo "Iterations: ${NUM_ITERATIONS}" | tee -a "$LOG_FILE"
    echo "Started:    $(date)" | tee -a "$LOG_FILE"
    echo "============================================" | tee -a "$LOG_FILE"
fi

# ── Helper: run a claude agent with retry ──────────────────────────────────
run_agent() {
    local agent_name="$1"
    local prompt="$2"
    local version="$3"
    local output_file="$4"

    echo "" | tee -a "$LOG_FILE"
    echo "── ${agent_name} (v${version}) ── $(date)" | tee -a "$LOG_FILE"
    echo "-------------------------------------------" | tee -a "$LOG_FILE"

    local attempt=0
    while [ $attempt -lt $MAX_RETRIES ]; do
        # Run claude with opus model
        claude --model opus --dangerously-skip-permissions -p "$prompt" 2>&1 | tee -a "$LOG_FILE"
        local exit_code=${PIPESTATUS[0]}

        # Check if the output file was created (success)
        if [ -f "$output_file" ]; then
            echo "" | tee -a "$LOG_FILE"
            echo "── ${agent_name} (v${version}) complete ──" | tee -a "$LOG_FILE"
            return 0
        fi

        # Agent failed — likely rate limit or error
        attempt=$((attempt + 1))
        echo "" | tee -a "$LOG_FILE"
        echo "!! ${agent_name} (v${version}) failed (attempt ${attempt}/${MAX_RETRIES}, exit code ${exit_code})" | tee -a "$LOG_FILE"

        if [ $attempt -lt $MAX_RETRIES ]; then
            echo "!! Waiting ${RETRY_WAIT}s before retry... ($(date))" | tee -a "$LOG_FILE"
            sleep $RETRY_WAIT
        fi
    done

    echo "!! ${agent_name} (v${version}) FAILED after ${MAX_RETRIES} attempts" | tee -a "$LOG_FILE"
    return 1
}

# ── Iteration 0: Initial draft (Clarity Agent) ────────────────────────────
CURRENT_TUTORIAL="${PAPER_DIR}/tutorial_v0.md"

if [ "$RESUME_FROM" -lt 0 ] 2>/dev/null; then
    INIT_PROMPT="$(cat "$CLARITY_PROMPT")

---

## Instructions for this run

You are creating the INITIAL tutorial draft. There is no existing tutorial yet.

**Source paper**: Read the file at \`${SOURCE_FILE}\`
**Reader's knowledge base**: Skim files in \`${KNOWLEDGE_DIR}/\` to calibrate your level (especially \`${KNOWLEDGE_DIR}/12d_review_v3/\` and \`${KNOWLEDGE_DIR}/scratch/\`)
**Output**: Write the tutorial to \`${CURRENT_TUTORIAL}\`

Read the source paper fully, then write the tutorial following the structure in your instructions above. Write it directly to the output file using the Write tool."

    run_agent "Clarity Agent (initial draft)" "$INIT_PROMPT" 0 "$CURRENT_TUTORIAL"

    if [ ! -f "$CURRENT_TUTORIAL" ]; then
        echo "Error: Initial tutorial was not created at ${CURRENT_TUTORIAL}" | tee -a "$LOG_FILE"
        exit 1
    fi
else
    echo "Skipping v0 (already exists)" | tee -a "$LOG_FILE"
fi

# ── Main loop: alternate fidelity and clarity ──────────────────────────────
for i in $(seq 1 "$NUM_ITERATIONS"); do
    FIDELITY_VERSION=$((2 * i - 1))
    CLARITY_VERSION=$((2 * i))

    # Determine the previous tutorial for this step
    PREV_FIDELITY=$((FIDELITY_VERSION - 1))
    PREV_TUTORIAL="${PAPER_DIR}/tutorial_v${PREV_FIDELITY}.md"

    # ── Fidelity Agent pass ──
    CURRENT_TUTORIAL="${PAPER_DIR}/tutorial_v${FIDELITY_VERSION}.md"

    if [ "$RESUME_FROM" -lt "$FIDELITY_VERSION" ] 2>/dev/null || [ "$RESUME_FROM" -lt 0 ] 2>/dev/null; then
        FIDELITY_RUN_PROMPT="$(cat "$FIDELITY_PROMPT")

---

## Instructions for this run

**Source paper**: Read the file at \`${SOURCE_FILE}\`
**Current tutorial draft**: Read the file at \`${PREV_TUTORIAL}\`
**Output**: Write the revised tutorial to \`${CURRENT_TUTORIAL}\`

Read both files. Audit the tutorial against the source paper. Print your critique, then write the revised tutorial to the output file using the Write tool."

        run_agent "Fidelity Agent" "$FIDELITY_RUN_PROMPT" "$FIDELITY_VERSION" "$CURRENT_TUTORIAL"

        if [ ! -f "$CURRENT_TUTORIAL" ]; then
            echo "Warning: Fidelity agent did not produce v${FIDELITY_VERSION}, using previous version" | tee -a "$LOG_FILE"
            cp "$PREV_TUTORIAL" "$CURRENT_TUTORIAL"
        fi
    else
        echo "Skipping v${FIDELITY_VERSION} (already exists)" | tee -a "$LOG_FILE"
    fi

    # ── Clarity Agent pass ──
    PREV_TUTORIAL="$CURRENT_TUTORIAL"
    CURRENT_TUTORIAL="${PAPER_DIR}/tutorial_v${CLARITY_VERSION}.md"

    if [ "$RESUME_FROM" -lt "$CLARITY_VERSION" ] 2>/dev/null || [ "$RESUME_FROM" -lt 0 ] 2>/dev/null; then
        CLARITY_RUN_PROMPT="$(cat "$CLARITY_PROMPT")

---

## Instructions for this run

You are REVISING an existing tutorial.

**Source paper**: Read the file at \`${SOURCE_FILE}\` (for reference)
**Current tutorial draft**: Read the file at \`${PREV_TUTORIAL}\`
**Reader's knowledge base**: Skim files in \`${KNOWLEDGE_DIR}/\` if you need to recalibrate
**Output**: Write the revised tutorial to \`${CURRENT_TUTORIAL}\`

Read the current draft. Audit it for clarity. Print your critique, then write the revised tutorial to the output file using the Write tool."

        run_agent "Clarity Agent" "$CLARITY_RUN_PROMPT" "$CLARITY_VERSION" "$CURRENT_TUTORIAL"

        if [ ! -f "$CURRENT_TUTORIAL" ]; then
            echo "Warning: Clarity agent did not produce v${CLARITY_VERSION}, using previous version" | tee -a "$LOG_FILE"
            cp "$PREV_TUTORIAL" "$CURRENT_TUTORIAL"
        fi
    else
        echo "Skipping v${CLARITY_VERSION} (already exists)" | tee -a "$LOG_FILE"
    fi
done

# ── Finalize ───────────────────────────────────────────────────────────────
cp "$CURRENT_TUTORIAL" "${PAPER_DIR}/tutorial.md"

echo "" | tee -a "$LOG_FILE"
echo "============================================" | tee -a "$LOG_FILE"
echo "Done! Final tutorial: ${PAPER_DIR}/tutorial.md" | tee -a "$LOG_FILE"
echo "Completed: $(date)" | tee -a "$LOG_FILE"
echo "============================================" | tee -a "$LOG_FILE"
echo ""
echo "All versions saved:"
ls -la "${PAPER_DIR}"/tutorial_v*.md "${PAPER_DIR}"/tutorial.md 2>/dev/null
