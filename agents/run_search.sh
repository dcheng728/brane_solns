#!/bin/bash
# Re-launches the 12d form field search agent until it finds a solution.
# The agent reads/writes agents/notes.md as persistent memory between runs.

cd "$(dirname "$0")/.."

echo "Starting 12d form field search loop..."
echo "Agent notes: agents/notes.md"
echo "Kill with Ctrl+C"
echo ""

RUN=1
while ! grep -q "SOLUTION FOUND" agents/notes.md 2>/dev/null; do
    echo "========================================"
    echo "  Run #$RUN"
    echo "========================================"
    claude --agent find-12d-form -p "Read agents/notes.md and continue searching from where you left off. Do not repeat experiments already recorded in the notes."
    RUN=$((RUN + 1))
    echo ""
    echo "Agent session ended. Checking for solution..."
    sleep 2
done

echo ""
echo "========================================"
echo "  SOLUTION FOUND! Check agents/notes.md"
echo "========================================"
