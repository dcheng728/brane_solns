#!/usr/bin/env python3
"""
process_paper.py — Two-agent literature summary generator.

Reads an arXiv source directory from `literature_tex/<paper_dir>/`, produces a
single summary file at `agents/outputs/literature_review/<paper_dir>.md` by
alternating two Claude Code agents on the same file:

    1. Clarity Agent  — writes/revises the summary for a reader who already
                        knows the field (see agents/tasks/clarity_agent.md)
    2. Fidelity Agent — audits the summary against the source paper for
                        accuracy (see agents/tasks/fidelity_agent.md)

Usage:
    python agents/process_paper.py <paper_dir_name> [num_iterations]

Example:
    python agents/process_paper.py arXiv-XXXX.XXXXXvN
    python agents/process_paper.py arXiv-XXXX.XXXXXvN 5

Behavior:
    - If the output .md does not exist, the clarity agent writes the initial
      draft. Then the fidelity ↔ clarity loop runs up to <num_iterations>
      times against that file.
    - If the output .md already exists, the initial draft is skipped and the
      revision loop runs up to <num_iterations> times against the existing
      file. Re-running with more iterations therefore deepens an existing
      summary.
    - Each iteration overwrites the same .md file. There are no version files.
    - Convergence: each agent may print "STATUS: CONVERGED" on its own line
      when it has no material changes to make. Two consecutive convergence
      signals exit the loop early — so num_iterations is a CAP, not a target.
    - Per-call rate-limit retries are handled automatically (50 attempts,
      5 min between retries).
"""

import argparse
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
LITERATURE_DIR = REPO_ROOT / "literature_tex"
AGENTS_DIR = REPO_ROOT / "agents"
KNOWLEDGE_DIR = REPO_ROOT / "writeups"
OUTPUT_DIR = AGENTS_DIR / "outputs" / "literature_review"
LOG_DIR = AGENTS_DIR / "outputs" / "log"

CLARITY_PROMPT_FILE = AGENTS_DIR / "tasks" / "clarity_agent.md"
FIDELITY_PROMPT_FILE = AGENTS_DIR / "tasks" / "fidelity_agent.md"

MAX_RETRIES = 50
RETRY_WAIT = 300  # seconds between retries on agent failure
CONVERGED_TOKEN = "STATUS: CONVERGED"


def now() -> str:
    return datetime.now().strftime("%a %b %d %H:%M:%S %Y")


class Logger:
    """Tee-style logger: writes each line to stdout and to a log file."""

    def __init__(self, path: Path, append: bool):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.fh = path.open("a" if append else "w")

    def write(self, line: str = "") -> None:
        print(line)
        self.fh.write(line + "\n")
        self.fh.flush()

    def close(self) -> None:
        self.fh.close()


def run_agent(
    log: Logger,
    agent_name: str,
    prompt: str,
    output_file: Path,
):
    """Invoke `claude -p <prompt>`, retrying on failure.

    Returns (success, converged):
      - success:   the call achieved a valid end state (file rewritten, OR file
                   left unchanged AND the agent declared convergence).
      - converged: the agent printed the CONVERGED_TOKEN on its own line.

    Failure modes (success=False) include: file missing entirely, or file
    unchanged with no convergence signal — both treated as agent errors and
    retried up to MAX_RETRIES times.
    """
    log.write("")
    log.write(f"── {agent_name} ── {now()}")
    log.write("-------------------------------------------")

    mtime_before = output_file.stat().st_mtime if output_file.is_file() else 0.0

    for attempt in range(1, MAX_RETRIES + 1):
        proc = subprocess.Popen(
            [
                "claude",
                "--model", "opus",
                "--dangerously-skip-permissions",
                "-p", prompt,
            ],
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        converged = False
        for line in proc.stdout:
            text = line.rstrip("\n")
            log.write(text)
            if text.strip() == CONVERGED_TOKEN:
                converged = True
        exit_code = proc.wait()

        if output_file.is_file():
            mtime_after = output_file.stat().st_mtime
            if mtime_after > mtime_before:
                log.write("")
                log.write(f"── {agent_name} complete ──")
                return True, converged
            if converged:
                # File unchanged, but agent intentionally declared convergence.
                log.write("")
                log.write(f"── {agent_name} converged (no changes) ──")
                return True, True

        log.write("")
        log.write(
            f"!! {agent_name} failed "
            f"(attempt {attempt}/{MAX_RETRIES}, exit code {exit_code})"
        )
        if attempt < MAX_RETRIES:
            log.write(f"!! Waiting {RETRY_WAIT}s before retry... ({now()})")
            time.sleep(RETRY_WAIT)

    log.write(f"!! {agent_name} FAILED after {MAX_RETRIES} attempts")
    return False, False


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Two-agent literature summary generator.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "paper_name",
        help="Directory name under literature_tex/ (e.g. arXiv-XXXX.XXXXXvN)",
    )
    parser.add_argument(
        "num_iterations",
        nargs="?",
        type=int,
        default=3,
        help=(
            "Cap on fidelity<->clarity rounds (default: 3). "
            "The loop exits early after two consecutive no-op passes "
            "(see convergence behavior in module docstring)."
        ),
    )
    parser.add_argument(
        "--max-chars",
        type=int,
        default=None,
        help=(
            "Hard cap on the summary length, in characters. "
            "Use shorter for tangential papers, longer for ones worth deep engagement. "
            "Default: leave to the agent's judgment (typically 800-2000 words / ~5000-13000 chars)."
        ),
    )
    args = parser.parse_args()

    paper_dir = LITERATURE_DIR / args.paper_name
    output_file = OUTPUT_DIR / f"{args.paper_name}.md"
    log_file = LOG_DIR / f"{args.paper_name}.log"

    if not paper_dir.is_dir():
        print(f"Error: Paper directory not found: {paper_dir}", file=sys.stderr)
        return 1

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    resuming = output_file.is_file()
    log = Logger(log_file, append=resuming)

    log.write("============================================")
    log.write(f"{'RESUMING' if resuming else 'Processing'}: {args.paper_name}")
    log.write(f"Source dir: {paper_dir}")
    log.write(f"Output:     {output_file}")
    log.write(f"Iterations: {args.num_iterations}")
    log.write(f"Max chars:  {args.max_chars if args.max_chars else 'unbounded'}")
    log.write(f"{'Resumed' if resuming else 'Started'}:    {now()}")
    log.write("============================================")

    clarity_prompt_text = CLARITY_PROMPT_FILE.read_text()
    fidelity_prompt_text = FIDELITY_PROMPT_FILE.read_text()

    source_instructions = (
        f"\n**Source paper directory**: `{paper_dir}/` (an unzipped arXiv source dump). "
        "List its contents to see what's there. arXiv submissions vary: there may be a single "
        "`main.tex`, a file named after the paper, multiple `.tex` files chained via `\\input`/"
        "`\\include`/`\\subfile`, a `.ltx` extension, or a `Makefile`/`latexmk` config indicating the "
        "build entry point. Identify the main TeX entry point (typically the file containing "
        "`\\documentclass` and `\\begin{document}`), then read it AND any files it pulls in. "
        "If a `README` or build script is present, consult it. Read enough of the source to write a "
        "faithful summary — do not stop at just the first file you find."
    )

    length_instruction = (
        f"\n**Length budget**: The final summary must be at most {args.max_chars} characters "
        "(measured by `len(file_contents)`). Stay under the cap; if you need to add accuracy fixes "
        "or new material, compress or trim less essential content to make room."
        if args.max_chars
        else ""
    )

    convergence_instruction = (
        f"\n**Convergence**: This is one pass in an iterative loop. If, after reading the source "
        "and the current draft, you have NO material changes to make — no factual issues to fix and "
        "no clarity/calibration improvements worth making — print `"
        f"{CONVERGED_TOKEN}` on its own line at the end of your output and SKIP the Write step. "
        "A genuine no-op pass is preferred over manufactured edits. If material issues exist, just "
        "produce the revision (do NOT print this token)."
    )

    if not resuming:
        init_prompt = (
            f"{clarity_prompt_text}\n"
            "\n---\n"
            "\n## Instructions for this run\n"
            "\nYou are creating the INITIAL summary draft. There is no existing summary yet.\n"
            f"{source_instructions}"
            f"\n**Reader's knowledge base**: Skim files in `{KNOWLEDGE_DIR}/` to gauge what the reader already knows "
            f"(especially `{KNOWLEDGE_DIR}/12d_review_v3/`)."
            f"\n**Output**: Write the summary to `{output_file}`."
            f"{length_instruction}\n"
            "\nRead the source paper, then write the summary following the structure in your "
            "instructions above. Write it directly to the output file using the Write tool."
        )

        ok, _ = run_agent(log, "Clarity Agent (initial draft)", init_prompt, output_file)
        if not ok or not output_file.is_file():
            log.write(f"Error: Initial summary was not created at {output_file}")
            log.close()
            return 1

    consecutive_converged = 0
    converged_early = False
    rounds_completed = 0

    for i in range(1, args.num_iterations + 1):
        fidelity_run_prompt = (
            f"{fidelity_prompt_text}\n"
            "\n---\n"
            "\n## Instructions for this run\n"
            f"{source_instructions}"
            f"\n**Current summary draft**: Read the file at `{output_file}`."
            f"\n**Output**: Overwrite the summary at `{output_file}` with your revised version."
            f"{length_instruction}"
            f"{convergence_instruction}\n"
            "\nRead the source paper and the current draft. Audit the summary against the paper. "
            "Print your critique to stdout. If material issues exist, write the revised summary to the "
            "output file using the Write tool (overwriting the previous version). Otherwise emit the "
            "convergence token and skip the Write step."
        )
        _, fid_converged = run_agent(
            log, f"Fidelity Agent (round {i}/{args.num_iterations})", fidelity_run_prompt, output_file
        )
        consecutive_converged = consecutive_converged + 1 if fid_converged else 0
        if consecutive_converged >= 2:
            log.write("")
            log.write(f"== Converged: 2 consecutive no-op passes after fidelity round {i}. ==")
            converged_early = True
            rounds_completed = i
            break

        clarity_run_prompt = (
            f"{clarity_prompt_text}\n"
            "\n---\n"
            "\n## Instructions for this run\n"
            "\nYou are REVISING an existing summary.\n"
            f"{source_instructions}"
            f"\n**Current summary draft**: Read the file at `{output_file}`."
            f"\n**Reader's knowledge base**: Skim files in `{KNOWLEDGE_DIR}/` if you need to recalibrate."
            f"\n**Output**: Overwrite the summary at `{output_file}` with your revised version."
            f"{length_instruction}"
            f"{convergence_instruction}\n"
            "\nRead the current draft. Audit it for clarity. Print your critique to stdout. If material "
            "issues exist, write the revised summary to the output file using the Write tool (overwriting "
            "the previous version). Otherwise emit the convergence token and skip the Write step."
        )
        _, cla_converged = run_agent(
            log, f"Clarity Agent (round {i}/{args.num_iterations})", clarity_run_prompt, output_file
        )
        consecutive_converged = consecutive_converged + 1 if cla_converged else 0
        rounds_completed = i
        if consecutive_converged >= 2:
            log.write("")
            log.write(f"== Converged: 2 consecutive no-op passes after clarity round {i}. ==")
            converged_early = True
            break

    log.write("")
    log.write("============================================")
    if converged_early:
        log.write(f"Done (converged after {rounds_completed} round(s)): {output_file}")
    else:
        log.write(f"Done (hit cap of {args.num_iterations} rounds): {output_file}")
    log.write(f"Completed: {now()}")
    log.write("============================================")
    log.close()

    return 0


if __name__ == "__main__":
    sys.exit(main())
