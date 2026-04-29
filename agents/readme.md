# Agents

Agent prompts and orchestrator scripts for the project.

## Literature summary loop

Generates project-relative summaries of arXiv papers, written for someone who already knows what's in `writeups/`. Two Claude Code agents iterate on a single `.md` file:

- **Clarity Agent** ([tasks/clarity_agent.md](tasks/clarity_agent.md)) — writes / revises for the reader's perspective; surfaces relevance to the 12d/IIB project goal.
- **Fidelity Agent** ([tasks/fidelity_agent.md](tasks/fidelity_agent.md)) — audits the summary against the source paper for accuracy.

### Workflow

1. On the arXiv abstract page, click **"Other formats" → "Source"** to download the paper's source as a tarball.
2. Unzip it into `literature_tex/<paper_dir>/` at the repo root. The directory should contain `main.tex`. (`literature_tex/` is gitignored.)
3. Run the orchestrator:
   ```bash
   python agents/process_paper.py <paper_dir> [num_iterations] [--max-chars N]
   ```
   For example:
   ```bash
   python agents/process_paper.py arXiv-2401.12345v1                       # 3 rounds (default)
   python agents/process_paper.py arXiv-2401.12345v1 5                     # 5 rounds
   python agents/process_paper.py arXiv-2401.12345v1 --max-chars 4000      # short summary (~600 words)
   python agents/process_paper.py arXiv-2401.12345v1 --max-chars 15000     # long summary (~2300 words)
   ```
   Use `--max-chars` to size the summary to the paper: tight cap for tangential papers, generous cap for ones worth deep engagement. The agent doesn't need to know which TeX file is the entry point — it'll inspect the directory and figure out the structure (multi-file `\input` chains, weird filenames, etc.).
4. Output appears at `agents/outputs/literature_review/<paper_dir>.md`. The full agent transcript is logged to `agents/outputs/log/<paper_dir>.log` (logs are gitignored via `*.log`).

Re-running on the same paper continues to refine the existing summary file in place — there are no `_v0`, `_v1` versions. To start over, delete the `.md` and re-run.

`num_iterations` is a **cap**, not a target: each agent can declare convergence (`STATUS: CONVERGED`) when it has no material changes to make, and the loop exits early after two consecutive no-op passes. The trailing log line tells you whether the run converged or hit the cap.

## Other agent scripts

- `run_review.sh` / `run_search.sh` / `run_tensor_simplify.py` / `run_tensorglow_tests.sh` — see `notes.md`.
