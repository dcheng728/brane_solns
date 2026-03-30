# Papers to Read — Adversarial Tutorial Generator

Drop arXiv paper source directories here, then run the script to generate a tutorial.

## Usage

```bash
# From the repo root:
bash papers_2_read/process_paper.sh <paper_dir_name> [num_iterations]

# Example:
bash papers_2_read/process_paper.sh <paper_dir_name>        # 3 rounds (default)
bash papers_2_read/process_paper.sh <paper_dir_name> 5      # 5 rounds
```

## How it works

Two Claude (Opus) agents alternate in a loop:

1. **Clarity Agent** writes the initial tutorial draft, calibrated to your knowledge base in `writeups/`
2. **Fidelity Agent** audits the draft against the source paper, fixes inaccuracies and fills gaps
3. **Clarity Agent** revises for readability and notation consistency
4. Repeat steps 2-3 for N iterations

## Directory structure

```
papers_2_read/
├── agents/
│   ├── fidelity_agent.md    # Accuracy-checking agent prompt
│   └── clarity_agent.md     # Readability-checking agent prompt
├── process_paper.sh         # Orchestrator script
├── readme.md
└── arXiv-XXXX.XXXXXvN/      # Paper directory (from arXiv source download)
    ├── main.tex              # Source paper (required)
    ├── tutorial.md           # Final tutorial (generated)
    ├── tutorial_v0.md        # Initial draft
    ├── tutorial_v1.md        # After fidelity pass 1
    ├── tutorial_v2.md        # After clarity pass 1
    ├── ...                   # More versions
    └── agent_log.txt         # Full agent output log
```

## Requirements

- `claude` CLI installed and authenticated
- Paper directory must contain `main.tex`
