"""
Tensor simplification loop using two agents: proposer + verifier.

Usage:
    python agents/run_tensor_simplify.py
    python agents/run_tensor_simplify.py --max-rounds 30
    python agents/run_tensor_simplify.py --proposer-model claude-sonnet-4-6

The working document is:
    writeups/eff_actions/simplifying_t8t8RRDPbDP.md

Output (appended to working document + separate log):
    writeups/eff_actions/simplifying_t8t8RRDPbDP.md  — derivation steps appended
    agents/results/tensor_simplify_<timestamp>.log    — full agent transcripts
"""

import asyncio
import argparse
import os
import sys
from datetime import datetime
from pathlib import Path

# Use OAuth (Max subscription), not API key
os.environ.pop("ANTHROPIC_API_KEY", None)

from claude_code_sdk import query, ClaudeCodeOptions, AssistantMessage, ResultMessage
from claude_code_sdk._errors import MessageParseError


# ═══════════════════════════════════════════════════════════════════════
# INPUT: Read from the working document
# ═══════════════════════════════════════════════════════════════════════

REPO_ROOT = Path(__file__).resolve().parent.parent
WORKING_DOC = REPO_ROOT / "writeups" / "eff_actions" / "simplifying_t8t8RRDPbDP.md"


def load_working_doc() -> str:
    """Read the full working document so agents have all prior context."""
    return WORKING_DOC.read_text(encoding="utf-8")


SYMMETRIES = r"""
### Tensor definitions
- $R_{mnrs}$: Riemann tensor of a 10d metric, with standard symmetries:
  $R_{mnrs} = -R_{nmrs} = -R_{mnsr} = R_{rsmn}$
- $K_{mn} \equiv D_m P_n$: symmetric in $(m,n)$, i.e. $K_{mn} = K_{nm}$
- $\bar{K}_{mn} \equiv D_m \bar{P}_n$: symmetric in $(m,n)$, i.e. $\bar{K}_{mn} = \bar{K}_{nm}$

### Composite object
- Define $\mathcal{K}^{mnrs} \equiv \bar{K}^{mn} K^{rs}$
- $\mathcal{K}^{mnrs}$ is symmetric in $(m,n)$ and symmetric in $(r,s)$
- For any contraction with a real tensor $T$:
  $\mathcal{K}^{mnrs} T_{mnrs} = \mathcal{K}^{rsmn} T_{mnrs}$
  (Proof: complex conjugation, since the effective action is real)

### Available identities
- First Bianchi identity: $R_{m[nrs]} = 0$, i.e.
  $R_{mnrs} + R_{mrsn} + R_{msnr} = 0$
- Riemann pair symmetry: $R_{mnrs} = R_{rsmn}$
- Contracting a symmetric tensor with an antisymmetric one gives zero

### Worked example: Bianchi application
Consider $R_{rtsu}$. Bianchi $R_{r[tsu]} = 0$ gives:
$R_{rtsu} + R_{rsut} + R_{rust} = 0$,
hence $R_{rtsu} = -R_{rsut} - R_{rust}$.

### Worked example: vanishing by symmetry
$K^{mn} R_{mnrs} = 0$ because $K^{mn}$ is symmetric in $(m,n)$
while $R_{mnrs}$ is antisymmetric in $(m,n)$.
"""


# ═══════════════════════════════════════════════════════════════════════
# AGENT PROMPTS
# ═══════════════════════════════════════════════════════════════════════

def load_task(name: str) -> str:
    """Load a task .md file, stripping YAML frontmatter."""
    path = Path(__file__).parent / "tasks" / f"{name}.md"
    text = path.read_text(encoding="utf-8")
    # Strip frontmatter (between --- markers)
    parts = text.split("---", 2)
    if len(parts) >= 3:
        return parts[2].strip()
    return text.strip()


def build_proposer_prompt(
    working_doc: str,
    symmetries: str,
    rejected: list[str] | None = None,
) -> str:
    instructions = load_task("tensor_proposer")
    prompt = f"""{instructions}

## Tensor symmetries and identities
{symmetries}

## Working document (contains the expression, prior derivations, and known identities)
{working_doc}

## Your task
Read the working document above. It contains the expression to simplify and any
identities/simplifications already established. Find ONE further simplification
of the expression that has not already been done. Apply it to the CURRENT state
of the expression (after all previously established simplifications).
"""
    if rejected:
        prompt += "\n## Previously rejected proposals (do NOT repeat these)\n"
        for i, r in enumerate(rejected, 1):
            prompt += f"\n### Rejection {i}\n{r}\n"
    return prompt


def build_verifier_prompt(
    working_doc: str,
    symmetries: str,
    proposal: str,
) -> str:
    instructions = load_task("tensor_verifier")
    return f"""{instructions}

## Tensor symmetries and identities
{symmetries}

## Working document (contains the expression and known identities)
{working_doc}

## Proposed simplification to verify
{proposal}
"""


# ═══════════════════════════════════════════════════════════════════════
# AGENT RUNNER
# ═══════════════════════════════════════════════════════════════════════

async def run_agent(prompt: str, model: str) -> str:
    """Run a Claude Code agent and collect its text output."""
    chunks: list[str] = []
    try:
        async for message in query(
            prompt=prompt,
            options=ClaudeCodeOptions(
                model=model,
                allowed_tools=[],
            ),
        ):
            if isinstance(message, AssistantMessage):
                for block in message.content:
                    if hasattr(block, "text"):
                        chunks.append(block.text)
            elif isinstance(message, ResultMessage):
                pass
    except MessageParseError:
        pass
    return "\n".join(chunks)


# ═══════════════════════════════════════════════════════════════════════
# OUTPUT
# ═══════════════════════════════════════════════════════════════════════

def extract_simplified_expression(proposal: str) -> str | None:
    """Pull out the SIMPLIFIED_EXPRESSION block from a proposal."""
    marker = "SIMPLIFIED_EXPRESSION:"
    if marker in proposal:
        after = proposal.split(marker, 1)[1]
        # Take everything until the next ``` fence or end of text
        if "```" in after:
            after = after.split("```")[0]
        return after.strip()
    return None


def is_done(proposal: str) -> bool:
    return "NO_MORE_SIMPLIFICATIONS" in proposal


def is_accepted(verdict: str) -> bool:
    # Look for VERDICT: ACCEPT (case-insensitive, flexible whitespace)
    for line in verdict.splitlines():
        stripped = line.strip().upper()
        if stripped.startswith("VERDICT:") and "ACCEPT" in stripped:
            return True
    return False


def append_step_to_working_doc(step: dict, step_num: int):
    """Append one accepted simplification step to the working document."""
    entry = [
        "",
        f"## Agent Simplification — Step {step_num}",
        f"*Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}*",
        "",
        "### Proposal",
        "",
        step["proposal"],
        "",
        "### Verification",
        "",
        step["verdict"],
        "",
        "---",
    ]
    with open(WORKING_DOC, "a", encoding="utf-8") as f:
        f.write("\n".join(entry) + "\n")


def write_log(entries: list[str], logpath: Path):
    """Append-style log of all agent interactions."""
    logpath.write_text("\n\n".join(entries), encoding="utf-8")


# ═══════════════════════════════════════════════════════════════════════
# MAIN LOOP
# ═══════════════════════════════════════════════════════════════════════

async def main():
    parser = argparse.ArgumentParser(description="Tensor simplification loop")
    parser.add_argument("--max-rounds", type=int, default=20)
    parser.add_argument("--max-retries", type=int, default=2,
                        help="Max retries per round on rejection")
    parser.add_argument("--proposer-model", default="claude-opus-4-6")
    parser.add_argument("--verifier-model", default="claude-opus-4-6")
    args = parser.parse_args()

    # Paths
    results_dir = REPO_ROOT / "agents" / "results"
    results_dir.mkdir(exist_ok=True)
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    logpath = results_dir / f"tensor_simplify_{ts}.log"

    steps: list[dict] = []
    log_entries: list[str] = []
    step_count = 0

    print(f"Starting tensor simplification loop")
    print(f"  Proposer model: {args.proposer_model}")
    print(f"  Verifier model: {args.verifier_model}")
    print(f"  Max rounds: {args.max_rounds}")
    print(f"  Working doc: {WORKING_DOC}")
    print(f"  Log: {logpath}")
    print()

    for round_num in range(1, args.max_rounds + 1):
        print(f"{'='*60}")
        print(f"  Round {round_num}")
        print(f"{'='*60}")

        # Re-read the working doc each round (it accumulates prior steps)
        working_doc = load_working_doc()

        rejected_proposals: list[str] = []
        accepted = False

        for retry in range(args.max_retries + 1):
            # ── Propose ──
            attempt = f"(attempt {retry + 1})" if retry > 0 else ""
            print(f"  [Proposer] Generating simplification {attempt}...")

            proposer_prompt = build_proposer_prompt(
                working_doc, SYMMETRIES, rejected_proposals or None
            )
            proposal = await run_agent(proposer_prompt, args.proposer_model)

            log_entries.append(
                f"# Round {round_num}, Attempt {retry + 1} — PROPOSER\n\n"
                f"## Response\n\n{proposal}"
            )

            # Check termination
            if is_done(proposal):
                print(f"  [Proposer] No more simplifications found.")
                write_log(log_entries, logpath)
                print(f"\nDone after {len(steps)} simplification(s).")
                print(f"  Working doc: {WORKING_DOC}")
                print(f"  Full log:    {logpath}")
                return

            # ── Verify ──
            print(f"  [Verifier] Checking proposal...")

            verifier_prompt = build_verifier_prompt(
                working_doc, SYMMETRIES, proposal
            )
            verdict = await run_agent(verifier_prompt, args.verifier_model)

            log_entries.append(
                f"# Round {round_num}, Attempt {retry + 1} — VERIFIER\n\n"
                f"## Response\n\n{verdict}"
            )

            if is_accepted(verdict):
                print(f"  [Verifier] ACCEPTED")

                step = {
                    "round": round_num,
                    "proposal": proposal,
                    "verdict": verdict,
                }
                steps.append(step)

                # Append accepted step to the working document
                append_step_to_working_doc(step, len(steps))
                write_log(log_entries, logpath)
                accepted = True
                break
            else:
                print(f"  [Verifier] REJECTED")
                rejected_proposals.append(
                    f"PROPOSAL:\n{proposal}\n\nREJECTION REASON:\n{verdict}"
                )

        if not accepted:
            print(f"  Round {round_num}: all attempts rejected, stopping.")
            break

    write_log(log_entries, logpath)
    print(f"\nFinished after {len(steps)} simplification(s).")
    print(f"  Working doc: {WORKING_DOC}")
    print(f"  Full log:    {logpath}")


if __name__ == "__main__":
    asyncio.run(main())
