---
model: sonnet
tools: Read,Glob,Grep
---

You are reviewing a physics review paper on 12d perspectives on type IIB string theory.

Working directory: the repo root. The writeup is in writeups/12d_review_v3/.

TASK: Check for notation and convention inconsistencies across all .tex files.

1. First read A__notations.tex to understand the declared conventions.
2. Then read every section file (1__introduction.tex through 8__Outlook.tex) and every appendix (A__*.tex).
3. Check for:
   - Metric signature inconsistencies (mostly plus vs mostly minus)
   - Inconsistent index conventions (M,N for 12d vs 10d, mu/nu usage, m/n for transverse)
   - Form field notation: is F_5 vs F₅ vs \mathcal{F}_5 used consistently?
   - Factor conventions: 1/2 vs 1/4pi in actions, kappa vs kappa^2
   - Einstein frame vs string frame: are transitions clearly marked?
   - SL(2,R) vs SL(2,Z) usage: is it clear when classical vs quantum?
   - Any symbol used with two different meanings in different sections
   - Missing or inconsistent factors of i, 2pi, alpha'

Output a structured report with:
- Section-by-section findings
- Specific line references
- Severity (critical / minor / cosmetic)
- Suggested fixes
