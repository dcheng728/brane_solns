"""
Adding equation
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import List, Optional, Tuple


class EquationNumberer:
    _num_suffix = re.compile(
        r"^(?P<body>.*?)(?P<ws>\s*)\\qquad\{(?P<lpar>\()?(?P<num>\d+)(?P<rpar>\))?\}\s*$"
    )

    def __init__(self, start: int = 1) -> None:
        self._start = start

    def process_text(self, text: str) -> str:
        newline = self._detect_newline(text)
        lines = text.splitlines(keepends=True)

        out: List[str] = []
        in_block = False
        in_html_comment = False
        block_lines: List[str] = []
        eq_number = self._start - 1

        for line in lines:
            if in_html_comment:
                out.append(line)
                if "-->" in line:
                    in_html_comment = False
                continue

            if "<!--" in line:
                out.append(line)
                if "-->" not in line:
                    in_html_comment = True
                continue

            if line.strip() == "$$":
                if not in_block:
                    in_block = True
                    block_lines = []
                    out.append(line)
                else:
                    eq_number += 1
                    out.extend(self._number_block(block_lines, eq_number, newline))
                    out.append(line)
                    in_block = False
                    block_lines = []
                continue

            if in_block:
                block_lines.append(line)
            else:
                out.append(line)

        if in_block:
            out.extend(block_lines)

        return "".join(out)

    def process_file(self, path: Path, in_place: bool = False) -> str:
        original = path.read_text(encoding="utf-8")
        updated = self.process_text(original)

        if in_place and updated != original:
            path.write_text(updated, encoding="utf-8", newline="")

        return updated

    def _number_block(self, block_lines: List[str], eq_number: int, newline: str) -> List[str]:
        if not block_lines:
            return [f"\\qquad{{({eq_number})}}{newline}"]

        last_nonempty = self._last_nonempty_index(block_lines)
        if last_nonempty is None:
            return [f"\\qquad{{({eq_number})}}{newline}"] + block_lines

        raw, eol = self._strip_eol(block_lines[last_nonempty])

        m = self._num_suffix.match(raw)
        if m:
            existing = int(m.group("num"))
            has_parens = (m.group("lpar") is not None) and (m.group("rpar") is not None)
            if existing == eq_number and has_parens:
                return block_lines
            body = m.group("body")
            ws = m.group("ws")
            block_lines[last_nonempty] = f"{body}{ws}\\qquad{{({eq_number})}}{eol}"
            return block_lines

        insert_eol = eol or newline
        insertion = f"\\qquad{{({eq_number})}}{insert_eol}"
        return block_lines[: last_nonempty + 1] + [insertion] + block_lines[last_nonempty + 1 :]

    @staticmethod
    def _detect_newline(text: str) -> str:
        return "\r\n" if "\r\n" in text else "\n"

    @staticmethod
    def _strip_eol(line: str) -> Tuple[str, str]:
        if line.endswith("\r\n"):
            return line[:-2], "\r\n"
        if line.endswith("\n"):
            return line[:-1], "\n"
        return line, ""

    @classmethod
    def _last_nonempty_index(cls, lines: List[str]) -> Optional[int]:
        for i in range(len(lines) - 1, -1, -1):
            if cls._strip_eol(lines[i])[0].strip() != "":
                return i
        return None


def main() -> int:
    ap = argparse.ArgumentParser(description="Add/normalize equation numbers in $$...$$ Markdown blocks.")
    ap.add_argument("path", type=Path)
    ap.add_argument("--start", type=int, default=1)
    ap.add_argument("--in-place", action="store_true")
    args = ap.parse_args()

    numberer = EquationNumberer(start=args.start)
    result = numberer.process_file(args.path, in_place=args.in_place)

    if not args.in_place:
        print(result, end="")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())