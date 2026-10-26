#!/usr/bin/env python3

"""Compare selected OpenFOAM dictionaries between two cases.

The script compares these files:
  - system/controlDict
  - system/fvSolution
  - system/fvSchemes
  - system/inputParameters
  - 0/mapOFBnd
  - 0/U
  - 0/p
  - 0/alpha
  - 0/Ur
  - 0/k
  - 0/omega
  - 0/nut
  
The script compares all files in the folder, but not in its subfolders:
  - constant

It parses OpenFOAM dictionary syntax, ignores comments and formatting-only
differences, and reports only keyword/value differences.

Examples:
      - cmpOFcases.py caseA caseB
      - cmpOFcases.py caseB
      - cmpOFcases.py caseA caseB -o myCompare.log

Usage:
    $FOAM_INST_DIR/BakerFoam/utils/forOFCase/cmpOFcases.py caseA caseB
    $FOAM_INST_DIR/BakerFoam/utils/forOFCase/cmpOFcases.py caseB
    $FOAM_INST_DIR/BakerFoam/utils/forOFCase/cmpOFcases.py caseA caseB -o myCompare.log

"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence


TARGET_FILES = (
    "system/fvSolution",
    "system/fvSchemes",
    "system/controlDict",
    "system/inputParameters",
    "0/mapOFBnd",
    "0/U",
    "0/p",
    "0/alpha",
    "0/Ur",
    "0/k",
    "0/omega",
    "0/nut",
)

BLOCK_MESH_FILE = "system/blockMeshDict"
CONSTANT_DIR = "constant"

SYMBOLS = {"{", "}", "(", ")", "[", "]", ";", ","}


@dataclass
class Node:
    kind: str
    value: str = ""
    children: Optional[Dict[str, "Node"]] = None


def strip_comments_keep_layout(text: str) -> str:
    result: List[str] = []
    i = 0
    length = len(text)
    in_string = False

    while i < length:
        char = text[i]
        next_char = text[i + 1] if i + 1 < length else ""

        if in_string:
            result.append(char)
            if char == '"' and (i == 0 or text[i - 1] != "\\"):
                in_string = False
            i += 1
            continue

        if char == '"':
            in_string = True
            result.append(char)
            i += 1
            continue

        if char == "/" and next_char == "/":
            while i < length and text[i] != "\n":
                i += 1
            continue

        if char == "/" and next_char == "*":
            i += 2
            while i + 1 < length and not (text[i] == "*" and text[i + 1] == "/"):
                if text[i] == "\n":
                    result.append("\n")
                i += 1
            i += 2
            continue

        result.append(char)
        i += 1

    return "".join(result)


def tokenize(text: str) -> List[str]:
    tokens: List[str] = []
    i = 0
    length = len(text)

    while i < length:
        char = text[i]

        if char in " \t\r":
            i += 1
            continue

        if char == "\n":
            tokens.append("\n")
            i += 1
            continue

        if char == '"':
            start = i
            i += 1
            while i < length:
                if text[i] == '"' and text[i - 1] != "\\":
                    i += 1
                    break
                i += 1
            tokens.append(text[start:i])
            continue

        if char == "#" and i + 1 < length and text[i + 1] == "{":
            tokens.append("#{")
            i += 2
            continue

        if char == "#" and i + 2 < length and text[i + 1] == "}" and text[i + 2] == ";":
            tokens.append("#};")
            i += 3
            continue

        if char == "#":
            start = i
            i += 1
            while i < length and re.match(r"[A-Za-z_]", text[i]):
                i += 1
            tokens.append(text[start:i])
            continue

        if char in SYMBOLS:
            tokens.append(char)
            i += 1
            continue

        start = i
        while i < length:
            if text[i] in " \t\r\n{}[];,":
                break

            if text[i] == "(":
                if i == start:
                    break

                depth = 0
                while i < length:
                    if text[i] == "(":
                        depth += 1
                    elif text[i] == ")":
                        depth -= 1
                        if depth == 0:
                            i += 1
                            break
                    i += 1
                continue

            i += 1

        tokens.append(text[start:i])

    return tokens


def normalize_value_tokens(tokens: Sequence[str]) -> str:
    compact: List[str] = []
    for token in tokens:
        if token == "\n":
            continue
        compact.append(token)

    joined = " ".join(compact)
    joined = re.sub(r"\s+([;,(){}\[\],])", r"\1", joined)
    joined = re.sub(r"([({\[] )", lambda m: m.group(1).strip(), joined)
    joined = re.sub(r"\s+", " ", joined).strip()
    return joined


class FoamDictParser:
    def __init__(self, tokens: Sequence[str]):
        self.tokens = list(tokens)
        self.index = 0

    def current(self) -> Optional[str]:
        if self.index >= len(self.tokens):
            return None
        return self.tokens[self.index]

    def advance(self) -> Optional[str]:
        token = self.current()
        if token is not None:
            self.index += 1
        return token

    def skip_newlines(self) -> None:
        while self.current() == "\n":
            self.index += 1

    def parse(self) -> Node:
        children = self.parse_entries(stop_token=None)
        return Node(kind="dict", children=children)

    def parse_entries(self, stop_token: Optional[str]) -> Dict[str, Node]:
        entries: Dict[str, Node] = {}

        while self.index < len(self.tokens):
            self.skip_newlines()
            token = self.current()

            if token is None:
                break

            if stop_token is not None and token == stop_token:
                self.advance()
                break

            if token == ";":
                self.advance()
                continue

            if token.startswith("#"):
                self.skip_directive()
                continue

            key = self.advance()
            if key is None:
                break

            self.skip_newlines()
            token = self.current()

            if token == "{":
                self.advance()
                entries[key] = Node(kind="dict", children=self.parse_entries(stop_token="}"))
                continue

            value_tokens = self.collect_value_tokens()
            entries[key] = Node(kind="value", value=normalize_value_tokens(value_tokens))

        return entries

    def skip_directive(self) -> None:
        while self.index < len(self.tokens):
            token = self.advance()
            if token == "\n":
                break

    def collect_value_tokens(self) -> List[str]:
        value_tokens: List[str] = []
        round_depth = 0
        square_depth = 0
        code_depth = 0

        while self.index < len(self.tokens):
            token = self.current()
            if token is None:
                break

            if token == "#;":
                self.advance()
                break

            if token == "#};" and code_depth > 0:
                value_tokens.append(self.advance())
                code_depth -= 1
                if code_depth == 0:
                    break
                continue

            if token == ";" and round_depth == 0 and square_depth == 0 and code_depth == 0:
                self.advance()
                break

            if token == "(":
                round_depth += 1
            elif token == ")" and round_depth > 0:
                round_depth -= 1
            elif token == "[":
                square_depth += 1
            elif token == "]" and square_depth > 0:
                square_depth -= 1
            elif token == "#{":
                code_depth += 1

            value_tokens.append(self.advance())

        return value_tokens


def parse_foam_dictionary(path: Path) -> Node:
    text = path.read_text(encoding="utf-8", errors="ignore")
    stripped = strip_comments_keep_layout(text)
    tokens = tokenize(stripped)
    return FoamDictParser(tokens).parse()


def compare_raw_text(path_a: Path, path_b: Path) -> List[str]:
    """Fallback comparison for files that are not parsed as dictionaries."""
    text_a = path_a.read_text(encoding="utf-8", errors="ignore")
    text_b = path_b.read_text(encoding="utf-8", errors="ignore")

    norm_a = "\n".join(
        line.strip() for line in strip_comments_keep_layout(text_a).splitlines() if line.strip()
    )
    norm_b = "\n".join(
        line.strip() for line in strip_comments_keep_layout(text_b).splitlines() if line.strip()
    )

    if norm_a == norm_b:
        return []

    return [
        "raw content differs after comment/whitespace normalization\n"
        "  case A: <different content>\n"
        "  case B: <different content>"
    ]


def list_relative_files(root: Path) -> List[Path]:
    if not root.is_dir():
        return []
    return sorted(p.relative_to(root) for p in root.iterdir() if p.is_file())


def prune_control_dict(node: Node) -> Node:
    """Remove heavy function blocks from controlDict comparison."""
    if node.kind != "dict" or not node.children:
        return node

    pruned_children: Dict[str, Node] = {}
    for key, child in node.children.items():
        if key == "functions":
            continue
        pruned_children[key] = child

    return Node(kind="dict", children=pruned_children)


def compare_nodes(
    relative_file: str,
    node_a: Optional[Node],
    node_b: Optional[Node],
    path_parts: Optional[List[str]] = None,
) -> List[str]:
    differences: List[str] = []
    path_parts = path_parts or []
    path_label = "/".join(path_parts) if path_parts else "<root>"

    if node_a is None and node_b is None:
        return differences

    if node_a is None:
        differences.append(
            f"{relative_file}: {path_label}\n"
            f"  case A: <missing>\n"
            f"  case B: {render_node(node_b)}"
        )
        return differences

    if node_b is None:
        differences.append(
            f"{relative_file}: {path_label}\n"
            f"  case A: {render_node(node_a)}\n"
            f"  case B: <missing>"
        )
        return differences

    if node_a.kind != node_b.kind:
        differences.append(
            f"{relative_file}: {path_label}\n"
            f"  case A: {render_node(node_a)}\n"
            f"  case B: {render_node(node_b)}"
        )
        return differences

    if node_a.kind == "value":
        if node_a.value != node_b.value:
            differences.append(
                f"{relative_file}: {path_label}\n"
                f"  case A: {node_a.value}\n"
                f"  case B: {node_b.value}"
            )
        return differences

    keys = sorted(set((node_a.children or {}).keys()) | set((node_b.children or {}).keys()))
    for key in keys:
        child_a = (node_a.children or {}).get(key)
        child_b = (node_b.children or {}).get(key)
        differences.extend(compare_nodes(relative_file, child_a, child_b, path_parts + [key]))

    return differences


def render_node(node: Optional[Node]) -> str:
    if node is None:
        return "<missing>"
    if node.kind == "value":
        return node.value
    keys = sorted((node.children or {}).keys())
    return "{ " + ", ".join(keys) + " }"


def sanitize_name(path: Path) -> str:
    return path.name.replace(" ", "_")


def split_common_prefix(name_a: str, name_b: str) -> tuple[str, str, str]:
    """Split two names into common prefix and unique suffixes.

    Example:
      annularFeifeiHighEccCoarse, annularFeifeiHighEccRef
      -> (annularFeifeiHighEcc, Coarse, Ref)
    """
    max_len = min(len(name_a), len(name_b))
    i = 0
    while i < max_len and name_a[i] == name_b[i]:
        i += 1

    prefix = name_a[:i]
    suffix_a = name_a[i:]
    suffix_b = name_b[i:]

    # Trim separator tails/heads when names use separators.
    prefix = prefix.rstrip("_-")
    suffix_a = suffix_a.lstrip("_-")
    suffix_b = suffix_b.lstrip("_-")

    # If one side equals the common prefix exactly, keep the prefix and mark
    # that side as "Base" to avoid repeating long names.
    if not suffix_a and suffix_b:
        return prefix, "Base", suffix_b

    if not suffix_b and suffix_a:
        return prefix, suffix_a, "Base"

    # If both suffixes are empty, names are identical.
    if not suffix_a and not suffix_b:
        return "", name_a, name_b

    return prefix, suffix_a, suffix_b


def default_output_name(case_a: Path, case_b: Path) -> str:
    """Build default output filename.

    Format:
      CMP_<prefix>_<suffixA>_VS_<suffixB>.log
      or CMP_<caseA>_VS_<caseB>.log when no useful prefix exists.
    """
    name_a = sanitize_name(case_a)
    name_b = sanitize_name(case_b)
    prefix, suffix_a, suffix_b = split_common_prefix(name_a, name_b)

    if prefix:
        return f"CMP_{prefix}_{suffix_a}_VS_{suffix_b}.log"

    return f"CMP_{name_a}_VS_{name_b}.log"


def block_mesh_differs(case_a: Path, case_b: Path) -> bool:
    file_a = case_a / BLOCK_MESH_FILE
    file_b = case_b / BLOCK_MESH_FILE

    if not file_a.exists() and not file_b.exists():
        return False

    if not file_a.exists() or not file_b.exists():
        return True

    text_a = file_a.read_text(encoding="utf-8", errors="ignore")
    text_b = file_b.read_text(encoding="utf-8", errors="ignore")
    return text_a != text_b


def build_log(
    case_a: Path,
    case_b: Path,
    results: Dict[str, List[str]],
    constant_results: Dict[str, List[str]],
    block_mesh_warning: bool,
) -> str:
    lines = [
        f"Case comparison",
        f"case A: {case_a}",
        f"case B: {case_b}",
        "",
    ]

    if block_mesh_warning:
        lines.extend(
            [
                f"WARNING: {BLOCK_MESH_FILE} files are different.",
                "They were compared as whole files, not parsed by dictionary keywords.",
                "",
            ]
        )

    total = 0
    for relative_file in TARGET_FILES:
        differences = results[relative_file]
        if not differences:
            continue
        lines.append(relative_file)
        lines.append("-" * len(relative_file))
        lines.extend(differences)
        lines.append("")
        total += len(differences)

    constant_keys = sorted(constant_results.keys())
    for relative_file in constant_keys:
        differences = constant_results[relative_file]
        if not differences:
            continue
        lines.append(relative_file)
        lines.append("-" * len(relative_file))
        lines.extend(differences)
        lines.append("")
        total += len(differences)

    if total == 0:
        lines.append("No keyword value differences found.")
    else:
        lines.append(f"Total differences: {total}")

    return "\n".join(lines) + "\n"


def compare_case_files(case_a: Path, case_b: Path) -> Dict[str, List[str]]:
    results: Dict[str, List[str]] = {}

    for relative_file in TARGET_FILES:
        file_a = case_a / relative_file
        file_b = case_b / relative_file

        if not file_a.exists() and not file_b.exists():
            results[relative_file] = []
            continue

        if not file_a.exists():
            results[relative_file] = [
                f"{relative_file}: <file>\n  case A: <missing file>\n  case B: present"
            ]
            continue

        if not file_b.exists():
            results[relative_file] = [
                f"{relative_file}: <file>\n  case A: present\n  case B: <missing file>"
            ]
            continue

        node_a = parse_foam_dictionary(file_a)
        node_b = parse_foam_dictionary(file_b)

        if relative_file == "system/controlDict":
            node_a = prune_control_dict(node_a)
            node_b = prune_control_dict(node_b)

        results[relative_file] = compare_nodes(relative_file, node_a, node_b)

    return results


def compare_constant_folder(case_a: Path, case_b: Path) -> Dict[str, List[str]]:
    results: Dict[str, List[str]] = {}

    root_a = case_a / CONSTANT_DIR
    root_b = case_b / CONSTANT_DIR

    if not root_a.exists() and not root_b.exists():
        return results

    if not root_a.exists():
        results[CONSTANT_DIR] = [
            f"{CONSTANT_DIR}: <folder>\n  case A: <missing folder>\n  case B: present"
        ]
        return results

    if not root_b.exists():
        results[CONSTANT_DIR] = [
            f"{CONSTANT_DIR}: <folder>\n  case A: present\n  case B: <missing folder>"
        ]
        return results

    files_a = set(list_relative_files(root_a))
    files_b = set(list_relative_files(root_b))

    for rel in sorted(files_a | files_b):
        relative_file = f"{CONSTANT_DIR}/{rel.as_posix()}"
        file_a = root_a / rel
        file_b = root_b / rel

        if not file_a.exists():
            results[relative_file] = [
                f"{relative_file}: <file>\n  case A: <missing file>\n  case B: present"
            ]
            continue

        if not file_b.exists():
            results[relative_file] = [
                f"{relative_file}: <file>\n  case A: present\n  case B: <missing file>"
            ]
            continue

        # Fast path for byte-identical files.
        if file_a.read_bytes() == file_b.read_bytes():
            results[relative_file] = []
            continue

        try:
            node_a = parse_foam_dictionary(file_a)
            node_b = parse_foam_dictionary(file_b)
            results[relative_file] = compare_nodes(relative_file, node_a, node_b)
        except Exception:
            results[relative_file] = compare_raw_text(file_a, file_b)

    return results


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compare OpenFOAM case dictionaries and save a differences-only log.",
    )
    parser.add_argument(
        "case_a",
        nargs="?",
        help=(
            "Relative or absolute path to the first case. "
            "If omitted and case_b is provided, current directory is used as first case."
        ),
    )
    parser.add_argument(
        "case_b",
        nargs="?",
        help="Relative or absolute path to the second case",
    )
    parser.add_argument(
        "-o",
        "--output",
        help=(
            "Optional log file path. "
            "Default: CMP_<commonPrefix>_<suffixA>_VS_<suffixB>.log in current directory"
        ),
    )
    return parser.parse_args()


def resolve_case_path(raw: str) -> Path:
    """Resolve case path with support for sibling-case invocation.

        This allows running from inside one case folder using:
            cmpOFcases.py caseA caseB
    where caseA and caseB are sibling directories.
    """
    candidate = Path(raw).expanduser()

    if candidate.is_absolute():
        return candidate.resolve()

    direct = (Path.cwd() / candidate).resolve()
    if direct.is_dir():
        return direct

    sibling = (Path.cwd().parent / candidate).resolve()
    if sibling.is_dir():
        return sibling

    # Fall back to direct-resolved path so error messages remain explicit.
    return direct


def main() -> int:
    args = parse_args()
    if args.case_a and args.case_b:
        case_a_raw = args.case_a
        case_b_raw = args.case_b
    elif args.case_a and not args.case_b:
        case_a_raw = str(Path.cwd())
        case_b_raw = args.case_a
    else:
        raise SystemExit("Usage: cmpOFcases.py <caseB> OR cmpOFcases.py <caseA> <caseB>")

    case_a = resolve_case_path(case_a_raw)
    case_b = resolve_case_path(case_b_raw)

    if not case_a.is_dir():
        raise SystemExit(f"First case path is not a directory: {case_a}")

    if not case_b.is_dir():
        raise SystemExit(f"Second case path is not a directory: {case_b}")

    if args.output:
        output_path = Path(args.output).expanduser().resolve()
    else:
        output_name = default_output_name(case_a, case_b)
        output_path = Path.cwd() / output_name

    results = compare_case_files(case_a, case_b)
    constant_results = compare_constant_folder(case_a, case_b)
    output_text = build_log(
        case_a,
        case_b,
        results,
        constant_results,
        block_mesh_differs(case_a, case_b),
    )
    output_path.write_text(output_text, encoding="utf-8")

    print(f"Wrote comparison log: {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
