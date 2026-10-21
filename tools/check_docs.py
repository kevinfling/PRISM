#!/usr/bin/env python3
"""Check local links, anchors, and known GitHub math pitfalls with the standard library.

This is a source check, not a complete TeX parser or a substitute for rendering.
"""
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
FILES = [ROOT / name for name in ("README.md", "CONTRIBUTING.md", "CHANGELOG.md")]
FILES += sorted((ROOT / "docs").rglob("*.md"))
FILES += sorted((ROOT / "results").rglob("*.md"))


def anchors(path):
    counts = {}
    found = set()
    fenced = False
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith("```"):
            fenced = not fenced
        if fenced:
            continue
        match = re.match(r"^#{1,6}\s+(.+?)\s*#*\s*$", line)
        if not match:
            continue
        heading = re.sub(r"[^\w\s-]", "", match.group(1).lower())
        slug = heading.replace(" ", "-")
        index = counts.get(slug, 0)
        counts[slug] = index + 1
        found.add(slug if index == 0 else f"{slug}-{index}")
    return found


def math_fragments(text):
    """Extract protected math and flag markup that can lose TeX escapes on GitHub."""
    expressions, errors = [], []
    fence, language, start, block = None, None, 0, []
    dollar_block = False
    for number, line in enumerate(text.splitlines(), 1):
        marker = re.match(r"^\s*(`{3,}|~{3,})(.*)$", line)
        if fence:
            if marker and marker[1][0] == fence[0] and len(marker[1]) >= len(fence) and not marker[2].strip():
                if language == "math":
                    expressions.append((start, "\n".join(block)))
                fence = None
            elif language == "math":
                block.append(line)
            continue
        if line.strip() == "$$":
            if dollar_block:
                expressions.append((start, "\n".join(block)))
            else:
                start, block = number, []
                errors.append(f"line {number}: use a fenced math block to preserve TeX escapes")
            dollar_block = not dollar_block
            continue
        if dollar_block:
            block.append(line)
            continue
        if marker:
            fence, language, start, block = marker[1], marker[2].strip(), number, []
            continue
        # Preserve math in dollar/backtick delimiters; ignore ordinary code spans.
        def token(match):
            if match[1] is not None:
                expressions.append((number, match[1]))
            return " " * len(match[0])

        remaining = re.sub(r"\$`([^`\n]+)`\$|(`+)(.*?)\2", token, line)
        if re.search(r"(?<!\\)\$", remaining):
            errors.append(f"line {number}: use dollar/backtick delimiters for inline math")
    if fence or dollar_block:
        errors.append(f"line {start}: unclosed code or math block")
    for number, expression in expressions:
        if re.search(r"\\operatorname\b", expression):
            errors.append(f"line {number}: GitHub rejects the operatorname macro; use mathop with mathrm")
        depth = 0
        for token in re.findall(r"\\[A-Za-z]+|\\.|[{}]", expression):
            if token == "{":
                depth += 1
            elif token == "}":
                depth -= 1
            if depth < 0:
                break
        if depth:
            errors.append(f"line {number}: unbalanced TeX braces")
    return expressions, errors


def main():
    errors = []
    checked, equations = 0, 0
    for source in FILES:
        if not source.is_file():
            errors.append(f"{source.relative_to(ROOT)}: missing required documentation file")
            continue
        text = source.read_text(encoding="utf-8")
        expressions, math_errors = math_fragments(text)
        equations += len(expressions)
        errors.extend(f"{source.relative_to(ROOT)}: {error}" for error in math_errors)
        if sum(line.startswith("```") for line in text.splitlines()) % 2:
            errors.append(f"{source.relative_to(ROOT)}: unbalanced code fences")
        for link in re.findall(r"(?<!!)\[[^\]]*\]\(([^\s)]+)\)", text):
            url = urlsplit(link)
            if url.scheme or url.netloc:
                continue
            checked += 1
            target = (source.parent / unquote(url.path)).resolve() if url.path else source
            try:
                target.relative_to(ROOT)
            except ValueError:
                errors.append(f"{source.relative_to(ROOT)}: link escapes repository: {link}")
                continue
            if not target.exists():
                errors.append(f"{source.relative_to(ROOT)}: missing link: {link}")
            elif url.fragment and target.suffix == ".md" and unquote(url.fragment) not in anchors(target):
                errors.append(f"{source.relative_to(ROOT)}: missing section: {link}")
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"Checked {checked} local links and {equations} math expressions across {len(FILES)} Markdown files")
    return 0


if __name__ == "__main__":
    sys.exit(main())
