#!/usr/bin/env python3
"""Generate llms.txt and llms-full.txt from Doxygen XML output and guide markdown.

Produces two files in docs/html/:
  - llms.txt       Concise index (library name, description, links)
  - llms-full.txt  Complete API reference + tutorial content

Requires: Doxygen XML output in docs/xml/ (GENERATE_XML = YES in Doxyfile).
Uses only Python 3 standard library.
"""

import os
import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
XML_DIR = REPO_ROOT / "docs" / "xml"
HTML_DIR = REPO_ROOT / "docs" / "html"
GUIDES_DIR = REPO_ROOT / "guides"
HEADER_XML = XML_DIR / "minidsp_8h.xml"

# Directories where \snippet source files live (matches EXAMPLE_PATH in Doxyfile)
EXAMPLE_PATHS = [
    REPO_ROOT / "examples",
    REPO_ROOT / "tools" / "audio_steg",
    REPO_ROOT / "tools" / "resample",
]


# ---------------------------------------------------------------------------
# XML text extraction
# ---------------------------------------------------------------------------

def xml_text(elem):
    """Recursively extract readable text from a Doxygen XML element."""
    if elem is None:
        return ""
    parts = []
    _xml_walk(elem, parts)
    return "".join(parts).strip()


def _xml_walk(elem, parts):
    """Walk XML tree, converting Doxygen elements to markdown text."""
    tag = elem.tag

    if tag == "para":
        _walk_children(elem, parts)
        parts.append("\n\n")
    elif tag == "formula":
        text = (elem.text or "").strip()
        if text.startswith("\\["):
            # Display math: \[...\] -> $$...$$
            inner = text[2:]
            if inner.endswith("\\]"):
                inner = inner[:-2]
            parts.append(f"\n$${inner.strip()}$$\n")
        elif text.startswith("$") and text.endswith("$"):
            # Inline math: $...$ -> keep as-is
            parts.append(text)
        else:
            parts.append(text)
    elif tag == "computeroutput":
        parts.append(f"`{elem.text or ''}`")
    elif tag == "bold":
        inner_parts = []
        _walk_children(elem, inner_parts)
        parts.append(f"**{''.join(inner_parts)}**")
    elif tag == "emphasis":
        inner_parts = []
        _walk_children(elem, inner_parts)
        parts.append(f"*{''.join(inner_parts)}*")
    elif tag == "ref":
        # Cross-reference: just use the text
        _walk_children(elem, parts)
    elif tag == "ulink":
        url = elem.get("url", "")
        inner_parts = []
        _walk_children(elem, inner_parts)
        parts.append(f"[{''.join(inner_parts)}]({url})")
    elif tag == "programlisting":
        code = _extract_code(elem)
        parts.append(f"\n```c\n{code}\n```\n")
    elif tag == "parameterlist":
        kind = elem.get("kind", "param")
        if kind == "param":
            parts.append("\n**Parameters:**\n")
        for item in elem.findall("parameteritem"):
            names = [n.text or "" for n in item.findall(".//parametername")]
            desc_elem = item.find("parameterdescription")
            desc = xml_text(desc_elem) if desc_elem is not None else ""
            for name in names:
                parts.append(f"- `{name}` - {desc}\n")
        parts.append("\n")
    elif tag == "simplesect":
        kind = elem.get("kind", "")
        inner_parts = []
        for child in elem:
            if child.tag != "title":
                _xml_walk(child, inner_parts)
        inner = "".join(inner_parts).strip()
        if kind == "return":
            parts.append(f"\n**Returns:** {inner}\n\n")
        elif kind == "note":
            parts.append(f"\n**Note:** {inner}\n\n")
        elif kind == "see":
            parts.append(f"\n**See also:** {inner}\n\n")
        elif kind == "par":
            title_elem = elem.find("title")
            title = title_elem.text if title_elem is not None else ""
            parts.append(f"\n**{title}** {inner}\n\n")
        else:
            parts.append(inner)
    elif tag == "itemizedlist":
        for li in elem.findall("listitem"):
            inner = xml_text(li)
            parts.append(f"- {inner}\n")
        parts.append("\n")
    elif tag == "orderedlist":
        for i, li in enumerate(elem.findall("listitem"), 1):
            inner = xml_text(li)
            parts.append(f"{i}. {inner}\n")
        parts.append("\n")
    elif tag in ("sp",):
        parts.append(" ")
    elif tag == "ndash":
        parts.append("--")
    elif tag == "mdash":
        parts.append("---")
    elif tag == "linebreak":
        parts.append("\n")
    elif tag == "anchor":
        pass  # skip anchors
    elif tag in ("sect1", "sect2", "sect3"):
        title_elem = elem.find("title")
        level = int(tag[-1])
        prefix = "#" * (level + 1)
        if title_elem is not None:
            parts.append(f"\n{prefix} {title_elem.text or ''}\n\n")
        _walk_children(elem, parts)
    elif tag == "title":
        pass  # handled by sect* parents
    elif tag in ("highlight", "codeline"):
        # These appear inside programlisting, handled by _extract_code
        _walk_children(elem, parts)
    else:
        # Default: walk children, preserving text
        _walk_children(elem, parts)

    # Tail text (text after closing tag, before next sibling)
    if elem.tail:
        parts.append(elem.tail)


def _walk_children(elem, parts):
    """Walk children of an element, including leading text."""
    if elem.text:
        parts.append(elem.text)
    for child in elem:
        _xml_walk(child, parts)


def _extract_code(programlisting):
    """Extract plain code text from a <programlisting> element."""
    lines = []
    for codeline in programlisting.findall("codeline"):
        line_parts = []
        _code_walk(codeline, line_parts)
        lines.append("".join(line_parts))
    return "\n".join(lines)


def _code_walk(elem, parts):
    """Walk a codeline/highlight tree extracting plain text."""
    if elem.tag == "sp":
        parts.append(" ")
        # sp elements have no children or meaningful text
    elif elem.tag in ("codeline", "highlight"):
        if elem.text:
            parts.append(elem.text)
        for child in elem:
            _code_walk(child, parts)
    elif elem.tag == "ref":
        # Cross-references in code: just use the text
        if elem.text:
            parts.append(elem.text)
    else:
        if elem.text:
            parts.append(elem.text)
    # Tail text belongs to the parent's sequence, not this element
    if elem.tail and elem.tag != "codeline":
        parts.append(elem.tail)


# ---------------------------------------------------------------------------
# API extraction from Doxygen XML
# ---------------------------------------------------------------------------

def extract_api(xml_path):
    """Extract all public API items from the minidsp.h XML file."""
    tree = ET.parse(xml_path)
    root = tree.getroot()
    compounddef = root.find("compounddef")

    sections = []
    for sectiondef in compounddef.findall("sectiondef"):
        kind = sectiondef.get("kind", "")
        members = []
        for memberdef in sectiondef.findall("memberdef"):
            mk = memberdef.get("kind", "")
            static = memberdef.get("static", "no")
            # Only include public, non-static items declared in minidsp.h
            loc = memberdef.find("location")
            if loc is not None:
                declfile = loc.get("declfile", "")
                if "minidsp.h" not in declfile:
                    continue
            if static == "yes":
                continue
            members.append(_extract_member(memberdef, mk))
        if members:
            sections.append((kind, members))
    return sections


def _extract_member(memberdef, kind):
    """Extract a single member definition into a dict."""
    name = (memberdef.findtext("name") or "").strip()
    result = {"kind": kind, "name": name}

    if kind == "function":
        ret_type = (memberdef.findtext("type") or "").strip()
        # Clean up ref tags in type
        type_elem = memberdef.find("type")
        if type_elem is not None:
            ret_type = "".join(type_elem.itertext()).strip()
        argsstring = (memberdef.findtext("argsstring") or "").strip()
        result["signature"] = f"{ret_type} {name}{argsstring}"
        result["definition"] = (memberdef.findtext("definition") or "").strip()
    elif kind == "enum":
        values = []
        for ev in memberdef.findall("enumvalue"):
            ev_name = (ev.findtext("name") or "").strip()
            ev_init = (ev.findtext("initializer") or "").strip()
            ev_brief = xml_text(ev.find("briefdescription"))
            values.append({"name": ev_name, "init": ev_init, "brief": ev_brief})
        result["values"] = values
    elif kind == "define":
        init = (memberdef.findtext("initializer") or "").strip()
        result["value"] = init
    elif kind == "typedef":
        definition = (memberdef.findtext("definition") or "").strip()
        result["definition"] = definition

    brief_elem = memberdef.find("briefdescription")
    result["brief"] = xml_text(brief_elem) if brief_elem is not None else ""

    detail_elem = memberdef.find("detaileddescription")
    result["detail"] = xml_text(detail_elem) if detail_elem is not None else ""

    return result


def format_api_section(sections):
    """Format extracted API sections into markdown text."""
    lines = []
    lines.append("# miniDSP API Reference\n\n")
    lines.append("C library for audio DSP. Header: `#include \"minidsp.h\"`\n\n")

    for kind, members in sections:
        for m in members:
            mk = m["kind"]
            name = m["name"]

            if mk == "function":
                sig = m.get("signature", name)
                lines.append(f"## `{sig}`\n\n")
                if m["brief"]:
                    lines.append(f"{m['brief']}\n\n")
                if m["detail"]:
                    detail = _clean_detail(m["detail"], m["brief"])
                    if detail:
                        lines.append(f"{detail}\n\n")
            elif mk == "enum":
                lines.append(f"## `enum {name}`\n\n")
                if m["brief"]:
                    lines.append(f"{m['brief']}\n\n")
                if m["detail"]:
                    lines.append(f"{m['detail']}\n\n")
                values = m.get("values", [])
                if values:
                    for v in values:
                        init = f" {v['init']}" if v["init"] else ""
                        brief = f" - {v['brief']}" if v["brief"] else ""
                        lines.append(f"- `{v['name']}{init}`{brief}\n")
                    lines.append("\n")
            elif mk == "define":
                val = m.get("value", "")
                if val:
                    lines.append(f"## `#define {name} {val}`\n\n")
                else:
                    lines.append(f"## `#define {name}`\n\n")
                if m["brief"]:
                    lines.append(f"{m['brief']}\n\n")
                if m["detail"]:
                    lines.append(f"{m['detail']}\n\n")
            elif mk == "typedef":
                defn = m.get("definition", "")
                if defn:
                    lines.append(f"## `{defn}`\n\n")
                else:
                    lines.append(f"## `{name}`\n\n")
                if m["brief"]:
                    lines.append(f"{m['brief']}\n\n")
                if m["detail"]:
                    lines.append(f"{m['detail']}\n\n")

            lines.append("---\n\n")

    return "".join(lines)


def _clean_detail(detail, brief):
    """Remove duplicate brief from detail text and clean up whitespace."""
    detail = detail.strip()
    brief = brief.strip()
    # Doxygen XML often duplicates the brief at the start of detailed
    if brief and detail.startswith(brief):
        detail = detail[len(brief):].strip()

    # Doxygen merges .h and .c doc-comments into one detaileddescription.
    # The .c version typically repeats the brief + adds implementation notes.
    # Remove the second copy: split on the brief text if it appears again.
    if brief and brief in detail:
        idx = detail.index(brief)
        # Keep only content before the duplicate brief
        detail = detail[:idx].strip()

    # Collapse excessive blank lines
    detail = re.sub(r"\n{3,}", "\n\n", detail)
    return detail


# ---------------------------------------------------------------------------
# Guide markdown processing
# ---------------------------------------------------------------------------

def get_guide_order():
    """Parse guides/tutorials.md for \\subpage order."""
    tutorials_path = GUIDES_DIR / "tutorials.md"
    if not tutorials_path.exists():
        return []
    text = tutorials_path.read_text()
    # Match \subpage page-id patterns
    pages = re.findall(r"\\subpage\s+([\w-]+)", text)
    return pages


def resolve_snippet(filename, snippet_id):
    """Find and extract code between //! [id] markers in a source file."""
    for base_dir in EXAMPLE_PATHS:
        filepath = base_dir / filename
        if filepath.exists():
            text = filepath.read_text()
            marker = f"//! [{snippet_id}]"
            parts = text.split(marker)
            if len(parts) >= 3:
                # Content is between first and second marker
                code = parts[1].strip()
                return code
    print(f"  Warning: unresolved snippet '{snippet_id}' in '{filename}'",
          file=sys.stderr)
    return None


def strip_doxygen_syntax(text):
    """Convert Doxygen-specific markdown syntax to standard markdown."""
    # Display math: \f[...\f] -> $$...$$
    text = re.sub(
        r"\\f\[(.*?)\\f\]",
        lambda m: f"\n$${m.group(1).strip()}$$\n",
        text,
        flags=re.DOTALL,
    )
    # Inline math: \f$...\f$ -> $...$
    text = re.sub(r"\\f\$(.*?)\\f\$", r"$\1$", text)

    # Remove \htmlonly ... \endhtmlonly blocks
    text = re.sub(
        r"\\htmlonly.*?\\endhtmlonly",
        "",
        text,
        flags=re.DOTALL,
    )

    # Resolve \snippet directives
    def resolve_snippet_match(m):
        filename = m.group(1)
        snippet_id = m.group(2)
        code = resolve_snippet(filename, snippet_id)
        if code is not None:
            return f"\n```c\n{code}\n```\n"
        return f"*(snippet: {filename} [{snippet_id}])*"

    text = re.sub(
        r"\\snippet\s+(\S+)\s+(.+)",
        resolve_snippet_match,
        text,
    )

    # Strip {#anchor-id} from headings
    text = re.sub(r"\s*\{#[\w-]+\}", "", text)

    # Convert \subpage page-id — description -> just the description text
    text = re.sub(
        r"\\subpage\s+[\w-]+\s*(?:—\s*)?",
        "",
        text,
    )

    # Convert \ref page-id to just the text
    text = re.sub(r"\\ref\s+([\w-]+)", r"\1", text)

    # Strip @p parameter references (just keep the name)
    text = re.sub(r"@p\s+(\w+)", r"`\1`", text)

    # Strip @c code references
    text = re.sub(r"@c\s+(\w+)", r"`\1`", text)

    # Remove Doxygen @note, @see, @par directives (keep content)
    text = re.sub(r"@note\s+", "**Note:** ", text)
    text = re.sub(r"@see\s+", "**See also:** ", text)

    return text


def process_guide(page_id):
    """Read and process a guide markdown file."""
    # Find the markdown file for this page ID
    # Page IDs use hyphens, filenames match
    md_path = GUIDES_DIR / f"{page_id}.md"
    if not md_path.exists():
        print(f"  Warning: guide file not found for '{page_id}'", file=sys.stderr)
        return None

    text = md_path.read_text()
    text = strip_doxygen_syntax(text)

    # Collapse excessive blank lines
    text = re.sub(r"\n{3,}", "\n\n", text)

    return text.strip()


# ---------------------------------------------------------------------------
# Output generation
# ---------------------------------------------------------------------------

def generate_llms_txt():
    """Generate the concise llms.txt index file."""
    return """# miniDSP

> A small C library for audio DSP: signal measurement, FFT spectrum analysis, biquad filtering, GCC-PHAT delay estimation, and audio I/O.

miniDSP provides building blocks for audio processing pipelines: signal generators, window functions, FFT-based spectrum analysis, mel filterbanks, MFCCs, FIR filtering, DTMF detection, audio steganography, sample rate conversion, and more.

## Documentation

- [Full API Reference and Tutorials](llms-full.txt): Complete function signatures, parameter docs, formulas, code examples, and tutorial content in a single file
- Source: https://github.com/wooters/miniDSP
- HTML Docs: https://wooters.github.io/miniDSP/
"""


def generate_llms_full_txt(api_text, guide_texts):
    """Generate the complete llms-full.txt file."""
    parts = []

    parts.append("# miniDSP - Complete Reference\n\n")
    parts.append("> A small C library for audio DSP\n\n")
    parts.append("This document contains the complete API reference and tutorial ")
    parts.append("content for the miniDSP library. It is generated automatically ")
    parts.append("from the source code documentation.\n\n")
    parts.append("Source: https://github.com/wooters/miniDSP\n\n")
    parts.append("---\n\n")

    # API Reference
    parts.append(api_text)
    parts.append("\n---\n\n")

    # Tutorials
    parts.append("# Tutorials\n\n")
    for title, text in guide_texts:
        parts.append(text)
        parts.append("\n\n---\n\n")

    return "".join(parts)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if not XML_DIR.exists():
        print(f"Error: Doxygen XML directory not found: {XML_DIR}", file=sys.stderr)
        print("Run 'doxygen Doxyfile' first (with GENERATE_XML = YES).",
              file=sys.stderr)
        sys.exit(1)

    if not HEADER_XML.exists():
        print(f"Error: {HEADER_XML} not found.", file=sys.stderr)
        sys.exit(1)

    HTML_DIR.mkdir(parents=True, exist_ok=True)

    # Extract API reference from XML
    print("Extracting API reference from Doxygen XML...")
    sections = extract_api(HEADER_XML)
    api_text = format_api_section(sections)

    # Count functions for reporting
    func_count = sum(
        1 for _, members in sections
        for m in members if m["kind"] == "function"
    )
    print(f"  Found {func_count} public functions")

    # Process guide markdown files
    print("Processing guide markdown files...")
    page_order = get_guide_order()
    guide_texts = []
    for page_id in page_order:
        text = process_guide(page_id)
        if text:
            guide_texts.append((page_id, text))
            print(f"  Processed: {page_id}")

    # Generate output files
    llms_txt = generate_llms_txt()
    llms_full_txt = generate_llms_full_txt(api_text, guide_texts)

    out_llms = HTML_DIR / "llms.txt"
    out_full = HTML_DIR / "llms-full.txt"

    out_llms.write_text(llms_txt)
    out_full.write_text(llms_full_txt)

    llms_size = len(llms_txt.encode("utf-8"))
    full_size = len(llms_full_txt.encode("utf-8"))
    print(f"\nGenerated:")
    print(f"  {out_llms}  ({llms_size:,} bytes)")
    print(f"  {out_full}  ({full_size:,} bytes)")


if __name__ == "__main__":
    main()
