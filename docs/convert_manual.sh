#!/bin/bash
#
# Convert a chapter of the LaTeX manual into a Documenter Markdown page.
#
#   Usage: docs/convert_manual.sh <chapter.tex> <output.md>
#   e.g.   docs/convert_manual.sh ch3-tutorial.tex docs/src/tutorial.md
#
# Requires pandoc. The manual source is expected under docs/manual/pcpop_manual/
# (that folder is git-ignored — keep the Overleaf source there locally).
#
# What the pipeline does:
#   - expands the custom LaTeX macros from commands.tex (fixing one broken \norm def),
#   - display math ($$...$$) -> ```math fences, \mathds -> \mathbb (KaTeX),
#   - strips \label / \nonumber and the author's colour-coded editorial notes,
#   - code environments -> ```julia blocks,
#   - citations [@key] / [@a; @b] -> Documenter [key](@cite) / [a, b](@cite).
#
# Manual follow-up still needed: \ref/\eqref cross-references, citations that
# carry a locator (e.g. "[@key Algorithm 4.1]"), and figures.
set -e

HERE="$(cd "$(dirname "$0")" && pwd)"
MAN="$HERE/manual/pcpop_manual"
CH="$1"
OUT="$2"

[ -n "$CH" ] && [ -n "$OUT" ] || { echo "usage: $0 <chapter.tex> <output.md>"; exit 1; }
command -v pandoc >/dev/null || { echo "pandoc not found (brew install pandoc)"; exit 1; }

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

# commands.tex has a broken \norm definition (\left\| with no \right\|) that
# breaks pandoc's parser; repair it before expanding macros.
sed 's/\\newcommand{\\norm}\[1\]{\\left\\|#1 }/\\newcommand{\\norm}[1]{\\left\\|#1\\right\\|}/' \
    "$MAN/commands.tex" > "$TMP/commands.tex"
cat "$TMP/commands.tex" "$MAN/$CH" > "$TMP/in.tex"

pandoc "$TMP/in.tex" --from latex --to markdown --wrap=none 2>/dev/null \
| perl -0777 -pe '
    s/\A.*?(^# )/$1/ms;                          # drop tcolorbox preamble before first heading
    s/\$\$(.*?)\$\$/\n```math\n$1\n```\n/gs;      # display math -> ```math fences
    s/\\mathds/\\mathbb/g;                        # KaTeX has no \mathds
    s/\\label\{[^}]*\}//g;                        # KaTeX chokes on \label
    s/\\nonumber//g;
' \
| perl -pe '
    s/\[\\\[.*?\\\]\]\{style="color[^}]*"\}//g;   # editorial notes with nested [\[..\]]
    s/\[[^\[\]]*\]\{style="color[^}]*"\}//g;       # simple editorial notes
    s{\[\@([\w:.-]+(?:\s*;\s*\@[\w:.-]+)*)\]}{ (my $k=$1) =~ s/\s*;\s*\@/, /g; "[$k](\@cite)"; }ge;  # citations
    s/\{reference-type[^}]*\}//g;                  # pandoc \ref/\eqref attributes
    s/^:::+.*$//;                                  # pandoc div fences
    s/^```\s*\{[^}]*style="code"[^}]*\}/```julia/; # tag code blocks julia
    s/\s*\{#[^}]*\}//g;                            # \label anchors
' > "$TMP/body.md"

{ printf '```@meta\nCurrentModule = PCPOP\n```\n\n'; cat "$TMP/body.md"; } > "$OUT"
echo "wrote $OUT ($(wc -l < "$OUT") lines)"
