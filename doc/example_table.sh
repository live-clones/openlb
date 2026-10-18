#!/usr/bin/env bash
#
# Pre-render helper for the OpenLB user guide.
#
# The tabular overview of the examples is maintained in a Google spreadsheet
# and exported from there as a web page, _sheet/Table_google.html. That
# export is a standalone document: it carries the spreadsheet frame -- the
# A, B, C column letters and the 1, 2, 3 row numbers -- and it links Google's
# own stylesheet, _sheet/resources/sheet.css, which is several megabytes of editor
# CSS of which the table needs three rules.
#
# This script turns the export into examples/_table_google.qmd, a raw HTML
# block that examples/index_web.qmd includes, so the overview shows up on the
# Examples page itself instead of behind a link to a separate document:
#
#   * the <thead> with the column letters is dropped, but its widths are kept
#     by rewriting them into a <colgroup> -- the cells themselves carry no
#     width, so without it the fixed layout of the sheet would be lost;
#   * the row-number <th> that starts every row is dropped;
#   * the cell classes the export defines inline (.s0, .s1, ...) are kept as
#     they are and the handful of base rules the table needs from Google's
#     stylesheet are inlined, which makes the result self-contained and
#     _sheet/resources/ unnecessary;
#   * the table is wrapped in a horizontally scrollable container, because it
#     is far wider than the body text column of the website.
#
# Re-exporting the spreadsheet therefore means replacing _sheet/Table_google.html
# and rendering again -- nothing else. _sheet/ is named with a leading underscore
# so that Quarto neither renders nor publishes the raw export, and the generated
# file is listed in .gitignore.
# The script is registered as a pre-render script in _quarto.yml, so a plain
# `quarto render` regenerates it automatically. It is idempotent and can also
# be run by hand.

set -euo pipefail

quarto_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export_html="$quarto_dir/_sheet/Table_google.html"
partial="$quarto_dir/examples/_table_google.qmd"

if [ ! -f "$export_html" ]; then
  echo "example_table: no spreadsheet export at $export_html" >&2
  exit 1
fi

perl -0777 -e '
  my $html = do { local $/; <STDIN> };

  # The export keeps the cell styles (.ritz .waffle .sN) in an inline <style>
  # block and the sheet itself in a single wrapper div.
  my ($style) = $html =~ m{<style[^>]*>(.*?)</style>}s
    or die "example_table: no <style> block in the export\n";
  my ($grid) = $html =~ m{(<div class="ritz grid-container".*</div>)}s
    or die "example_table: no grid container in the export\n";

  # The column widths live on the header cells, so they have to be salvaged
  # into a <colgroup> before the header row goes.
  my ($thead) = $grid =~ m{<thead>(.*?)</thead>}s
    or die "example_table: no table header in the export\n";
  my @widths;
  while ($thead =~ m{<th\b[^>]*\bstyle="width:\s*(\d+)px;?"[^>]*class="column-headers-background"}g) {
    push @widths, $1;
  }
  @widths or die "example_table: no column widths in the table header\n";
  my $colgroup = "<colgroup>"
    . join("", map { "<col style=\"width:${_}px\">" } @widths)
    . "</colgroup>";

  $grid =~ s{<thead>.*?</thead>}{}s;                        # column letters
  $grid =~ s{<th\b[^>]*class="row-headers-background".*?</th>}{}gs;  # row numbers
  $grid =~ s{(<table\b[^>]*>)}{$1$colgroup}s;

  print <<"HTML";
```{=html}
<style>
/* Base rules the exported sheet expects from Google\x27s stylesheet, which is
   not shipped with the user guide. */
.ritz .waffle{border-collapse:collapse;table-layout:fixed;font-family:docs-Calibri,Calibri,Arial,sans-serif;font-size:11pt;color:#000000;}
.ritz .waffle td{overflow:hidden;border:1px solid #c4c7c5;border-width:0 1px 1px 0;vertical-align:bottom;line-height:inherit;background-color:#ffffff;padding:0 3px;}
.ritz .waffle .softmerge{overflow:visible;}
/* The sheet is much wider than the text column of the page. */
.ritz.grid-container{overflow-x:auto;max-width:100%;}
$style
</style>
$grid
```
HTML
' < "$export_html" > "$partial"

echo "example_table: wrote $(basename "$partial")"
