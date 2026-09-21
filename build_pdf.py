#!/usr/bin/env python3
"""Build THESIS.pdf from THESIS.md. Research deposit helper."""

from pathlib import Path

import markdown
from weasyprint import HTML

ROOT = Path(__file__).resolve().parent
src = (ROOT / "THESIS.md").read_text(encoding="utf-8")
body = markdown.markdown(
    src,
    extensions=["extra", "sane_lists", "tables", "nl2br"],
)
html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<title>Thesis #16. Phytochemical prioritisation under claim–evidence gates</title>
<style>
@page {{
  size: A4;
  margin: 22mm 18mm 24mm 18mm;
  @bottom-center {{
    content: "Ogbonna KE. Thesis #16. Computational research only. Not a medical device.  " counter(page);
    font-size: 8.5pt;
    color: #444;
    font-family: "DejaVu Serif", "Liberation Serif", "Times New Roman", Times, serif;
  }}
}}
html, body {{
  font-family: "DejaVu Serif", "Liberation Serif", "Times New Roman", Times, serif;
  font-size: 10.6pt;
  line-height: 1.42;
  color: #111;
}}
h1 {{
  font-size: 15.5pt;
  line-height: 1.25;
  margin: 0 0 0.8em 0;
  font-weight: 700;
}}
h2 {{
  font-size: 12.4pt;
  margin: 1.35em 0 0.45em 0;
  border-bottom: 0.4pt solid #333;
  padding-bottom: 0.15em;
}}
h3 {{
  font-size: 11.1pt;
  margin: 1.05em 0 0.35em 0;
}}
p {{ margin: 0.45em 0 0.55em 0; }}
strong {{ font-weight: 700; }}
code, pre {{
  font-family: "DejaVu Sans Mono", "Liberation Mono", monospace;
  font-size: 8.6pt;
}}
pre {{
  background: #f4f4f4;
  padding: 0.6em 0.7em;
  white-space: pre-wrap;
}}
table {{
  border-collapse: collapse;
  width: 100%;
  margin: 0.7em 0 1em 0;
  font-size: 8.4pt;
}}
th, td {{
  border: 0.4pt solid #555;
  padding: 0.22em 0.32em;
  vertical-align: top;
}}
th {{ background: #eee; text-align: left; }}
blockquote {{
  margin: 0.6em 0 0.6em 0.4em;
  padding-left: 0.8em;
  border-left: 2pt solid #555;
}}
ul, ol {{ margin: 0.3em 0 0.6em 1.3em; }}
li {{ margin: 0.15em 0; }}
hr {{ border: none; border-top: 0.4pt solid #999; margin: 1.2em 0; }}
img {{
  max-width: 100%;
  height: auto;
  display: block;
  margin: 0.6em auto 0.8em auto;
}}
.eq {{
  margin: 0.45em 0 0.45em 1.4em;
  font-style: italic;
}}
</style>
</head>
<body>
{body}
</body>
</html>
"""
out = ROOT / "THESIS.pdf"
HTML(string=html, base_url=str(ROOT)).write_pdf(out)
print(f"wrote {out} ({out.stat().st_size} bytes)")
