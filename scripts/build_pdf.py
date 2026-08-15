#!/usr/bin/env python3
"""
Generate a high-quality PDF from docs/TUTORIAL.md using Headless Chromium,
Marked.js, Highlight.js, and MathJax 3.
"""

import os
import subprocess
import json

def build_pdf():
    tutorial_path = os.path.abspath("docs/TUTORIAL.md")
    with open(tutorial_path, "r", encoding="utf-8") as f:
        md_content = f.read()

    # Escape markdown content for embedding into JS template literal
    escaped_md = json.dumps(md_content)

    html_template = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Clifford Algebra Tutorial</title>
  
  <!-- MathJax for rendering LaTeX math formulas -->
  <script>
  window.MathJax = {{
    tex: {{
      inlineMath: [['$', '$'], ['\\\\(', '\\\\)']],
      displayMath: [['$$', '$$'], ['\\\\[', '\\\\]']]
    }},
    svg: {{ fontCache: 'global' }}
  }};
  </script>
  <script src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js"></script>

  <!-- Marked.js for markdown parsing -->
  <script src="https://cdn.jsdelivr.net/npm/marked/marked.min.js"></script>

  <!-- Highlight.js for code syntax highlighting -->
  <link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/github.min.css">
  <script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/highlight.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/languages/haskell.min.js"></script>
  <script src="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/languages/bash.min.js"></script>

  <!-- Google Fonts for typography -->
  <link rel="preconnect" href="https://fonts.googleapis.com">
  <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
  <link href="https://fonts.googleapis.com/css2?family=Fira+Code:wght@400;500&family=Inter:wght@300;400;500;600;700&family=Lora:ital,wght@0,400;0,600;1,400&display=swap" rel="stylesheet">

  <style>
    @page {{
      size: A4 portrait;
      margin: 20mm 18mm 22mm 18mm;
      @bottom-right {{
        content: counter(page);
      }}
    }}

    body {{
      font-family: 'Inter', -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
      color: #1a1a2e;
      background-color: #ffffff;
      line-height: 1.65;
      font-size: 10.5pt;
      margin: 0;
      padding: 0;
    }}

    h1, h2, h3, h4 {{
      color: #0f172a;
      font-weight: 700;
      line-height: 1.3;
      page-break-after: avoid;
    }}

    h1 {{
      font-size: 20pt;
      border-bottom: 2px solid #3b82f6;
      padding-bottom: 8px;
      margin-top: 0;
      margin-bottom: 16px;
    }}

    h2 {{
      font-size: 14pt;
      border-bottom: 1px solid #e2e8f0;
      padding-bottom: 6px;
      margin-top: 24px;
      margin-bottom: 12px;
    }}

    h3 {{
      font-size: 12pt;
      margin-top: 18px;
      margin-bottom: 8px;
    }}

    p {{
      margin-top: 0;
      margin-bottom: 10px;
      text-align: justify;
    }}

    code {{
      font-family: 'Fira Code', monospace;
      font-size: 9pt;
      background-color: #f1f5f9;
      padding: 2px 5px;
      border-radius: 4px;
      color: #0f172a;
    }}

    pre {{
      background-color: #f8fafc !important;
      border: 1px solid #e2e8f0;
      border-left: 4px solid #3b82f6;
      border-radius: 6px;
      padding: 12px 14px;
      margin: 12px 0;
      overflow-x: auto;
      page-break-inside: avoid;
    }}

    pre code {{
      background-color: transparent !important;
      padding: 0;
      font-size: 9pt;
      line-height: 1.5;
    }}

    table {{
      width: 100%;
      border-collapse: collapse;
      margin: 14px 0;
      font-size: 9.5pt;
      page-break-inside: avoid;
    }}

    th, td {{
      padding: 8px 12px;
      border: 1px solid #cbd5e1;
      text-align: left;
    }}

    th {{
      background-color: #f1f5f9;
      font-weight: 600;
      color: #0f172a;
    }}

    tr:nth-child(even) {{
      background-color: #f8fafc;
    }}

    blockquote {{
      margin: 12px 0;
      padding: 8px 16px;
      border-left: 4px solid #3b82f6;
      background-color: #eff6ff;
      color: #1e3a8a;
      border-radius: 0 6px 6px 0;
    }}

    ul, ol {{
      margin-top: 0;
      margin-bottom: 10px;
      padding-left: 20px;
    }}

    li {{
      margin-bottom: 4px;
    }}

    .MathJax_SVG {{
      display: inline-block !important;
      vertical-align: middle !important;
    }}

    mjx-container[display="true"] {{
      margin: 12px 0 !important;
    }}
  </style>
</head>
<body>
  <div id="content"></div>

  <script>
    const markdownText = {escaped_md};
    
    // Configure marked to use highlight.js
    marked.setOptions({{
      highlight: function(code, lang) {{
        const language = hljs.getLanguage(lang) ? lang : 'plaintext';
        return hljs.highlight(code, {{ language }}).value;
      }},
      gfm: true,
      breaks: false
    }});

    document.getElementById('content').innerHTML = marked.parse(markdownText);

    // Trigger MathJax rendering
    if (window.MathJax) {{
      MathJax.typesetPromise();
    }}
  </script>
</body>
</html>
"""

    html_path = os.path.abspath("docs/TUTORIAL.html")
    pdf_path = os.path.abspath("docs/TUTORIAL.pdf")

    with open(html_path, "w", encoding="utf-8") as f:
        f.write(html_template)

    print(f"Generated HTML at: {html_path}")
    print("Converting HTML to PDF via Headless Chromium...")

    cmd = [
        "chromium",
        "--headless=new",
        "--disable-gpu",
        "--no-sandbox",
        "--run-all-compositor-stages-before-draw",
        "--virtual-time-budget=5000",
        f"--print-to-pdf={pdf_path}",
        f"file://{html_path}"
    ]

    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode == 0 and os.path.exists(pdf_path):
        size_kb = os.path.getsize(pdf_path) / 1024
        print(f"Successfully generated PDF: {pdf_path} ({size_kb:.1f} KB)")
        return pdf_path
    else:
        print(f"Error rendering PDF: {res.stderr}")
        return None

if __name__ == "__main__":
    build_pdf()
