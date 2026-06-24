# Sphinx configuration for the CRISPRware documentation site.
project = "CRISPRware"
author = "Eric Malekos"
copyright = "2026, Eric Malekos"
release = "0.2"

extensions = [
    "myst_parser",            # author pages in Markdown
    "sphinx_copybutton",      # copy button on code blocks
    "sphinx_design",          # cards / grids / tabs
    "sphinxcontrib.mermaid",  # pipeline diagrams
]

myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "attrs_inline",
    "tasklist",
    "fieldlist",
]
myst_heading_anchors = 3

source_suffix = {".md": "markdown"}
exclude_patterns = ["build", "_build", "_code", "Thumbs.db", ".DS_Store"]

html_theme = "sphinx_rtd_theme"
html_title = "CRISPRware"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "navigation_depth": 3,
    "collapse_navigation": False,
    "sticky_navigation": True,
    "titles_only": False,
}

mermaid_init_js = (
    "mermaid.initialize({startOnLoad:true, theme:'neutral', "
    "flowchart:{curve:'basis', nodeSpacing:45, rankSpacing:45, useMaxWidth:true}});"
)
