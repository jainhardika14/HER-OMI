import markdown, re
md = open("STUDY_GUIDE.md").read()
body = markdown.markdown(md, extensions=["tables", "fenced_code", "md_in_html"])
css = open("guide.css").read()
html = f'<!doctype html><html><head><meta charset="utf-8"><title>HER-OMI study guide</title><style>{css}</style></head><body>{body}</body></html>'
open("STUDY_GUIDE.html", "w").write(html)
