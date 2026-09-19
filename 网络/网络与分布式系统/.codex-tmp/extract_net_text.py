from pathlib import Path
import json
import re

import pdfplumber


ROOT = Path(__file__).resolve().parents[1]
TARGETS = [ROOT / "作业", ROOT / "Lecture"]
OUT_DIR = ROOT / ".codex-tmp" / "networds"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def clean_text(text: str) -> str:
    text = text.replace("\x00", " ")
    text = re.sub(r"-\s*\n\s*", "", text)
    text = re.sub(r"[ \t]+", " ", text)
    return text


records = []
all_text_parts = []

for base in TARGETS:
    for path in sorted(base.rglob("*")):
        if not path.is_file():
            continue
        if path.suffix.lower() == ".pdf":
            pages = []
            with pdfplumber.open(path) as pdf:
                for page in pdf.pages:
                    pages.append(page.extract_text() or "")
            text = clean_text("\n".join(pages))
        elif path.suffix.lower() in {".md", ".txt"}:
            text = clean_text(path.read_text(encoding="utf-8", errors="ignore"))
        else:
            continue
        rel = path.relative_to(ROOT).as_posix()
        records.append({"file": rel, "chars": len(text)})
        all_text_parts.append(f"\n\n===== {rel} =====\n{text}")

(OUT_DIR / "documents.json").write_text(
    json.dumps(records, ensure_ascii=False, indent=2), encoding="utf-8"
)
(OUT_DIR / "all_text.txt").write_text("\n".join(all_text_parts), encoding="utf-8")

print(json.dumps(records, ensure_ascii=False, indent=2))
