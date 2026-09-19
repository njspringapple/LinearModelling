from collections import Counter
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
all_text = (ROOT / ".codex-tmp" / "networds" / "all_text.txt").read_text(
    encoding="utf-8", errors="ignore"
)

match = re.search(
    r"===== Lecture/Kapitel 5 Hardwarenah\.pdf =====\n(?P<text>.*?)(?:\n\n=====|\Z)",
    all_text,
    re.S,
)
text = match.group("text") if match else ""

patterns = [
    r"\b[A-ZÄÖÜ][A-Za-zÄÖÜäöüß0-9/-]{2,}\b",
    r"\b[A-Z]{2,}[A-Z0-9/-]*\b",
    r"\b[a-zäöüß]+(?:-[A-Za-zÄÖÜäöüß0-9]+)+\b",
]

terms = []
for pattern in patterns:
    terms.extend(re.findall(pattern, text))

stop = {
    "Prof", "Kranzlmüller", "Rechnernetze", "Systeme", "Sommersemester",
    "Kapitel", "Hardwarenah", "Abbildung", "Seite", "München", "Daniel",
    "Diefenthaler", "Fabian", "Dreer", "Ludwig-Maximilians-Universität",
}

counter = Counter(t for t in terms if t not in stop and not t.isdigit())
for term, count in counter.most_common(240):
    if count >= 2:
        print(f"{count:4} {term}")
