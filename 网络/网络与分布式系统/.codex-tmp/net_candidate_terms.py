from collections import Counter
from pathlib import Path
import re


ROOT = Path(__file__).resolve().parent
text = (ROOT / "networds" / "all_text.txt").read_text(encoding="utf-8", errors="ignore")

token_re = re.compile(r"\b[A-ZÄÖÜ][A-Za-zÄÖÜäöüß0-9/-]{2,}\b")
tokens = token_re.findall(text)

stop = {
    "Die", "Der", "Das", "Eine", "Ein", "Sie", "Ihre", "Beispiel", "Abbildung",
    "Tabelle", "Kapitel", "Sommersemester", "Sommer", "LMU", "Fakultät",
    "Institut", "LFE", "Folgende", "Hinweis", "Aufgabe", "Aufgaben",
}

counter = Counter(t for t in tokens if t not in stop and not t.isdigit())
for term, count in counter.most_common(300):
    if count >= 2:
        print(f"{count:4} {term}")
