from collections import Counter, defaultdict
from pathlib import Path
import csv
import re


ROOT = Path(__file__).resolve().parents[1]
TEXT = (ROOT / ".codex-tmp" / "networds" / "all_text.txt").read_text(
    encoding="utf-8", errors="ignore"
)
EXISTING = {
    row[0]
    for row in csv.reader((ROOT / "networds.csv").open("r", encoding="utf-8", newline=""))
    if row
}

sections = re.split(r"\n===== (?P<file>.+?) =====\n", TEXT)
by_file = {}
for index in range(1, len(sections), 2):
    by_file[sections[index]] = sections[index + 1]

patterns = {
    "upper": re.compile(r"\b[A-Z][A-Z0-9]{1,}(?:/[A-Z0-9]+)?\b"),
    "hyphen": re.compile(r"\b[A-Za-zÄÖÜäöüß]+(?:-[A-Za-zÄÖÜäöüß0-9]+)+\b"),
    "compound": re.compile(
        r"\b[A-ZÄÖÜ][A-Za-zÄÖÜäöüß]*(?:adresse|adressen|algorithmus|bereich|bit|bits|block|code|codierung|dienst|dienste|domäne|einheit|fehler|fenster|fluss|funktion|header|kanal|knoten|kontrolle|korrektur|layer|leitung|medium|modell|netz|netzwerk|paket|phase|protokoll|rahmen|rate|router|schicht|schnittstelle|segment|signal|socket|system|tabelle|technik|verfahren|vermittlung|verzögerung|zugriff)\b"
    ),
    "english_phrase": re.compile(
        r"\b(?:Address|Application|Binary|Carrier|Collision|Data|Distance|Dynamic|Fast|Flow|Internet|Link|Media|Network|Open|Packet|Physical|Request|Selective|Service|Sliding|Slow|Source|Stop|Transport|Well)\s+[A-Z][A-Za-z]+(?:\s+[A-Z][A-Za-z]+)?\b"
    ),
}

stop = {
    "Prof", "Kranzlmüller", "Rechnernetze", "Systeme", "Sommersemester",
    "Kapitel", "Abbildung", "Seite", "München", "Daniel", "Fabian",
    "Dreer", "Ludwig-Maximilians-Universität", "RNVS", "Alice", "Bob",
    "Carol", "Moodle", "Übungsblatt", "Aufgabe", "Aufgaben",
}

for name, text in by_file.items():
    print(f"\n## {name}")
    seen = defaultdict(Counter)
    for kind, pattern in patterns.items():
        for match in pattern.findall(text):
            term = " ".join(match.split())
            if term in stop or term in EXISTING or len(term) < 3:
                continue
            seen[kind][term] += 1
    for kind, counter in seen.items():
        useful = [
            (term, count)
            for term, count in counter.most_common(80)
            if count >= (1 if kind in {"upper", "hyphen", "english_phrase"} else 2)
        ]
        if useful:
            print(f"\n{kind}:")
            print(", ".join(f"{term}({count})" for term, count in useful[:60]))
