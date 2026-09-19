from collections import Counter
from pathlib import Path
import csv


ROOT = Path(__file__).resolve().parents[1]
path = ROOT / "networds.csv"

with path.open("r", encoding="utf-8", newline="") as handle:
    rows = list(csv.reader(handle))

bad = [(index + 1, len(row), row) for index, row in enumerate(rows) if len(row) != 3]
terms = [row[0] for row in rows if row]
duplicates = [term for term, count in Counter(terms).items() if count > 1]

required = {
    "ARP", "MAC-Adresse", "VLAN", "CSMA/CD", "CRC", "Hammingabstand",
    "Lichtwellenleiter", "Manchester-Codierung", "Slotted Aloha", "Switch",
}
missing = sorted(required - set(terms))

print(f"rows={len(rows)}")
print(f"bad_rows={len(bad)}")
print(f"duplicates={len(duplicates)}")
print(f"missing_required={len(missing)}")
if bad:
    print("bad=" + repr(bad[:5]))
if duplicates:
    print("duplicates_list=" + repr(duplicates[:20]))
if missing:
    print("missing_list=" + repr(missing))
