from pathlib import Path
import re


ROOT = Path(__file__).resolve().parents[1]
text = (ROOT / ".codex-tmp" / "networds" / "all_text.txt").read_text(
    encoding="utf-8", errors="ignore"
)
match = re.search(
    r"===== Lecture/Kapitel 5 Hardwarenah\.pdf =====\n(?P<text>.*?)(?:\n\n=====|\Z)",
    text,
    re.S,
)
chapter = match.group("text") if match else ""
lines = [line.strip() for line in chapter.splitlines()]

skip = {
    "Rechnernetze und verteilte Systeme",
    "(Prof. Dr. D. Kranzlmüller)",
    "Sommersemester 2026",
    "Ludwig-Maximilians-Universität München",
    "Prof. Dr. D. Kranzlmüller",
    "Anwendung",
    "Transport",
    "Vermittlung",
    "Netzanschluss",
}

for i, line in enumerate(lines):
    if not line or line in skip:
        continue
    if re.fullmatch(r"\d+", line):
        continue
    if line.startswith(("•", "→", "=", "⎫", "⎬", "⎭")):
        continue
    if len(line) > 80:
        continue
    if any(key in line for key in [
        "Sicherungsschicht", "Schicht 2", "MAC", "ARP", "LAN", "Aloha",
        "CSMA", "MACA", "Switch", "VLAN", "Fehler", "Parität", "CRC",
        "Hamming", "Ethernet", "Bitübertragung", "Abtast", "Codierung",
        "Modulation", "Medien", "Dämpfung", "Lichtwellenleiter", "Zusammenfassung",
        "Aufgaben", "Adressierung", "Kollision", "Vielfachzugriff", "Rahmen",
    ]):
        print(f"{i+1:04}: {line}")
