# Lightweight reproduce keyword matching for Klausur 2019 Frage 10
from pathlib import Path
import re
raw=Path(r'D:\lmustudy\网络与分布式系统\考卷\2019.md').read_text(encoding='utf-8', errors='replace')
heads=list(re.finditer(r'(?m)^(###\s+[^\n]+|##\s+(?!I\.|II\.|III\.|IV\.|V\.|VI\.|VII\.|VIII\.|IX\.|X\.|1 |2 |3 |4 |5 |6 |7 |8 |9 )[A-Za-zÄÖÜäöü].*)', raw))
for i,h in enumerate(heads):
    body=raw[h.start():(heads[i+1].start() if i+1<len(heads) else len(raw))]
    title=body.splitlines()[0].lstrip('#').strip()
    if 'Frage 10' in title:
        hay=title+'\n'+body
        keys=['Verzoeger','Verz枚ger','Verzöger','Paketvermittlung','Leitungsvermittlung','Signalver','Warteschlangen','Traffic Intensity','Intensit']
        print(title)
        for k in keys:
            if k.lower() in hay.lower(): print('match', k)
        print(body[:800].encode('unicode_escape').decode())
