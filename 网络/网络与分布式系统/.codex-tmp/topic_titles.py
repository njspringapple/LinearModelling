from pathlib import Path
import re
p=Path(r'D:\lmustudy\网络与分布式系统\分知识点汇总')
for f in sorted(p.glob('[0-9][0-9]_*.md')):
    print('\n##', f.name.encode('unicode_escape').decode())
    s=f.read_text(encoding='utf-8-sig')
    for line in re.findall(r'^### 题目 \d+: (.*)$', s, flags=re.M):
        print('-', line[:120].encode('unicode_escape').decode())
