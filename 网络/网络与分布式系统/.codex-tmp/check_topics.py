from pathlib import Path
import re
p=Path(r'D:\lmustudy\网络与分布式系统\分知识点汇总')
missing=[]
for f in sorted(p.glob('[0-9][0-9]_*.md')):
    s=f.read_text(encoding='utf-8-sig')
    q=len(re.findall(r'(?m)^### 题目\s+\d+:', s))
    src=s.count('来源说明：')
    zh=s.count('#### 题目中文翻译 / 中文题意')
    de=s.count('#### 德文原题')
    sol=s.count('#### 解答')
    details=s.count('<details')
    imgs=s.count('pictures/')
    print(f.name.encode('unicode_escape').decode(), 'q', q, 'src', src, 'zh', zh, 'de', de, 'sol', sol, 'details', details, 'img_refs', imgs)
    if not (q==src==zh==de==sol) or details:
        missing.append((f.name,q,src,zh,de,sol,details))
img_missing=[]
for f in p.glob('[0-9][0-9]_*.md'):
    s=f.read_text(encoding='utf-8-sig')
    for m in re.findall(r'!\[[^\]]*\]\((pictures/[^)]+)\)', s):
        if not (p/m).exists(): img_missing.append((f.name,m))
print('mismatch', missing)
print('missing_images', img_missing[:10], 'count', len(img_missing))
