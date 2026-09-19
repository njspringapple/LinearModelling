from pathlib import Path
import importlib.util
spec=importlib.util.spec_from_file_location('gen', r'D:\lmustudy\网络与分布式系统\.codex-tmp\gen_topics.py')
# Avoid import because script executes. Instead inspect generated delay file snippets around suspicious question titles.
p=Path(r'D:\lmustudy\网络与分布式系统\分知识点汇总\06_延迟_分组交换与电路交换.md')
s=p.read_text(encoding='utf-8-sig')
for target in ['Frage 10', 'Frage 19', 'Frage 22']:
    i=s.find(target)
    print('\nTARGET', target, 'pos', i)
    print(s[i-300:i+500].encode('unicode_escape').decode() if i!=-1 else 'not found')
