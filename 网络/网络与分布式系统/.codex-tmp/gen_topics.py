from pathlib import Path
import re, hashlib, traceback
base = Path(r'D:\lmustudy\网络与分布式系统')
outdir = base / '分知识点汇总'
outdir.mkdir(exist_ok=True)

def clean_lines(text):
    return '\n'.join(line for line in text.splitlines() if not line.startswith('Seite ') and not line.startswith('--- PAGE'))

def norm(text):
    return re.sub(r'\s+', ' ', text).strip().lower()

def sanitize_markdown_fragment(text):
    """Keep rich markdown renderable, but prevent imported headings from polluting Obsidian outline."""
    lines = []
    in_fence = False
    for line in text.strip().splitlines():
        if line.strip().startswith('```'):
            in_fence = not in_fence
            lines.append(line)
            continue
        if not in_fence and re.match(r'^\s{0,3}#{1,6}\s+', line):
            label = re.sub(r'^\s{0,3}#{1,6}\s+', '', line).strip()
            lines.append(f'**{label}**')
        elif not in_fence and line.strip() in {'<details open>', '</details>', '<details>'}:
            continue
        elif not in_fence and line.strip().startswith('<summary>') and line.strip().endswith('</summary>'):
            continue
        else:
            lines.append(line)
    return '\n'.join(lines).strip()

def keyword_hit(hay, key):
    if re.fullmatch(r'[A-Z0-9]{2,6}', key):
        return re.search(r'(?i)(?<![A-Za-z0-9])' + re.escape(key) + r'(?![A-Za-z0-9])', hay) is not None
    if key in {'IPv4', 'IPv6', 'ICMPv6', 'IP-Header', 'SSHv2'}:
        return re.search(r'(?i)(?<![A-Za-z0-9])' + re.escape(key) + r'(?![A-Za-z0-9])', hay) is not None
    return key.lower() in hay.lower()

assignment_zh = {
    'Uebungsblatt 00, Aufgabe 1': '设计互联网连接结构：比较全互联网络和每个节点最多连接 5 个下级节点的层次结构，分别计算 8 个、300 个以及 N 个参与者所需连接数。',
    'Uebungsblatt 00, Aufgabe 2': '复习位值计数系统：写出十六进制数字，转换若干十进制数到十六进制、八进制和二进制，并分析 `2^32 - 1` 的二进制位数。',
    'Uebungsblatt 00, Aufgabe 3': '在二进制、八进制、十进制和十六进制中进行乘法和幂运算，理解基数幂在对应进制中的表示。',
    'Uebungsblatt 01, Aufgabe 1': '用两位古代学者通过信鸽远程下棋的例子，说明分布式系统中消息丢失、重复、延迟、乱序或损坏会带来的问题，并提出避免方法。',
    'Uebungsblatt 01, Aufgabe 2': '复习图和树：画 4 个节点的全互联图，列举 4 个节点树的形状，并计算二叉树在给定高度下的最少/最多节点数。',
    'Uebungsblatt 01, Aufgabe 3': '研究二叉树中的寻址：根据节点名或边名写路径，判断路径是否能说明节点类型，并把二进制路径转换成十进制、十六进制和字节元组。',
    'Uebungsblatt 01, Aufgabe 4': '把二进制数看成完整二叉树中的路径，计算给定前缀下有多少叶子节点。',
    'Uebungsblatt 02, Aufgabe 1': '用电话订披萨设计一个协议：画订餐顺序图，加入底层配送服务，区分控制数据和有效载荷，并讨论改用 Messenger 时分层如何变化。',
    'Uebungsblatt 02, Aufgabe 2': '判断给定系统是计算机网络还是分布式系统，并说明理由，例如 MWN、Messenger、WWW、SuperMUC-NG 等。',
    'Uebungsblatt 02, Aufgabe 3': '阅读 RFC 768，说明 UDP 描述了什么、核心特性是什么，并分析如果用 UDP 订披萨会遇到哪些可靠性问题。',
    'Uebungsblatt 02, Aufgabe 4': '练习 Linux 命令行、ping 和 traceroute：查看目录、理解 man page、测量往返时延并解释 traceroute 输出。',
    'Uebungsblatt 03, Aufgabe 1': '补全 OSI 七层模型，说明各层主要任务，分析分层架构优缺点，并比较 OSI 应用层与 Internet 模型应用层。',
    'Uebungsblatt 03, Aufgabe 2': '根据 Python `sendto(bytes, addr_tuple)` 示例，识别 Nutzdaten、控制信息、PDU/SDU/ICI，并解码字节消息。',
    'Uebungsblatt 03, Aufgabe 3': '解释 PDU 如何由 SDU 加 PCI 形成，说明对等实体、层间封装关系，并判断图中的 Dienstschnitt、Protokollschnitt 或 Systemschnitt。',
    'Uebungsblatt 03, Aufgabe 4': '阅读 OSI FAQ，分析 OSI 模型中有争议的层、OSI 与 TCP/IP 的主要批评点，以及 OSI 模型的教学价值。',
    'Uebungsblatt 04, Aufgabe 1': '比较无连接通信和面向连接通信，给出例子，并说明无连接通信在哪些情况下有优势。',
    'Uebungsblatt 04, Aufgabe 2': '分析 2-way handshake 的失败场景，画 3-way handshake 与连接释放状态/时序图，并解释为什么任何有限握手都不能绝对保证成功。',
    'Uebungsblatt 04, Aufgabe 3': '在 Stop-and-Wait 协议中使用序列号：画无错传输、ACK 丢失和数据损坏的时序图，并说明双方如何检测错误。',
    'Uebungsblatt 05, Aufgabe 1': '根据传输层伪代码识别 UDP 风格的 PDU、PCI、ICI、下层接口信息，并判断协议是否面向连接。',
    'Uebungsblatt 05, Aufgabe 2': '给定 TCP 已接收字节数、段长度和端口，计算后续段的序列号、ACK 号和端口，并画 ACK 丢失时序图。',
    'Uebungsblatt 05, Aufgabe 3': '在带发送窗口的协议中分析无错传输、消息丢失、NACK、乱序缓存和累计确认优化。',
    'Uebungsblatt 05, Aufgabe 4': '用 Wireshark 分析 TCP 抓包：找出三次握手、连接释放、RTD、绝对/相对序列号，并判断 SSHv2 是否使用 TCP。',
    'Uebungsblatt 06, Aufgabe 1': '画完整 TCP 请求-响应交换，包括三次握手、数据、ACK 和连接释放，并计算有连接/无连接情况下的时间差。',
    'Uebungsblatt 06, Aufgabe 2': '在 TCP 中使用选择确认 SACK：给定 8 个段和丢失模式，写出累计 ACK 与 SACK block。',
    'Uebungsblatt 06, Aufgabe 3': '根据 TCP Reno 的拥塞窗口图，识别 Slow Start、Congestion Avoidance、丢包轮次、threshold、Fast Recovery 和 Tahoe/Reno 差异。',
    'Uebungsblatt 06, Aufgabe 4': '推导 TCP 拥塞避免阶段的丢包率公式，并由丢包率近似 TCP 平均吞吐率。',
    'Uebungsblatt 07, Aufgabe 1': '说明处理、排队、发送和传播四类延迟出现在哪里，并计算家庭网络中一个 1500B 包从 A 到 B 的总延迟。',
    'Uebungsblatt 07, Aufgabe 2': '比较固定速率长期应用在电路交换和分组交换中的适用性，并分析即使平均容量足够时为何仍可能丢包。',
    'Uebungsblatt 07, Aufgabe 3': '在卫星链路中计算 Sliding Window 不同窗口大小下的有效传输率和信道利用率，并求满利用所需最小窗口。',
    'Uebungsblatt 07, Aufgabe 4': '在 TCP Slow Start 中给定 RTT、数据量、段大小和速率，画传输时序并比较有无 Slow Start 的总传输时间。',
    'Uebungsblatt 08, Aufgabe 1': '比较 CIDR 与分类地址，给 `131.42.0.0/16` 按不同主机需求划分子网，并写出掩码、可用范围和广播地址。',
    'Uebungsblatt 08, Aufgabe 2': '根据 ISP 层次结构，把 `160.229.0.0/16` 分配给不同路由器和下级子网。',
    'Uebungsblatt 08, Aufgabe 3': '对给定四路由器图从 A 运行 Dijkstra/SPF，画中间步骤、最终路由表，并分析链路故障后的最短路径树。',
    'Uebungsblatt 08, Aufgabe 4': '把 RFC1918 私有 IPv4 地址范围写成前缀形式，证明 `172.16.0.0/12` 可包含 16 个 `/16` 网络，并讨论私有地址优缺点。',
    'Uebungsblatt 08, Aufgabe 5': '在 NAT 场景中说明私网主机访问公网服务器时 IP/端口如何被改写，NAT 表如何转发返回包，并讨论 NAT 与安全/IPv6 的关系。',
    'Uebungsblatt 09, Aufgabe 1': '在九个路由器和子网 G 的拓扑中，按距离向量协议逐轮填写每个路由器通告到 G 的 hop 数，直到稳定。',
    'Uebungsblatt 09, Aufgabe 2': '说明 AS、IGP、EGP、BGP、路径向量、Routing Policy、Transit 和 Peering 的含义与差异。',
    'Uebungsblatt 09, Aufgabe 3': '在 IPv6 ISP 拓扑中分配链路前缀和路由器接口地址，并为 Router B 写路由表和默认路由。',
    'Uebungsblatt 09, Aufgabe 4': '判断 IPv6 地址是否合法，并在完整写法和最短写法之间转换。',
    'Uebungsblatt 09, Aufgabe 5': '在给定 MTU 路径中计算 IPv4 分片长度、标志和 offset，并比较 IPv6 下 Packet Too Big 与源主机分片。',
    'Uebungsblatt 09, Aufgabe 6': '分析距离向量路由中的 Count-to-Infinity：链路断开后距离如何逐步增加，以及 Split Horizon 如何缓解。',
    'Uebungsblatt 10, Aufgabe 1': '给三个 IPv4 子网和两个路由器分配 IP/MAC，分析 E 到 B 的转发过程中每一跳的源/目的 IP 与 MAC，并讨论 ARP 表为空时的流程。',
    'Uebungsblatt 10, Aufgabe 2': '设计实验区分未知设备是 Hub 还是 Switch，并测量交换机转发表老化时间。',
    'Uebungsblatt 10, Aufgabe 3': '对消息 RNVS 做 7-bit/8-bit ASCII 编码，计算 Hamming 距离、二维奇偶校验、Internet Checksum 和 CRC。',
    'Uebungsblatt 10, Aufgabe 4': '给定生成多项式 `G=x^3+1`，计算消息的 CRC 校验位，并验证接收比特串是否正确。',
    'Uebungsblatt 10, Aufgabe 5': '分析 CSMA/CD 中两主机同时发送导致碰撞、jam signal、退避和重新发送的时间关系与公式。',
    'Uebungsblatt 10, Aufgabe 6': '解释以太网最小帧长对碰撞检测的意义，并根据速率、距离和传播速度计算最小消息长度。',
}

assignment_topics = {
    'Uebungsblatt 00, Aufgabe 1': ['03_进制_图树与前缀地址.md'],
    'Uebungsblatt 00, Aufgabe 2': ['03_进制_图树与前缀地址.md'],
    'Uebungsblatt 00, Aufgabe 3': ['03_进制_图树与前缀地址.md'],
    'Uebungsblatt 01, Aufgabe 1': ['01_基础_协议与分布式系统.md'],
    'Uebungsblatt 01, Aufgabe 2': ['03_进制_图树与前缀地址.md'],
    'Uebungsblatt 01, Aufgabe 3': ['03_进制_图树与前缀地址.md'],
    'Uebungsblatt 01, Aufgabe 4': ['03_进制_图树与前缀地址.md'],
    'Uebungsblatt 02, Aufgabe 1': ['01_基础_协议与分布式系统.md', '02_分层模型_OSI_Internet_PDU与接口.md'],
    'Uebungsblatt 02, Aufgabe 2': ['01_基础_协议与分布式系统.md'],
    'Uebungsblatt 02, Aufgabe 3': ['04_传输层_UDP_TCP_序列号与滑动窗口.md'],
    'Uebungsblatt 02, Aufgabe 4': ['10_DNS_HTTP_应用层协议与Wireshark.md'],
    'Uebungsblatt 03, Aufgabe 1': ['02_分层模型_OSI_Internet_PDU与接口.md'],
    'Uebungsblatt 03, Aufgabe 2': ['02_分层模型_OSI_Internet_PDU与接口.md', '04_传输层_UDP_TCP_序列号与滑动窗口.md'],
    'Uebungsblatt 03, Aufgabe 3': ['02_分层模型_OSI_Internet_PDU与接口.md'],
    'Uebungsblatt 03, Aufgabe 4': ['02_分层模型_OSI_Internet_PDU与接口.md'],
    'Uebungsblatt 04, Aufgabe 1': ['01_基础_协议与分布式系统.md'],
    'Uebungsblatt 04, Aufgabe 2': ['01_基础_协议与分布式系统.md', '04_传输层_UDP_TCP_序列号与滑动窗口.md'],
    'Uebungsblatt 04, Aufgabe 3': ['04_传输层_UDP_TCP_序列号与滑动窗口.md'],
    'Uebungsblatt 05, Aufgabe 1': ['04_传输层_UDP_TCP_序列号与滑动窗口.md'],
    'Uebungsblatt 05, Aufgabe 2': ['04_传输层_UDP_TCP_序列号与滑动窗口.md'],
    'Uebungsblatt 05, Aufgabe 3': ['04_传输层_UDP_TCP_序列号与滑动窗口.md'],
    'Uebungsblatt 05, Aufgabe 4': ['04_传输层_UDP_TCP_序列号与滑动窗口.md', '10_DNS_HTTP_应用层协议与Wireshark.md'],
    'Uebungsblatt 06, Aufgabe 1': ['04_传输层_UDP_TCP_序列号与滑动窗口.md'],
    'Uebungsblatt 06, Aufgabe 2': ['05_TCP流量控制_拥塞控制_Reno与SACK.md'],
    'Uebungsblatt 06, Aufgabe 3': ['05_TCP流量控制_拥塞控制_Reno与SACK.md'],
    'Uebungsblatt 06, Aufgabe 4': ['05_TCP流量控制_拥塞控制_Reno与SACK.md'],
    'Uebungsblatt 07, Aufgabe 1': ['06_延迟_分组交换与电路交换.md'],
    'Uebungsblatt 07, Aufgabe 2': ['06_延迟_分组交换与电路交换.md'],
    'Uebungsblatt 07, Aufgabe 3': ['04_传输层_UDP_TCP_序列号与滑动窗口.md', '06_延迟_分组交换与电路交换.md'],
    'Uebungsblatt 07, Aufgabe 4': ['05_TCP流量控制_拥塞控制_Reno与SACK.md', '06_延迟_分组交换与电路交换.md'],
    'Uebungsblatt 08, Aufgabe 1': ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md'],
    'Uebungsblatt 08, Aufgabe 2': ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md'],
    'Uebungsblatt 08, Aufgabe 3': ['08_路由_Dijkstra_距离向量_BGP与自治系统.md'],
    'Uebungsblatt 08, Aufgabe 4': ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md'],
    'Uebungsblatt 08, Aufgabe 5': ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md'],
    'Uebungsblatt 09, Aufgabe 1': ['08_路由_Dijkstra_距离向量_BGP与自治系统.md'],
    'Uebungsblatt 09, Aufgabe 2': ['08_路由_Dijkstra_距离向量_BGP与自治系统.md'],
    'Uebungsblatt 09, Aufgabe 3': ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md', '08_路由_Dijkstra_距离向量_BGP与自治系统.md'],
    'Uebungsblatt 09, Aufgabe 4': ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md'],
    'Uebungsblatt 09, Aufgabe 5': ['09_分片_MTU_IPv4与IPv6.md'],
    'Uebungsblatt 09, Aufgabe 6': ['08_路由_Dijkstra_距离向量_BGP与自治系统.md'],
    'Uebungsblatt 10, Aufgabe 1': ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md', '12_以太网_Hub_Switch_CSMACD与物理层.md'],
    'Uebungsblatt 10, Aufgabe 2': ['12_以太网_Hub_Switch_CSMACD与物理层.md'],
    'Uebungsblatt 10, Aufgabe 3': ['11_差错检测_校验和_CRC与汉明距离.md'],
    'Uebungsblatt 10, Aufgabe 4': ['11_差错检测_校验和_CRC与汉明距离.md'],
    'Uebungsblatt 10, Aufgabe 5': ['12_以太网_Hub_Switch_CSMACD与物理层.md'],
    'Uebungsblatt 10, Aufgabe 6': ['12_以太网_Hub_Switch_CSMACD与物理层.md'],
}

exam_topics = {}
def map_exam(year, nums, topic_names):
    for num in nums:
        exam_topics[f'Klausur {year}, Frage {num}'] = topic_names

map_exam(2017, [1, 2, 11, 12, 13], ['02_分层模型_OSI_Internet_PDU与接口.md'])
map_exam(2017, [3, 7], ['01_基础_协议与分布式系统.md', '04_传输层_UDP_TCP_序列号与滑动窗口.md'])
map_exam(2017, [4], ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md', '09_分片_MTU_IPv4与IPv6.md'])
map_exam(2017, [5, 10, 16, 17, 18, 19, 20], ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md'])
map_exam(2017, [6], ['08_路由_Dijkstra_距离向量_BGP与自治系统.md'])
map_exam(2017, [8, 9, 25, 26], ['05_TCP流量控制_拥塞控制_Reno与SACK.md'])
map_exam(2017, [14, 15], ['10_DNS_HTTP_应用层协议与Wireshark.md'])
map_exam(2017, [21, 22], ['09_分片_MTU_IPv4与IPv6.md'])
map_exam(2017, [32], ['11_差错检测_校验和_CRC与汉明距离.md'])

map_exam(2018, [1, 11, 12, 13], ['02_分层模型_OSI_Internet_PDU与接口.md'])
map_exam(2018, [2], ['01_基础_协议与分布式系统.md'])
map_exam(2018, [3, 30, 31, 32], ['12_以太网_Hub_Switch_CSMACD与物理层.md'])
map_exam(2018, [4, 5, 33], ['11_差错检测_校验和_CRC与汉明距离.md'])
map_exam(2018, [6, 8, 21, 22, 23, 24], ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md'])
map_exam(2018, [7, 25, 28, 29], ['12_以太网_Hub_Switch_CSMACD与物理层.md'])
map_exam(2018, [9, 26, 27], ['05_TCP流量控制_拥塞控制_Reno与SACK.md'])
map_exam(2018, [10, 14, 15], ['10_DNS_HTTP_应用层协议与Wireshark.md'])
map_exam(2018, [16, 17, 18], ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md', '10_DNS_HTTP_应用层协议与Wireshark.md'])
map_exam(2018, [19, 20], ['09_分片_MTU_IPv4与IPv6.md'])

map_exam(2019, [1, 2], ['06_延迟_分组交换与电路交换.md'])
map_exam(2019, [3, 15, 16, 17, 18], ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md'])
map_exam(2019, [4], ['08_路由_Dijkstra_距离向量_BGP与自治系统.md'])
map_exam(2019, [5, 7, 8, 9, 10, 11], ['10_DNS_HTTP_应用层协议与Wireshark.md'])
map_exam(2019, [6], ['04_传输层_UDP_TCP_序列号与滑动窗口.md'])
map_exam(2019, [12], ['07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md', '10_DNS_HTTP_应用层协议与Wireshark.md'])
map_exam(2019, [13, 14], ['09_分片_MTU_IPv4与IPv6.md'])
map_exam(2019, [19, 20, 21], ['05_TCP流量控制_拥塞控制_Reno与SACK.md'])
map_exam(2019, [22, 23, 24], ['12_以太网_Hub_Switch_CSMACD与物理层.md'])
map_exam(2019, [25], ['11_差错检测_校验和_CRC与汉明距离.md'])

def strip_markdown_for_translation(line):
    line = line.strip()
    line = re.sub(r'^\s*[-*+]\s*', '', line)
    line = re.sub(r'^[✅☑☒☐▪•]+\s*', '', line)
    line = line.replace('*', '').replace('_', '').replace('`', '')
    line = re.sub(r'\s{2,}$', '', line)
    return line.strip()

def extract_chinese_question_text(text):
    """Use the existing bilingual exam question text when it is present."""
    part = re.split(r'(?im)^\s*\*\*(?:L枚sung|Lösung|参考答案|答案|解析|瑙ｆ瀽|鍙傝€冪瓟妗).*', text, maxsplit=1)[0]
    zh_lines = []
    for raw_line in part.splitlines():
        if re.match(r'^\s*#+\s+', raw_line):
            continue
        line = strip_markdown_for_translation(raw_line)
        if not line or line in {'---'}:
            continue
        cjk_count = len(re.findall(r'[\u4e00-\u9fff]', line))
        if cjk_count >= 2:
            zh_lines.append(line)
    cleaned = []
    for line in zh_lines:
        if cleaned and cleaned[-1] == line:
            continue
        cleaned.append(line)
    return '\n'.join(cleaned).strip()

def chinese_hint(block, topic_title):
    if block['source'] in assignment_zh:
        return assignment_zh[block['source']]
    extracted = extract_chinese_question_text(block.get('question', block.get('match', block['text'])))
    if extracted:
        return extracted
    title = block['title']
    text = block.get('match', block['text'])
    hay = (title + '\n' + text).lower()
    if 'welche' in hay or 'multiple choice' in hay:
        return f'本题是关于“{topic_title}”的判断/选择题：需要根据德文题干和选项判断哪些陈述正确，并结合后面的答案解析复习相关概念。'
    if 'berechnen' in hay or 'calculate' in hay:
        return f'本题是关于“{topic_title}”的计算题：需要从题干给出的参数出发，按公式或协议规则计算结果，并查看后续解答核对步骤。'
    if 'zeichnen' in hay:
        return f'本题要求围绕“{topic_title}”画图或补图，例如时序图、分片图、路由图或协议结构图；后续解答给出画法和关键标注。'
    if 'nennen' in hay or 'geben sie' in hay:
        return f'本题要求列举或说明“{topic_title}”中的概念、协议、字段或区别；后续解答给出要点。'
    return f'本题围绕“{topic_title}”展开；请先阅读德文原题，再结合下方解答理解题意和考点。'

topics = [
    {'name':'01_基础_协议与分布式系统.md','title':'基础、协议与分布式系统','keys':['Pizzadienst','Rechnernetze und verteilte Systeme','verteilte Systeme','verbindungslose','verbindungsorientierte','Verbindungsaufbau','Handshake','Grundlagen']},
    {'name':'02_分层模型_OSI_Internet_PDU与接口.md','title':'分层模型：OSI、Internet、PDU 与接口','keys':['Protokollschichtung','OSI-Modell','OSI-Referenzmodell','Schichtenmodell','PDU','SDU','PCI','Schnittbildung','Internetmodell','Peer-Entity','Dienstschnitt','Protokollschnitt','Systemschnitt']},
    {'name':'03_进制_图树与前缀地址.md','title':'进制、图树与前缀地址','keys':['Stellenwertsystem','Zahlensystem','Baeum','B盲um','Binaer','Bin盲r','Praefix','Pr盲fix','vollvermascht','Anforderungen des Internets']},
    {'name':'04_传输层_UDP_TCP_序列号与滑动窗口.md','title':'传输层：UDP、TCP、序列号与滑动窗口','keys':['RFC 768','UDP','Transportschicht','sendto','TCP Sequenznummern','Sequenznummern','Sendefenster','Sliding-Window','Socket']},
    {'name':'05_TCP流量控制_拥塞控制_Reno与SACK.md','title':'TCP 流量控制、拥塞控制、Reno 与 SACK','keys':['Flusssteuerung','Staukontrolle','Ueberlast','Überlast','Reno','Tahoe','Slow Start','SACK','selektive','Durchsatz','Verlustrate','CongWin','Congestion']},
    {'name':'06_延迟_分组交换与电路交换.md','title':'延迟、分组交换与电路交换','keys':['Verzoeger','Verz枚ger','Verzöger','Paketvermittlung','Leitungsvermittlung','Signalver','Warteschlangen','Traffic Intensity','Intensit']},
    {'name':'07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md','title':'IP 地址：CIDR、IPv4/IPv6、NAT 与 ARP','keys':['Addressierung in Rechnernetzen','Adressierung /','CIDR','Subnetz','Subnetting','IPv4','IPv6','NAT','ARP','Private IP','RFC1918','IP-Header','Default Route','Network Address Translation','Zusammenspiel von IPv4 und ARP']},
    {'name':'08_路由_Dijkstra_距离向量_BGP与自治系统.md','title':'路由：Dijkstra、距离向量、BGP 与自治系统','keys':['Routing','Dijkstra','Link-State','Distanzvektor','Distance-Vector','Count to Infinity','Autonome Systeme','BGP','OSPF','Peering','Transit','Wegewahl']},
    {'name':'09_分片_MTU_IPv4与IPv6.md','title':'分片：MTU、IPv4 与 IPv6','keys':['Fragmentierung','MTU','Packet too big','ICMPv6','Fragment Offset']},
    {'name':'10_DNS_HTTP_应用层协议与Wireshark.md','title':'DNS、HTTP、应用层协议与 Wireshark','keys':['DNS','HTTP','IMAP','SMTP','FTP','Anwendungsschicht','Wireshark','ping','traceroute','SSHv2']},
    {'name':'11_差错检测_校验和_CRC与汉明距离.md','title':'差错检测：校验和、CRC 与汉明距离','keys':['Fehlererkennung','Fehlerkorrektur','CRC','Checksum','Checksumme','Parit','Hamming','Block Check Character','BCC','UDP错误检测']},
    {'name':'12_以太网_Hub_Switch_CSMACD与物理层.md','title':'以太网、Hub/Switch、CSMA/CD 与物理层','keys':['Ethernet','CSMA','Hub','Switch','MAC','Rahmen','Mindestrahmen','Bituebertragung','Bit眉bertragung','Lichtwellenleiter','elektrische Leiter','Signalen','Bitstroemen','Bitströmen']},
    {'name':'13_ADSL_ISDN_PPP_协议栈与传输介质.md','title':'ADSL、ISDN、PPP、协议栈与传输介质','keys':['ADSL','ISDN','PPP','ATM','Splitter','Protokollstapel','Lichtwellenleiter','elektrische Leiter','Medien']},
]
summary = {
    t['name']: f"- {t['title']} 的题目汇总。\n- 下方每道题均包含完整题块摘录和来源说明。" for t in topics
}
summary.update({
'01_基础_协议与分布式系统.md':'- 协议规定消息格式、顺序、语义和错误处理。\n- 分布式系统关注多个节点协同以及通信不可靠带来的状态一致性问题。\n- 无连接通信每个包独立发送；面向连接通信先建立逻辑状态。',
'02_分层模型_OSI_Internet_PDU与接口.md':'- OSI 七层与 Internet 模型是高频基础。\n- PDU = SDU + PCI；发送向下封装，接收向上解封装。\n- Dienstschnitt、Protokollschnitt、Systemschnitt 要结合图判断。',
'03_进制_图树与前缀地址.md':'- 进制转换、完全图边数、树高度和前缀子树是地址学习的基础。\n- 前缀越长，范围越小；前缀越短，聚合范围越大。',
'04_传输层_UDP_TCP_序列号与滑动窗口.md':'- UDP 无连接、报文式、无可靠性保证。\n- TCP 面向连接、字节流、用 Seq/ACK/窗口实现可靠传输。\n- ACK 号表示下一个期望字节。',
'05_TCP流量控制_拥塞控制_Reno与SACK.md':'- 流量控制保护接收方，拥塞控制保护网络。\n- Slow Start 指数增长，Congestion Avoidance 线性增长。\n- Reno 用 Fast Retransmit/Fast Recovery 处理三重复 ACK。',
'06_延迟_分组交换与电路交换.md':'- 主要延迟：处理、排队、发送、传播。\n- 发送延迟 L/R；传播延迟 d/v。\n- Traffic intensity I=aL/R，用于判断队列是否稳定。',
'07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md':'- CIDR 用可变前缀进行地址分配和聚合。\n- ARP 解析下一跳 MAC；NAT/NAPT 改写地址/端口。\n- IPv6 地址 128 bit，路由器不做分片。',
'08_路由_Dijkstra_距离向量_BGP与自治系统.md':'- Link-State 用 Dijkstra/SPF；Distance-Vector 用 Bellman-Ford 思想。\n- BGP 是 AS 之间的 Path-Vector 协议。\n- Count-to-Infinity 是距离向量经典问题。',
'09_分片_MTU_IPv4与IPv6.md':'- MTU 限制每条链路最大 IP 包大小。\n- IPv4 路由器可分片，offset 单位是 8 字节。\n- IPv6 中间路由器不分片，而返回 Packet Too Big。',
'10_DNS_HTTP_应用层协议与Wireshark.md':'- DNS 是应用层协议，常用 UDP 53，也可用 TCP 53。\n- HTTP/IMAP/SSH 等是应用层协议。\n- Wireshark 题重点是识别协议层次和报文序列。',
'11_差错检测_校验和_CRC与汉明距离.md':'- 奇偶校验、BCC、Internet Checksum、CRC 都是错误检测机制。\n- Hamming 距离决定检测/纠错能力。\n- CRC 使用模 2 除法，生成多项式最高次数决定校验位数。',
'12_以太网_Hub_Switch_CSMACD与物理层.md':'- Hub 广播到所有端口；Switch 学习 MAC 后定向转发。\n- CSMA/CD 需要最小帧长以保证碰撞可检测。\n- 物理层关注信号、介质、带宽、编码和传播。',
'13_ADSL_ISDN_PPP_协议栈与传输介质.md':'- ADSL/ISDN/PPP/ATM 题常考协议栈封装和单位换算。\n- PPP 的 LCP 管链路，NCP 管网络层协议配置。'
})

blocks = []
solution_map = {}
solution_dir = base / '作业解答'
if solution_dir.exists():
    for sol_path in sorted(solution_dir.glob('blatt-*_中德对照解答.md')):
        m = re.search(r'blatt-(\d+)_', sol_path.name)
        if not m:
            continue
        blatt_no = m.group(1)
        sol_raw = sol_path.read_text(encoding='utf-8', errors='replace')
        sol_heads = list(re.finditer(r'(?m)^##\s+(\d+)\.\s+', sol_raw))
        for i, head in enumerate(sol_heads):
            start = head.start()
            end = sol_heads[i+1].start() if i+1 < len(sol_heads) else len(sol_raw)
            task_no = head.group(1)
            solution_map[(blatt_no, task_no)] = sol_raw[start:end].strip()

assign_dir = base / '.codex-tmp' / 'assignment_text'
assignment_prefixes = {
    '00': [('1', 'Anforderungen'), ('2', 'Das Stellenwertsystem'), ('3', 'Rechnen')],
    '01': [('1', 'Grundlagen'), ('2', 'Von Netzen'), ('3', 'Adressierung'), ('4', 'Zahlen')],
    '02': [('1', 'Der Pizzadienst'), ('2', 'Rechnernetze'), ('3', 'RFC 768'), ('4', 'Einf')],
    '03': [('1', 'Protokollschichtung'), ('2', 'Datenpakete'), ('3', 'Bestandteile'), ('4', 'ISO')],
    '04': [('1', 'Verbindungslose'), ('2', 'Verbindungsaufbau'), ('3', 'Sequenznummern')],
    '05': [('1', 'Transportschicht'), ('2', 'TCP Sequenznummern'), ('3', 'Sequenznummern'), ('4', '3-Way')],
    '06': [('1', 'TCP-Verbindung'), ('2', 'Selektive'), ('3', 'TCP Reno'), ('4', 'Approximierung')],
    '07': [('1', 'Verz'), ('2', 'Paket-'), ('3', 'Fenster'), ('4', 'Staukontrolle')],
    '08': [('1', 'Addressierung'), ('2', 'Hierarchische'), ('3', 'Link-State'), ('4', 'Private'), ('5', 'Network')],
    '09': [('1', 'Distanz'), ('2', 'Autonome'), ('3', 'Wegewahl'), ('4', 'IPv6'), ('5', 'Fragmentierung'), ('6', 'Count')],
    '10': [('1', 'Zusammenspiel'), ('2', 'Was ist'), ('3', 'Fehlererkennung'), ('4', 'CRC'), ('5', 'CSMA'), ('6', 'Ethernet')],
}
if assign_dir.exists():
    for p in sorted(assign_dir.glob('blatt-*_uebungsblatt.txt')):
        m = re.search(r'blatt-(\d+)_', p.name)
        if not m:
            continue
        blatt = m.group(1)
        raw = clean_lines(p.read_text(encoding='utf-8', errors='replace'))
        starts = []
        for task_no, prefix in assignment_prefixes.get(blatt, []):
            pat = re.compile(r'(?m)^' + re.escape(task_no) + r'\.\s+' + re.escape(prefix))
            mat = pat.search(raw)
            if mat:
                starts.append((task_no, mat.start()))
        starts.sort(key=lambda item: item[1])
        for i, (expected_no, start) in enumerate(starts):
            end = starts[i+1][1] if i+1 < len(starts) else len(raw)
            text = raw[start:end].strip()
            mm = re.match(r'(\d+)\.', text)
            if not mm:
                continue
            title = text.splitlines()[0]
            task_no = mm.group(1)
            answer = solution_map.get((blatt, task_no), '')
            if answer:
                answer_text = answer
            else:
                answer_text = '未在作业解答目录中匹配到对应解答。'
            blocks.append({
                'source': f'Uebungsblatt {blatt}, Aufgabe {task_no}',
                'title': title,
                'question': text,
                'answer': answer_text,
                'text': text + '\n\n' + answer_text,
                'match': text,
                'kind': '作业',
            })

for p in [base/'考卷'/'2015.md', base/'考卷'/'2017.md', base/'考卷'/'2018.md', base/'考卷'/'2019.md']:
    if not p.exists():
        continue
    year = re.search(r'(20\d\d)', p.name).group(1)
    raw = p.read_text(encoding='utf-8', errors='replace')
    heads = list(re.finditer(r'(?m)^(###\s+[^\n]+|##\s+(?!I\.|II\.|III\.|IV\.|V\.|VI\.|VII\.|VIII\.|IX\.|X\.|1 |2 |3 |4 |5 |6 |7 |8 |9 )[A-Za-zÄÖÜäöü].*)', raw))
    for i, head in enumerate(heads):
        start = head.start()
        end = heads[i+1].start() if i+1 < len(heads) else len(raw)
        body = raw[start:end].strip()
        if len(body) < 40:
            continue
        title = body.splitlines()[0].lstrip('#').strip()
        q = re.search(r'Frage\s*([0-9]+)', title, flags=re.I)
        source = f'Klausur {year}, Frage {q.group(1)}' if q else f'Klausur {year}, Abschnitt {title}'
        cut = re.search(r'(?im)^\*\*(?:L枚sung|Lösung|参考答案|答案|解析|瑙ｆ瀽|鍙傝€冪瓟妗).*', body)
        question_for_match = body[:cut.start()].strip() if cut else body
        if len(question_for_match) < 40:
            question_for_match = body
        answer_text = body[cut.start():].strip() if cut else body
        blocks.append({
            'source': source,
            'title': title,
            'question': question_for_match,
            'answer': answer_text,
            'text': body,
            'match': question_for_match,
            'kind': '考卷',
        })

assigned = {topic['name']: [] for topic in topics}
for block in blocks:
    if block['kind'] == '作业':
        for mapped_topic in assignment_topics.get(block['source'], []):
            assigned[mapped_topic].append(block)
        continue
    if block['source'] in exam_topics:
        for mapped_topic in exam_topics[block['source']]:
            assigned[mapped_topic].append(block)
        continue
    hay = block['title'] + '\n' + block.get('match', block['text'])
    for topic in topics:
        if any(keyword_hit(hay, key) for key in topic['keys']):
            assigned[topic['name']].append(block)

for topic_name, hint in [
    ('06_延迟_分组交换与电路交换.md','Klausur 2019, Frage 1'),
    ('06_延迟_分组交换与电路交换.md','Klausur 2019, Frage 2'),
    ('07_IP地址_CIDR_IPv4_IPv6_NAT与ARP.md','Klausur 2019, Frage 3'),
    ('08_路由_Dijkstra_距离向量_BGP与自治系统.md','Klausur 2019, Frage 4'),
    ('04_传输层_UDP_TCP_序列号与滑动窗口.md','Klausur 2019, Frage 6'),
    ('10_DNS_HTTP_应用层协议与Wireshark.md','Klausur 2019, Frage 7'),
]:
    for block in blocks:
        if hint == block['source'] and block not in assigned[topic_name]:
            assigned[topic_name].append(block)

for topic in topics:
    topic_name = topic['name']
    collected = []
    by_digest = {}
    for block in assigned[topic_name]:
        digest = hashlib.sha1(norm(block['text']).encode('utf-8', 'ignore')).hexdigest()
        if digest not in by_digest:
            by_digest[digest] = {'block': block, 'sources': []}
            collected.append(by_digest[digest])
        if block['source'] not in by_digest[digest]['sources']:
            by_digest[digest]['sources'].append(block['source'])
    lines = [f"# {topic['title']}\n", '## 知识点总结\n', summary[topic_name] + '\n', '## 完整题目与解答汇总\n']
    if not collected:
        lines.append('暂无自动匹配到的题目。\n')
    for idx, item in enumerate(collected, 1):
        block = item['block']
        sources = item['sources']
        lines.append(f"### 题目 {idx}: {block['title']}\n")
        lines.append(f"**类型：** {block['kind']}  \n")
        lines.append('**来源说明：** ' + '；'.join(sources) + '  \n')
        if len(sources) > 1:
            lines.append('**重复处理：** 完全重复题目已合并保留一份。  \n')
        lines.append('\n#### 题目中文翻译 / 中文题意\n')
        lines.append(chinese_hint(block, topic['title']) + '\n')
        lines.append('#### 德文原题\n')
        lines.append('```text\n' + block.get('question', block['match']).strip() + '\n```\n')
        lines.append('#### 解答\n')
        lines.append(sanitize_markdown_fragment(block.get('answer', block['text'])) + '\n')
        lines.append('**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。\n')
        lines.append('\n---\n')
    (outdir / topic_name).write_text('\n'.join(lines), encoding='utf-8-sig')

readme = ['# 分知识点汇总\n', '本目录按知识点汇总作业与考卷题目。每个文件内包含知识点总结、完整题目与解答摘录、来源说明；完全重复题目在同一知识点内合并。\n', '## 文件列表\n']
for topic in topics:
    readme.append(f"- [{topic['title']}]({topic['name']})")
(outdir / 'README.md').write_text('\n'.join(readme) + '\n', encoding='utf-8-sig')
print('generated', len(topics), 'topic files')
for topic in topics:
    print(topic['name'], len(assigned[topic['name']]))
