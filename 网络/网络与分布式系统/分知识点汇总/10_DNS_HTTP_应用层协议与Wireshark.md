# DNS、HTTP、应用层协议与 Wireshark

## 知识点总结

- DNS 是应用层协议，常用 UDP 53，也可用 TCP 53。
- HTTP/IMAP/SSH 等是应用层协议。
- Wireshark 题重点是识别协议层次和报文序列。

## 完整题目与解答汇总

### 题目 1: 4. Einführung in textbasiertes Arbeiten mit Linux

**类型：** 作业  

**来源说明：** Uebungsblatt 02, Aufgabe 4  


#### 题目中文翻译 / 中文题意

练习 Linux 命令行、ping 和 traceroute：查看目录、理解 man page、测量往返时延并解释 traceroute 输出。

#### 德文原题

```text
4. Einführung in textbasiertes Arbeiten mit Linux
Linux bietet eine Vielzahl an Programmen, die der praktischen Veranschaulichung der in der Vorlesung
vermittelten Inhalte dienen. Im Rahmen dieser Aufgaben lernen Sie grundlegende Tools kennen.
Falls noch nicht geschehen, machen Sie sich daher mit dem grundlegende Umgang der Kommandozeile
unter Linux vertraut.
(a) Melden Sie sich zunächst mit Ihrer Benutzerkennung und Ihrem Passwort an einem Rechner des
CIP-Pools an und öffnen Sie eine Konsole. Sollte Sie die Aufgaben von außerhalb des Universitäts-
gebäudes erledigen wollen, können Sie sich alternativ via SSH von einem beliebigen Rechner Ihrer
Wahl einloggen. Voraussetzung ist ein installierter SSH Client. Eine ausführliche Anleitung bietet
Ihnen die Seite der Rechnerbetriebsgruppe1.
i. Ermitteln Sie den absoluten Pfad Ihres Home-Verzeichnisses und zeigen Sie dessen Inhalt an!
ii. Wechseln Sie in das Wurzelverzeichnis und dann zurück in Ihr Home-Verzeichnis!
iii. Was ist eine „man-Page”? Hinweis: Benutzen Sie den Befehl man man!
iv. Mit welchem Parameter zeigt ls auch versteckte Dateien an? Hinweis: man-Page: [ls(1)]!
(b) Der ping-Befehl schickt Anfragen zu dem per Hostname oder IP-Adresse spezifizierten Zielrechner,
um festzustellen ob der Zielrechner erreichbar ist. Mit dem Erhalt einer Antwort zeigt ping die
RTD (roundtrip delay) an. Beachten Sie die man-page des Befehls: (man ping).
i. Was versteht man unter roundtrip delay?
ii. Versuchen Sie den Host „www.nm.ifi.lmu.de” mit dem Programm ping zu erreichen! Dabei
sollen 10 Anfragen im Abstand von 2 Sekunden und je 100 Bytes Nutzdaten verschickt werden.
iii. Wie sind die einzelnen Spalten in der Ausgabe des ping-Befehls zu interpretieren?
(c) Der traceroute-Befehl zeigt den Pfad von der Quelle bis zur Senke durch ein Rechnernetz und
misst die RTD zu jedem einzelnen Knoten auf diesem Pfad. Beachten Sie die man-page des Befehls:
(man traceroute).
i. Interpretieren Sie die Ausgabe von traceroute zum Zielrechner „www.nm.ifi.lmu.de”! Welche
Informationen beinhaltet die erste Zeile der Ausgabe?
ii. In den darauffolgenden Zeilen stehen je drei Werte, meist in Millisekunden angegeben. Wofür
stehen diese Werte?
iii. Die häufige Überprüfung des Pfades zu einem bestimmten Zielrechner mit traceroute zeigt
manchmal andere Einträge mit einem verschiedenen Pfad. Was kann diese Beobachtung bedeu-
ten?
1https://www.rz.ifi.lmu.de/infos/ssh_de.html
```

#### 解答

**4. Textbasiertes Arbeiten mit Linux / Linux 命令行**

**(a) Grundbefehle**

**中文说明：** 这些命令是 Linux 文本工作环境的基础：`pwd` 看当前位置，`ls` 列出文件，`cd` 切换目录，`man` 查看帮助手册。

| Aufgabe | Befehl |
|---|---|
| Home-Pfad anzeigen | `pwd` |
| Inhalt anzeigen | `ls` |
| Wurzelverzeichnis | `cd /` |
| zurueck ins Home | `cd ~` |
| Man-Page | Handbuchseite zu Befehlen, z.B. `man man` |
| versteckte Dateien mit `ls` | `ls -a` |

**(b) ping**

**DE:** Roundtrip delay (RTD) ist die Zeit vom Senden einer Anfrage bis zum Empfang der Antwort.

**中文：** 往返时延 RTD/RTT 是从发出请求到收到响应之间经过的总时间，包含去程和回程。

Beispielbefehl:

```bash
ping -c 10 -i 2 -s 100 www.nm.ifi.lmu.de
```

Typische Spalten: Anzahl Bytes, Zielhost/IP, ICMP-Sequenznummer, TTL, Zeit/RTD.

**(c) traceroute**

```bash
traceroute www.nm.ifi.lmu.de
```

**DE:** Die erste Zeile nennt Ziel, Ziel-IP, maximale Hop-Zahl und Paketgroesse. Danach zeigt jede Zeile einen Hop; die drei Zeitwerte sind Messungen fuer drei Probe-Pakete. Unterschiedliche Pfade koennen Lastverteilung, dynamisches Routing oder geaenderte Netzbedingungen bedeuten.

**中文：** 第一行通常给出目标、目标 IP、最大跳数和探测包大小。之后每一行代表一跳，三个时间值是三次探测的往返时间。多次 traceroute 路径不同，可能表示负载均衡、动态路由变化或网络状态变化。

**Wissen / 知识点：** `ping` misst Erreichbarkeit und RTT; `traceroute` nutzt TTL/Hop-Limit, um Zwischenrouter sichtbar zu machen。

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 2: 4. 3-Way-Handshake und Sequenznummern bei TCP

**类型：** 作业  

**来源说明：** Uebungsblatt 05, Aufgabe 4  


#### 题目中文翻译 / 中文题意

用 Wireshark 分析 TCP 抓包：找出三次握手、连接释放、RTD、绝对/相对序列号，并判断 SSHv2 是否使用 TCP。

#### 德文原题

```text
4. 3-Way-Handshake und Sequenznummern bei TCP
Protokollkonzepte wie 3-Way-Handshaking und Sequenznummern sind Mechanismen für das Verbin-
dungsmanagement und für die zuverlässige Kommunikation. Diese Mechanismen werden z.B. im Trans-
mission Control Protocol (TCP) eingesetzt. Zur Bearbeitung dieser Aufgabe wird die Trace-Datei
trace3.pcap bereitgestellt, die mitgeschnittenen TCP-Verkehr enthält.
Die Datei lässt sich z.B. mit dem freien Programm Wireshark1 öffnen, das den mitgeschnittenen TCP-
Verkehr grafisch aufbereitet anzeigen und filtern kann.
Wireshark stellt die PDUs verschiedener Schichten tabellarisch dar. In dieser Aufgabe soll es um TCP-
PDUs (Segmente) gehen. Sie sind daran zu erkennen, dass in der Protocol -Spalte TCP steht. Die Infor-
mationen der TCP-PCI werden von Wireshark übersichtlich aufbereitet und je Segment dargestellt.
(a) Identifizieren Sie die zum 3-Way-Handshaking Vorgang gehörenden Segmente in trace3.pcap.
(b) Identifizieren Sie die zum Verbindungsabbau gehörigen Segmente.
(c) Berechnen Sie aus den Paketen des 3-Way-Handshake das so genannte Round Trip Delay (RTD).
Das ist die Zeit, die vom Versenden eines Segments bis zum Erhalt einer Antwort vergeht.
(d) Welche absoluten und relativen TCP-Sequenznummern besitzen diese Segmente?
(e) In dem Mitschnitt werden (in der Standardkonfiguration von Wireshark) auch PDUs vom Protokoll
der Anwendungsschicht (SSHv2 ) angezeigt. Nutzt dieses Protokoll auch TCP? Begründen Sie Ihre
Vermutung kurz.
1https://www.wireshark.org/
```

#### 解答

**4. TCP in Wireshark**

**DE:** Ohne die Datei `trace3.pcap` im Arbeitsverzeichnis lassen sich konkrete Paketnummern und Zeiten nicht bestimmen. Vorgehen:

**中文：** 当前工作目录中没有 `trace3.pcap`，所以不能给出具体包号、时间戳和绝对序列号。可以按下面方法在 Wireshark 中分析：

| Teil | Vorgehen |
|---|---|
| 3-Way-Handshake | Filter `tcp.flags.syn == 1`, Segmente `SYN`, `SYN ACK`, `ACK` |
| Verbindungsabbau | Segmente mit `FIN` und abschliessenden ACKs |
| RTD | Zeit zwischen SYN und SYN+ACK oder zwischen SYN+ACK und ACK |
| absolute/relative SeqNr | In Wireshark TCP-Details; relative Nummern starten meist bei 0 |
| SSHv2 nutzt TCP? | Ja, wenn SSHv2-PDUs in TCP-Segmenten mit TCP-Port 22 oder entsprechender Verbindung gekapselt sind |

**Wissen / 知识点：** TCP 的序号按字节编号；ACK 号表示“下一个期望字节”。UDP 是报文式，TCP 是字节流式。

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 3: Welches Transportschicht mit DNS?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Welches Transportschicht mit DNS?  


#### 题目中文翻译 / 中文题意

DNS使用哪种传输层协议？

#### 德文原题

```text
### Welches Transportschicht mit DNS?

**DNS使用哪种传输层协议？**
```

#### 解答

**Lösung / 答案：**

- **UDP端口53：** 普通DNS查询（默认）
- **TCP端口53：** 区域传输、大响应（>512字节）

---

**III. wie in 2014 (HTTP, ICMP, OSPF, IMAP)**

这似乎是指2014年考过的协议分类题。

|协议|层级|说明|
|---|---|---|
|HTTP|应用层|超文本传输协议|
|ICMP|网络层|Internet控制消息协议|
|OSPF|网络层|开放最短路径优先（链路状态路由）|
|IMAP|应用层|邮件访问协议|

---

**IV. Elektrische Leiter v. Lichtwellenleiter**

**电导体 vs 光纤**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 4: Frage (b) / 问题(b)

**类型：** 考卷  

**来源说明：** Klausur 2017, Abschnitt Frage (b) / 问题(b)  


#### 题目中文翻译 / 中文题意

该图显示的是接口。
图示：显示HTTP、TCP、IP、Ethernet、(WAN)等协议层

#### 德文原题

```text
### Frage (b) / 问题(b)

**Diese Abbildung zeigt den _____ Schnitt.**  
**该图显示的是_____接口。**

图示：显示HTTP、TCP、IP、Ethernet、(WAN)等协议层
```

#### 解答

**参考答案 / Lösung:** **Protokollschnitt / 协议接口**

**Begründung / 理由:**

- 显示了对等实体之间的通信
- 同层协议之间的逻辑通信
- Kommunikation zwischen Peer-Entities
- Logische Kommunikation zwischen gleichrangigen Protokollen

---

**4 Domain Name System (6 Punkte)**

**4 域名系统（6分）**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 5: Frage 14 / 第14题

**类型：** 考卷  

**来源说明：** Klausur 2017, Frage 14  


#### 题目中文翻译 / 中文题意

本题围绕“DNS、HTTP、应用层协议与 Wireshark”展开；请先阅读德文原题，再结合下方解答理解题意和考点。

#### 德文原题

```text
### Frage 14 / 第14题

**(a) Wie viele Anfragen zur oben gezeigten DNS-Anfrage waren nötig, um den Hostnamen [www.ifi.lmu.de](http://www.ifi.lmu.de/) aufzulösen?**
```

#### 解答

**解析主机名 [www.ifi.lmu.de](http://www.ifi.lmu.de/) 需要多少次DNS查询？**

**参考答案 / Lösung:** **4**

分析查询过程：

1. 查询根服务器 → 获取 .de 的NS记录
2. 查询 .de 服务器 → 获取 lmu.de 的NS记录
3. 查询 lmu.de 服务器 → 获取 [www.ifi.lmu.de](http://www.ifi.lmu.de/) 的CNAME
4. 查询获取最终A记录

---

**(b) Die URL [http://www.ifi.lmu.de/](http://www.ifi.lmu.de/) soll in einem Web-Browser angezeigt werden. Geben Sie die IPv4-Adresse des Rechners an, an den die HTTP-Anfrage gestellt wird!**  
**URL [http://www.ifi.lmu.de/](http://www.ifi.lmu.de/) 要在浏览器中显示。给出HTTP请求发送到的IPv4地址！**

**参考答案 / Lösung:** **141.84.94.49**

从DNS输出中可以看到：

- [www.ifi.lmu.de](http://www.ifi.lmu.de/) → CNAME → salerno.tcs.ifi.lmu.de
- salerno.tcs.ifi.lmu.de → A → 141.84.94.49

---

**(c) Ist die DNS-Anfrage rekursiv oder iterativ? Begründen Sie Ihre Antwort!**  
**DNS查询是递归的还是迭代的？请说明理由！**

**参考答案 / Lösung:** **Iterativ / 迭代**

**Begründung / 理由:**

- 客户端自己发送请求，不依赖DNS服务器代为查询
- 可以看到客户端分别向根服务器、.de服务器、lmu.de服务器发送查询
- Der Client würde man sehen, dass er selbst die Requests sendet und nicht der DNS-Server für ihn

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 6: Frage 15 / 第15题

**类型：** 考卷  

**来源说明：** Klausur 2017, Frage 15  


#### 题目中文翻译 / 中文题意

在ISO/OSI参考模型的第四层，DNS使用哪种传输协议...
...用于区域传输？

#### 德文原题

```text
### Frage 15 / 第15题

**Welches Transportprotokoll wird auf Schicht IV des ISO/OSI-Referenzmodells bei DNS...**  
**在ISO/OSI参考模型的第四层，DNS使用哪种传输协议...**

**(a) ...für Zonentransfers eingesetzt?**  
**...用于区域传输？**
```

#### 解答

**参考答案 / Lösung:** **TCP**

区域传输数据量大，需要可靠传输。

**(b) ...für DNS-Anfragen empfohlen?**  
**...推荐用于DNS查询？**

**参考答案 / Lösung:** **UDP**

普通DNS查询数据量小，UDP更高效。

---

**5 Adressierung in Rechnernetzen**

**5 计算机网络中的寻址**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 7: Frage 10 / 第10题

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 10  


#### 题目中文翻译 / 中文题意

以下关于DNS的哪些陈述是正确的？
DNS将给定的主机名映射到最多一个IP地址。
DNS是应用层协议。
DNS将IP地址映射到MAC地址。

#### 德文原题

```text
### Frage 10 / 第10题

**Welche Aussagen treffen auf DNS zu?**  
**以下关于DNS的哪些陈述是正确的？**

- ○ DNS bildet einen gegebenen Hostnamen auf höchstens eine IP-Adresse ab.
    - DNS将给定的主机名映射到最多一个IP地址。
- ☒ DNS ist ein Protokoll der Anwendungsschicht.
    - DNS是应用层协议。
- ○ DNS bildet IP-Adressen auf MAC-Adressen ab.
    - DNS将IP地址映射到MAC地址。
```

#### 解答

**解析：**

- ✗ 第一项错误：一个域名可以对应多个IP地址（负载均衡）
- ✓ 第二项正确：DNS工作在应用层（第7层）
- ✗ 第三项错误：IP到MAC的映射是ARP的功能，不是DNS

---

**II. ISO OSI-Schichtmodell / ISO OSI层模型**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 8: Frage 14 / 第14题

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 14  


#### 题目中文翻译 / 中文题意

本题围绕“DNS、HTTP、应用层协议与 Wireshark”展开；请先阅读德文原题，再结合下方解答理解题意和考点。

#### 德文原题

```text
### Frage 14 / 第14题

**(a) Wie viele DNS-Anfragen waren nötig, um den Hostnamen [www.ifi.lmu.de](http://www.ifi.lmu.de/) wie oben abgebildet aufzulösen?**
```

#### 解答

**解析主机名 [www.ifi.lmu.de](http://www.ifi.lmu.de/) 需要多少次DNS查询？**

**Lösung / 答案：** **4**

**分析查询过程：**

1. 查询根服务器（.）→ 获取 .de 的NS记录
2. 查询 .de 服务器 → 获取 lmu.de 的NS记录
3. 查询 lmu.de 服务器 → 获取 ifi.lmu.de 的NS记录
4. 查询 ifi.lmu.de 服务器 → 获取 [www.ifi.lmu.de](http://www.ifi.lmu.de/) 的记录

---

**(b) Die URL [http://www.ifi.lmu.de/](http://www.ifi.lmu.de/) soll in einem Web-Browser angezeigt werden. Geben Sie die IPv4-Adresse des Rechners an, an den die HTTP-Anfrage gestellt wird.**  
**URL [http://www.ifi.lmu.de/](http://www.ifi.lmu.de/) 要在浏览器中显示。给出HTTP请求发送到的IPv4地址。**

**Lösung / 答案：** **141.84.94.144**

从DNS输出中可以看到：

- [www.ifi.lmu.de](http://www.ifi.lmu.de/) → CNAME → stellenbosch.tcs.ifi.lmu.de
- stellenbosch.tcs.ifi.lmu.de → A → 141.84.94.144

---

**(c) Ist die oben angeführte DNS-Anfrage rekursiv, iterativ oder hybrid? Warum?**  
**上述DNS查询是递归的、迭代的还是混合的？为什么？**

**Lösung / 答案：** **Iterativ / 迭代**

**Begründung / 理由：**

- 使用了 `+trace` 选项，客户端自己逐级查询
- 可以看到客户端分别向根服务器、.de服务器、lmu.de服务器发送查询
- 每次查询后客户端收到指向下一级服务器的响应，然后自己发起下一次查询
- 如果是递归查询，DNS服务器会代为完成所有查询

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 9: Frage 15 / 第15题

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 15  


#### 题目中文翻译 / 中文题意

在ISO/OSI参考模型的第四层，DNS使用哪种传输协议…
…用于区域传输？

#### 德文原题

```text
### Frage 15 / 第15题

**Welches Transportprotokoll wird auf Schicht 4 des ISO/OSI-Referenzmodells bei DNS …**  
**在ISO/OSI参考模型的第四层，DNS使用哪种传输协议…**

**(a) … für Zonentransfer eingesetzt?**  
**…用于区域传输？**
```

#### 解答

**Lösung / 答案：** **TCP**

区域传输（Zone Transfer）涉及大量数据传输，需要TCP的可靠传输保证。

**(b) … für DNS-Anfragen empfohlen?**  
**…推荐用于DNS查询？**

**Lösung / 答案：** **UDP**

普通DNS查询数据量小（通常小于512字节），UDP更高效，延迟更低。

---

**IV. Zusammenspiel verschiedener Protokolle / 不同协议的协作**

**场景描述：**

- 网络由两个以太网组成，通过组件X连接
- 客户端通过浏览器访问名为www的Web服务器
- 客户端知道DNS服务器的IP地址，只知道主机名www
- DNS服务器知道所有主机名和对应的IP地址

**网络参数：**

- DNS-Server: 10.10.8.2/24, MAC: 00:30:05:79:55:0A
- WWW-Server: 10.10.8.3/24, MAC: 00:30:05:79:55:DD
- Client: MAC: 00:30:05:79:55:55
- Port X1: MAC: 00:30:05:79:55:B1
- Port X2: MAC: 00:30:05:79:55:B2

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 10: Frage 16 / 第16题

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 16  


#### 题目中文翻译 / 中文题意

假设组件X是Switch，客户端IP为10.10.8.4（24位网络ID）
客户端发送包含DNS请求的帧到哪个MAC地址？

#### 德文原题

```text
### Frage 16 / 第16题

**假设组件X是Switch，客户端IP为10.10.8.4（24位网络ID）**

**(a) An welche MAC-Adresse sendet der Client Rahmen, die DNS-Anfragen enthalten?**  
**客户端发送包含DNS请求的帧到哪个MAC地址？**
```

#### 解答

**Lösung / 答案：** **00:30:05:79:55:0A**（DNS服务器的MAC地址）

**解释：** 当X是Switch时，客户端和DNS服务器在同一子网，客户端直接发送到DNS服务器的MAC地址。

**(b) 同上问题**

**Lösung / 答案：** **00:30:05:79:55:0A**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 11: Frage 17 / 第17题

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 17  


#### 题目中文翻译 / 中文题意

假设组件X是Router
客户端发送HTTP请求到哪个MAC地址？

#### 德文原题

```text
### Frage 17 / 第17题

**假设组件X是Router**

- Client: 192.168.1.2/28
- Port X1: 192.168.1.1/28
- Port X2: 10.10.8.1/24

**(a) An welche MAC-Adresse sendet der Client HTTP-Anfragen?**  
**客户端发送HTTP请求到哪个MAC地址？**
```

#### 解答

**Lösung / 答案：** **00:30:05:79:55:B1**（路由器端口X1的MAC地址）

**解释：** WWW服务器（10.10.8.3）与客户端（192.168.1.2）不在同一子网，需要通过路由器转发。客户端将帧发送到默认网关（路由器）的MAC地址。

**(b) An welche IPv4-Adresse sendet der Client DNS-Anfragen?**  
**客户端发送DNS请求到哪个IPv4地址？**

**Lösung / 答案：** **10.10.8.2**（DNS服务器的IP地址）

**解释：** IP地址不变，只是二层MAC地址通过路由器转换。

**(c) An welche MAC-Adresse versendet der Router einen Rahmen, mit der Ziel-IP 10.10.8.3 und dem Ziel-UDP-Port 53?**  
**路由器将目标IP为10.10.8.3、目标UDP端口为53的帧发送到哪个MAC地址？**

**Lösung / 答案：** **00:30:05:79:55:0A**（DNS服务器的MAC地址）

**解释：** 目标IP是10.10.8.3，但UDP端口53是DNS端口，所以这个包实际上可能是DNS响应或请求。根据目标IP，应该发送到10.10.8.3对应的MAC地址，即00:30:05:79:55:DD（WWW服务器）。

**注意：** 这道题可能有歧义。如果严格按IP地址来看，10.10.8.3对应WWW服务器（00:30:05:79:55:DD）。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 12: Frage 18 / 第18题

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 18  


#### 题目中文翻译 / 中文题意

为什么路由器的端口需要IPv4地址才能在网络间传输数据，而交换机的端口不需要？

#### 德文原题

```text
### Frage 18 / 第18题

**Warum benötigen die Ports eines Routers IPv4-Adressen um Daten zwischen den Netzen übertragen zu können, die Ports eines Switches jedoch nicht?**  
**为什么路由器的端口需要IPv4地址才能在网络间传输数据，而交换机的端口不需要？**
```

#### 解答

**Lösung / 答案：**

**路由器工作在第3层（网络层）**，需要参与IP路由决策。每个端口必须有IP地址，因为：

- 路由器需要作为不同子网的网关
- 主机将数据包发送到路由器的IP地址（下一跳）
- 路由器需要处理ICMP消息（如TTL超时）

**交换机工作在第2层（数据链路层）**，只根据MAC地址转发帧：

- 交换机对帧是"透明"的
- 不参与IP层的决策
- 只需要MAC地址表来转发帧

---

**V. Fragmentierung / 分片**

**网络拓扑：**

- E1 ←(2500B)→ R1 ←(280B)→ R2 ←(1500B)→ E2
- Kanal A: 2500 Bytes MTU
- Kanal B: 280 Bytes MTU
- Kanal C: 1500 Bytes MTU

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 13: Frage 5 / 第5题

**类型：** 考卷  

**来源说明：** Klausur 2019, Frage 5  


#### 题目中文翻译 / 中文题意

列举两个应用层协议。

#### 德文原题

```text
### Frage 5 / 第5题

**Nennen Sie zwei Protokolle der Anwendungsschicht. (2分)**  
**列举两个应用层协议。**
```

#### 解答

**Lösung / 答案：**（任选两个）

- **HTTP** (Hypertext Transfer Protocol) - 超文本传输协议
- **HTTPS** (HTTP Secure)
- **FTP** (File Transfer Protocol) - 文件传输协议
- **SMTP** (Simple Mail Transfer Protocol) - 简单邮件传输协议
- **DNS** (Domain Name System) - 域名系统
- **SSH** (Secure Shell)
- **Telnet**
- **DHCP** (Dynamic Host Configuration Protocol)
- **SNMP** (Simple Network Management Protocol)

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 14: Frage 7 / 第7题

**类型：** 考卷  

**来源说明：** Klausur 2019, Frage 7  


#### 题目中文翻译 / 中文题意

根据Internet模型，捕获中包含哪些应用层协议？

#### 德文原题

```text
### Frage 7 / 第7题

**(a) Welche Protokolle der Anwendungsschicht gemäß dem Internetmodell sind in dem Mitschnitt enthalten? (2分)**  
**根据Internet模型，捕获中包含哪些应用层协议？**
```

#### 解答

**Lösung / 答案：** **DNS, HTTP**

从表中可以看到DNS查询（#1-4）和HTTP请求/响应（#8, #10）。

---

**(b) Welche PDUs waren am TCP drei-Wege-Handschlag beteiligt? (1分)**  
**哪些PDU参与了TCP三次握手？**

**Lösung / 答案：** **#5, #6, #7**

- #5: SYN（客户端→服务器）
- #6: SYN ACK（服务器→客户端）
- #7: ACK（客户端→服务器）

---

**(c) Über welches Protokoll der Schicht 4 wurde PDU 1 übertragen? (1分)**  
**PDU 1通过哪个第4层协议传输？**

**Lösung / 答案：** **UDP**

DNS查询通常使用UDP端口53。

---

**(d) Wie lautet die IPv6-Adresse des Hosts [www.gnu.org](http://www.gnu.org/)? (1分)**  
**[www.gnu.org的IPv6地址是什么？](http://www.gnu.xn--orgipv6%3F-pk0mv2b333drfan52qh19b/)**

**Lösung / 答案：** **2001:470:142:3::5**

从#3和#4的DNS响应中可以看到AAAA记录。

---

**(e) Über welches Protokoll der Schicht 3 (inklusive Version) findet die HTTP-Anfrage und -Antwort statt? (1分)**  
**HTTP请求和响应通过哪个第3层协议（包括版本）进行？**

**Lösung / 答案：** **IPv4**

HTTP通信（#8-#11）使用的IP地址是209.51.188.148（IPv4格式）。

---

**(f) An welcher PDU kann man erkennen, dass die HTTP-Anfrage erfolgreich war? (1分)**  
**从哪个PDU可以看出HTTP请求成功？**

**Lösung / 答案：** **#10**

#10包含 "HTTP/1.1 200 OK"，状态码200表示请求成功。

---

**(g) Wie viele Bytes enthält die HTTP-PDU der Antwort (#10)? Hinweis: die PDU wurde per Ethernet übertragen. Eine Ethernet-PCI ist 14 Bytes lang. Im TCP-Header sind 12 Byte für Optionen genutzt. (1分)**  
**HTTP响应（#10）的PDU包含多少字节？提示：PDU通过以太网传输。以太网PCI为14字节。TCP头部中有12字节用于选项。**

**Lösung / 答案：**

**计算过程：**

- 总长度（Länge）= 407 字节
- 以太网头部 = 14 字节
- IPv4头部 = 20 字节
- TCP头部 = 20 + 12 = 32 字节（标准20字节 + 12字节选项）

HTTP数据长度 = 407 - 14 - 20 - 32 = **341 字节**

---

**III. Domain Name System (DNS) (6分)**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 15: Frage 8 / 第8题

**类型：** 考卷  

**来源说明：** Klausur 2019, Frage 8  


#### 题目中文翻译 / 中文题意

本题围绕“DNS、HTTP、应用层协议与 Wireshark”展开；请先阅读德文原题，再结合下方解答理解题意和考点。

#### 德文原题

```text
### Frage 8 / 第8题

**(a) Wie viele DNS-Anfragen waren nötig, um den Hostnamen [www.ifi.lmu.de](http://www.ifi.lmu.de/) wie oben abgebildet aufzulösen? (1分)**
```

#### 解答

**解析主机名 [www.ifi.lmu.de](http://www.ifi.lmu.de/) 需要多少次DNS查询？**

**Lösung / 答案：** **4**

查询过程：

1. 查询根服务器（.）→ 获取 .de 的NS记录
2. 查询 .de 服务器 → 获取 lmu.de 的NS记录
3. 查询 lmu.de 服务器 → 获取 ifi.lmu.de 的信息
4. 查询 ifi.lmu.de 服务器 → 获取 [www.ifi.lmu.de](http://www.ifi.lmu.de/) 的最终记录

---

**(b) Die URL [http://www.ifi.lmu.de/](http://www.ifi.lmu.de/) soll in einem Web-Browser angezeigt werden. Geben Sie die IPv4-Adresse des Rechners an, an den die HTTP-Anfrage gestellt wird. (1分)**  
**在浏览器中显示URL [http://www.ifi.lmu.de/](http://www.ifi.lmu.de/) 时，HTTP请求发送到哪个IPv4地址？**

**Lösung / 答案：** **141.84.94.144**

从DNS输出：

- [www.ifi.lmu.de](http://www.ifi.lmu.de/) → CNAME → stellenbosch.tcs.ifi.lmu.de
- stellenbosch.tcs.ifi.lmu.de → A → 141.84.94.144

---

**(c) Ist die oben aufgeführte DNS-Anfrage rekursiv, iterativ oder hybrid? Warum? (2分)**  
**上述DNS查询是递归的、迭代的还是混合的？为什么？**

**Lösung / 答案：** **Iterativ / 迭代**

**理由：**

- 使用了 `dig +trace` 命令
- 客户端自己逐级向各DNS服务器发送查询
- 每次收到响应后，客户端自己决定下一步查询哪个服务器
- 如果是递归查询，DNS服务器会代为完成整个解析过程

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 16: Frage 9 / 第9题

**类型：** 考卷  

**来源说明：** Klausur 2019, Frage 9  


#### 题目中文翻译 / 中文题意

在Internet模型的第4层，DNS使用哪种传输协议...
...用于区域传输？

#### 德文原题

```text
### Frage 9 / 第9题

**Welches Transportprotokoll wird auf Schicht 4 des Internetmodells bei DNS...**  
**在Internet模型的第4层，DNS使用哪种传输协议...**

**(a) ... für Zonentransfers eingesetzt? (1分)**  
**...用于区域传输？**
```

#### 解答

**Lösung / 答案：** **TCP**

区域传输涉及大量数据，需要TCP的可靠传输。

**(b) ... für DNS-Anfragen empfohlen? (1分)**  
**...推荐用于DNS查询？**

**Lösung / 答案：** **UDP**

普通DNS查询数据量小，UDP更高效。

---

**IV. Zusammenspiel verschiedener Protokolle / 不同协议的协作 (7分)**

**网络参数：**

- DNS-Server: 10.10.8.2/24, MAC: 00:30:05:79:55:0A
- WWW-Server: 10.10.8.3/24, MAC: 00:30:05:79:55:DD
- Client: MAC: 00:30:05:79:55:55
- Port X1: MAC: 00:30:05:79:55:B1
- Port X2: MAC: 00:30:05:79:55:B2

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 17: Frage 10 / 第10题

**类型：** 考卷  

**来源说明：** Klausur 2019, Frage 10  


#### 题目中文翻译 / 中文题意

假设组件X是Switch，客户端IP为10.10.8.4（24位网络ID）
客户端发送包含DNS请求的帧到哪个MAC地址？

#### 德文原题

```text
### Frage 10 / 第10题

**假设组件X是Switch，客户端IP为10.10.8.4（24位网络ID）**

**(a) An welche MAC-Adresse sendet der Client Rahmen, die DNS-Anfragen enthalten? (1分)**  
**客户端发送包含DNS请求的帧到哪个MAC地址？**
```

#### 解答

**Lösung / 答案：** **00:30:05:79:55:0A**

客户端和DNS服务器在同一子网（10.10.8.0/24），Switch透明转发，客户端直接发送到DNS服务器的MAC地址。

**(b) An welche IPv4-Adresse sendet der Client Pakete, die HTTP-Anfragen enthalten? (1分)**  
**客户端发送包含HTTP请求的数据包到哪个IPv4地址？**

**Lösung / 答案：** **10.10.8.3**

WWW服务器的IP地址。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 18: Frage 11 / 第11题

**类型：** 考卷  

**来源说明：** Klausur 2019, Frage 11  


#### 题目中文翻译 / 中文题意

假设组件X是Router
客户端发送HTTP请求到哪个MAC地址？

#### 德文原题

```text
### Frage 11 / 第11题

**假设组件X是Router**

- Client: 192.168.1.2/28
- Port X1: 192.168.1.1/28
- Port X2: 10.10.8.1/24

**(a) An welche MAC-Adresse sendet der Client HTTP-Anfragen? (1分)**  
**客户端发送HTTP请求到哪个MAC地址？**
```

#### 解答

**Lösung / 答案：** **00:30:05:79:55:B1**

WWW服务器（10.10.8.3）与客户端（192.168.1.2）不在同一子网，需要通过默认网关（路由器X1端口）转发。

**(b) An welche IPv4-Adresse sendet der Client DNS-Anfragen? (1分)**  
**客户端发送DNS请求到哪个IPv4地址？**

**Lösung / 答案：** **10.10.8.2**

DNS服务器的IP地址不变，IP地址是端到端的。

**(c) An welche MAC-Adresse versendet der Router einen Rahmen, mit der Ziel-IP 10.10.8.3 und dem Ziel-UDP-Port 53? (1分)**  
**路由器将目标IP为10.10.8.3、目标UDP端口为53的帧发送到哪个MAC地址？**

**Lösung / 答案：** **00:30:05:79:55:DD**

目标IP是10.10.8.3，即WWW服务器，所以发送到其MAC地址。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 19: Frage 12 / 第12题

**类型：** 考卷  

**来源说明：** Klausur 2019, Frage 12  


#### 题目中文翻译 / 中文题意

为什么路由器的端口需要IPv4地址才能在网络间传输数据，而交换机的端口不需要？

#### 德文原题

```text
### Frage 12 / 第12题

**Warum benötigen die Ports eines Routers IPv4-Adressen um Daten zwischen den Netzen übertragen zu können, die Ports eines Switches jedoch nicht? (2分)**  
**为什么路由器的端口需要IPv4地址才能在网络间传输数据，而交换机的端口不需要？**
```

#### 解答

**Lösung / 答案：**

**路由器（第3层设备）：**

- 工作在网络层，需要处理IP数据包
- 每个端口连接不同的IP子网
- 需要IP地址作为该子网的网关
- 主机将数据包发送到路由器的IP地址（下一跳）
- 需要处理ICMP消息（如TTL超时等）

**交换机（第2层设备）：**

- 工作在数据链路层，只处理以太网帧
- 根据MAC地址表转发帧
- 对IP层是"透明"的
- 所有端口属于同一个广播域
- 不参与IP路由决策

---

**V. Fragmentierung / 分片 (7分)**

**网络拓扑：**

- E1 ←(2500B)→ R1 ←(280B)→ R2 ←(1500B)→ E2

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---
