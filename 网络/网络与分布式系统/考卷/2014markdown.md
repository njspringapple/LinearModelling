## 1. ISO OSI-Schichtenmodell / ISO OSI 分层模型

### Aufgabe / 题目

**Deutsch:**  
Ergänzen Sie die Namen der Schichten im ISO OSI-Schichtenmodell in Deutsch und Englisch und geben Sie je Schicht eine charakteristische Aufgabe an!

**中文：**  
请补全 ISO OSI 分层模型中各层的德文和英文名称，并为每一层给出一个典型任务。

---

### Lösung / 解答

**Hinweis / 提示：**  
Nachfolgend sind z. T. mehr als eine charakteristische Aufgabe je Schicht aufgeführt. Entsprechend der Aufgabenstellung wäre lediglich eine charakteristische Aufgabe notwendig, um diese Aufgabe vollumfänglich zu erfüllen.

下面有些层列出了不止一个典型任务。根据题目要求，每层只需要写出一个典型任务即可完整完成本题。

| Schicht / 层 | Name / 名称                                            | Charakteristische Aufgabe / 典型任务                                                                                                                                                  |
| ----------: | ---------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
|           7 | Anwendungsschicht, Application Layer / 应用层           | Allgemein verwendbare Dienste werden standardisiert und als Dienste und Protokolle spezifiziert / 标准化通用服务，并将其定义为服务和协议                                                             |
|           6 | Darstellungsschicht, Presentation Layer / 表示层        | Datenmodellierung in Objekten, Aushandeln der konkreten Transfersyntax, Abbilden lokaler konkreter Syntax, z. B. Basic Encoding Rules, BER / 对象中的数据建模，协商具体传输语法，将本地具体语法进行映射，例如 BER |
|           5 | Kommunikationssteuerungsschicht, Session Layer / 会话层 | Benutzeridentifizierung, Dialogführung und Synchronisation innerhalb einer Sitzung / 用户识别、会话中的对话控制与同步                                                                             |
|           4 | Transportschicht, Transport Layer / 传输层              | Verbindung zwischen zwei Prozessen, netzunabhängiger Transport von Nachrichten zwischen zwei Endsystemen / 两个进程之间的连接，在两个端系统之间进行与网络无关的消息传输                                         |
|           3 | Vermittlungsschicht, Network Layer / 网络层             | Wegewahl und Vermittlung / 路径选择与转发                                                                                                                                                |
|           2 | Sicherungsschicht, Data Link Layer / 数据链路层           | Zusammenfassung von Bits zu Blöcken/Frames, Fehlererkennung, ggf. Fehlerkorrektur, Medium Access Control / 将比特组合成块或帧，错误检测，必要时错误纠正，介质访问控制                                          |
|           1 | Bitübertragungsschicht, Physical Layer / 物理层         | Darstellung von Daten auf Medium, transparente Übertragung von Bits / 在介质上表示数据，透明传输比特                                                                                             |

---

## 2. Übertragungsraten / 传输速率

### Aufgabe / 题目

**Deutsch:**  
Über ein Medium mit einer Bandbreite von 1 MHz wird mit einer 2-Stufencodierung übertragen. Wie viele Bits pro Sekunde können maximal übertragen werden, wenn...

**中文：**  
通过一个带宽为 1 MHz 的介质，使用 2 级编码进行传输。问在以下情况下，每秒最多可以传输多少比特？

**Hinweis / 提示：**  
Geben Sie für jede Teilaufgabe jeweils eine Rechnung bzw. Formel und ein Ergebnis an!

请对每个小题分别给出计算过程或公式以及结果。

---

### (a) Kein Rauschen / 无噪声，理想介质

**Deutsch:**  
... kein Rauschen vorkommt, also ideales Medium?

**中文：**  
如果没有噪声，即为理想介质？

**Lösung / 解答：**

Ansatz mit Nyquist / 使用 Nyquist 公式：

$$
C = 2 \cdot B \cdot \log_2 M
$$

$$
C = 2 \cdot 1\,MHz \cdot \log_2 2
= 2 \cdot 10^6\,bit/s
$$

**Ergebnis / 结果：**

$$
C = 2\,Mbit/s
$$

---

### (b) Signal-Rausch-Verhältnis / 信噪比

**Deutsch:**  
... ein Verhältnis zwischen Signal und Rauschen von \(S/N = 1023\) vorherrscht?

**中文：**  
如果信号与噪声之比为 \(S/N = 1023\)？

**Lösung / 解答：**

Ansatz mit Shannon / 使用 Shannon 公式： - 香浓定理

$$
C = B \cdot \log_2(1 + S/N)
$$

$$
C = 1\,MHz \cdot \log_2(1 + 1023)
$$

$$
C = 1\,MHz \cdot \log_2(1024)
= 10^7\,bit/s
$$

**Ergebnis / 结果：**

$$
C = 10\,Mbit/s
$$

---

## 3. Codierungsverfahren / 编码方法

---

### (a) Manchester-Codierung / 曼彻斯特编码

**Deutsch:**  
Geben Sie das in Manchestercodierung dargestellte Bitmuster an!

**中文：**  
请给出图中曼彻斯特编码所表示的比特模式。

**Lösung / 解答：**

Bitmuster / 比特序列：

```text
0 1 1 0 1 0 0 1 1 0 0 1
````

或写作：

```text
011010011001
```

---

### (b) Übertragungsrate und Baud-Rate / 传输速率与波特率

**Deutsch:**  
Angenommen, das obige Bitmuster für Manchestercodierung wird in 1 ms übertragen.

**中文：**  
假设上述曼彻斯特编码的比特序列在 1 ms 内完成传输。

---

#### i. Übertragungsrate / 传输速率

**Deutsch:**  
Wie hoch ist die Übertragungsrate?

**中文：**  
传输速率是多少？

**Lösung / 解答：**

$\text{Übertragungsrate} = \frac{12\,Bit}{1\,ms} = 12000\,bit/s$

**Ergebnis / 结果：**

$12000\,bit/s$

---

#### ii. Baud-Rate / 波特率

**Deutsch:**  
Wie hoch ist die Baud-Rate des Signals?

**中文：**  
该信号的波特率是多少？

**Lösung / 解答：**

Bei Manchester-Codierung gilt / 对曼彻斯特编码：

$\text{Baud-Rate} = 2 \cdot \text{Übertragungsrate}$

$\text{Baud-Rate} = 2 \cdot 12000 = 24000\,Baud$

**Ergebnis / 结果：**

$24000\,Baud$
---

## 4. Ethernet, CSMA / 以太网与 CSMA

### Aufgabe / 题目

**Deutsch:**  
Gegeben sei ein Ethernet mit einer Übertragungsrate von 1 GBit/s, einer Leitungslänge von 1000 m und einer Signalgeschwindigkeit von 2⋅108 m/s2 \cdot 10^8\,m/s2⋅108m/s. Berechnen Sie die minimale Rahmengröße, bei der CSMA/CD als Vielfachzugriffsverfahren noch einsetzbar wäre. Geben Sie das Ergebnis in Bytes/Oktetten sowie den Rechenweg an!

**中文：**  
给定一个以太网，其传输速率为 1 GBit/s，线路长度为 1000 m，信号传播速度为 $2 \cdot 10^8\,m/s$。请计算在 CSMA/CD 多路访问机制仍可使用的情况下所需的最小帧大小。请以 Byte/Octet 给出结果，并写出计算过程。

**Hinweis / 提示：**

$1\,GBit = 10^9\,Bits$

---

### Lösung / 解答

#### Signallaufzeit hin und zurück / 信号往返传播时间

$T_{rtt-signal} = \frac{1\,km}{200000\,km/s} \cdot 2 = 10^{-5}\,s$

#### Minimale Framegröße / 最小帧大小

$S_{frame-size} = 10^{-5}\,s \cdot 10^9\,bit/s = 10^4\,bit$

$10^4\,bit = 10000\,bit$

$\frac{10000}{8} = 1250\,byte$

**Ergebnis / 结果：**

$1250\,Byte$

---

## 5. Internet Protocol / 互联网协议

---

### (a) Klassenbasierte IPv4-Adressierung / 基于类别的 IPv4 地址分配

#### i. Bestandteile einer IPv4-Adresse / IPv4 地址的组成部分

**Deutsch:**  
Ursprünglich wurde der Adressraum für Internetadressen in Klassen aufgeteilt. Aus welchen zwei Teilen besteht demzufolge eine IPv4-Adresse?

**中文：**  
最初，互联网地址空间被划分为不同类别。根据这种分类方式，一个 IPv4 地址由哪两个部分组成？

**Lösung / 解答：**

```text
Netz-ID, Host-ID
```

中文：

```text
网络 ID，主机 ID
```

---

#### ii. Vorteil und Nachteil / 优点与缺点

**Deutsch:**  
Nennen Sie einen Vorteil und einen Nachteil der klassenbasierten Adressvergabe.

**中文：**  
请说出基于类别的地址分配方式的一个优点和一个缺点。

**Lösung / 解答：**

**Vorteil / 优点：**  
Anhand Netz-ID/Präfix können schnell Routing-Entscheidungen getroffen werden; Adressraum leichter zu verwalten.

根据网络 ID 或前缀可以快速做出路由决策；地址空间更容易管理。

**Nachteil / 缺点：**  
Großer Teil der Adressen bleibt unbenutzt.

大量地址会被浪费或未被使用。

---

#### iii. CIDR / 无类别域间路由

**Deutsch:**  
Mit CIDR wurde ein flexibleres Schema für die Vergabe von Adressräumen benutzt. Worin besteht der Unterschied zur klassenbasierten Aufteilung des Adressraums?

**中文：**  
CIDR 使用了一种更灵活的地址空间分配方案。它与基于类别的地址空间划分有什么区别？

**Lösung / 解答：**

Die Länge der Netz-ID ist variabel und nicht an Klassen gebunden.

网络 ID 的长度是可变的，不再受固定类别限制。

---

#### iv. Maximale Netz-ID-Länge / 最大网络 ID 长度

**Deutsch:**  
Wie lang in Bits darf eine Netz-ID für ein IPv4-basiertes Subnetz mit 58 Hosts höchstens sein?

**中文：**  
对于一个需要容纳 58 台主机的 IPv4 子网，网络 ID 最多可以是多少位？

**Lösung / 解答：**

```text
26
```

**Erklärung / 说明：**

58 Hosts benötigen mindestens 6 Host-Bits, denn:

$2^6 - 2 = 62$

IPv4 hat 32 Bits:

$32 - 6 = 26$

**Ergebnis / 结果：**

$/26$

---

#### v. Netzmaske / 子网掩码

**Deutsch:**  
Wie lautet die Netzmaske für das Subnetz 192.168.218.0/28? Machen Sie ihre Angabe in der Form r.s.p.q mit r,s,p,q∈{0,…,255}r,s,p,q \in \{0,\dots,255\}r,s,p,q∈{0,…,255}, d. h. in dezimaler Schreibweise.

**中文：**  
子网 192.168.218.0/28 的子网掩码是多少？请以十进制点分形式 r.s.p.q给出，其中 $r,s,p,q \in \{0,\dots,255\}$。

**Lösung / 解答：**

```text
255.255.255.240
```

---

### (b) Fragmentierung / 分片

**Deutsch:**  
Nennen Sie einen Fall, in dem IPv4-Pakete fragmentiert werden müssen!

**中文：**  
请举出一种 IPv4 数据包必须被分片的情况。

**Lösung / 解答：**

Wenn Paketlänge > MTU auf dem Pfad.

当数据包长度大于路径上的 MTU 时。

---

### (c) ICMP-Meldungen / ICMP 消息

**Deutsch:**  
Im Internet kann mittels des Internet Control Message Protocol, ICMP, signalisiert werden, dass kein Weg zum Ziel eines IP-Paketes ermittelt werden kann, also destination unreachable. Nennen Sie zwei weitere Meldungen, die mittels ICMP gesendet bzw. empfangen werden können!

**中文：**  
在互联网中，可以通过 ICMP 协议通知无法找到到达某个 IP 数据包目标的路径，即 destination unreachable。请再列举两个可以通过 ICMP 发送或接收的消息。

**Lösung / 解答：**

Zum Beispiel / 例如：

- echo request / 回显请求
- time exceeded / 超时

---

### (d) Routing-Protokolle / 路由协议

**Deutsch:**  
Zwischen autonomen Systemen werden andere Routing-Protokolle eingesetzt als innerhalb. Nennen Sie einen Grund dafür mit kurzer Erklärung!

**中文：**  
自治系统之间使用的路由协议不同于自治系统内部使用的路由协议。请给出一个原因并简要说明。

**Lösung / 解答：**

Zum Beispiel / 例如：

- verschiedene Metriken/Entscheidungskriterien für Wegewahl
- 路径选择时采用不同的度量标准或决策准则

---

## 6. Transmission Control Protocol / 传输控制协议 TCP

### Aufgabe / 题目

**Deutsch:**  
Das Diagramm zeigt die Zustände und Zustandsübergänge in einem TCP-basierten Client. Ergänzen Sie den Text in den weißen Flächen!

**中文：**  
该图展示了一个基于 TCP 的客户端中的状态和状态转换。请补全白色框中的文字。

---

### Lösung / 解答

| Zustand / 状态           | Übergang / 转换                | Bedeutung / 含义 |
| ---------------------- | ---------------------------- | -------------- |
| CLOSED                 | SYN senden                   | 发送 SYN，发起连接    |
| SYN_SENT               | SYNACK empfangen             | 接收到 SYN-ACK    |
| SYN_SENT → ESTABLISHED | ACK senden                   | 发送 ACK，连接建立    |
| ESTABLISHED            | FIN senden                   | 发送 FIN，开始关闭连接  |
| FIN_WAIT_1             | ACK empfangen, nichts senden | 收到 ACK，不发送数据   |
| FIN_WAIT_2             | FIN empfangen                | 收到对方 FIN       |
| TIME_WAIT              | ACK senden                   | 发送 ACK         |
| TIME_WAIT → CLOSED     | warten                       | 等待后进入 CLOSED   |

---

### Ergänzte Einträge / 补全项

```text
SYN senden
SYNACK empfangen
ACK senden
ESTABLISHED
FIN senden
ACK empfangen (nichts senden)
FIN empfangen
ACK senden
warten
```

---

## 7. E-Mail / 电子邮件

---

### (a) E-Mail-Protokolle / 电子邮件协议

**Deutsch:**  
Beschriften Sie alle Pfeile in der Zeichnung mit den entsprechenden E-Mail-Protokollen.

**中文：**  
请用相应的电子邮件协议标注图中的所有箭头。

---

### Lösung / 解答

| Verbindung / 连接                                      | Protokoll / 协议 |
| ---------------------------------------------------- | -------------- |
| Sender User Agent → Mail Server des Senders          | SMTP           |
| Mail Server des Senders → Mail Server des Empfängers | SMTP           |
| Mail Server des Empfängers → Empfänger User Agent    | IMAP oder POP3 |

中文说明：

|连接|协议|
|---|---|
|发送方用户代理 → 发送方邮件服务器|SMTP|
|发送方邮件服务器 → 接收方邮件服务器|SMTP|
|接收方邮件服务器 → 接收方用户代理|IMAP 或 POP3|

---

### (b) Dienstgüteparameter / 服务质量参数

**Deutsch:**  
Internet E-Mail ist empfindlich gegen den Dienstgüteparameter „Datenverlust“ des Transportnetzes. Nennen Sie zwei Dienstgüteparameter, gegen die E-Mail unempfindlich ist und begründen Sie.

**中文：**  
互联网电子邮件对传输网络中的“数据丢失”这一服务质量参数较为敏感。请列举两个电子邮件不敏感的服务质量参数，并说明原因。

---

### Lösung / 解答

Zum Beispiel / 例如：

- **Latenz**, da kein Echtzeitdienst  
    **延迟**，因为电子邮件不是实时服务
    
- **Jitter**, da kein Echtzeitdienst  
    **抖动**，因为电子邮件不是实时服务
    

---

## 8. Kommunikationsablauf / 通信过程

### Aufgabe / 题目

**Deutsch:**  
Das in der Abbildung skizzierte Netz besteht aus zwei Ethernets, die so mit einem Router verbunden sind, dass IPv4-Pakete zwischen ihnen vermittelt werden. Auf dem Client wird ein Browser-Programm ausgeführt, das eine Verbindung zu einem Webserver namens `www` aufbaut, um ein HTML-Dokument abzurufen.

**中文：**  
图中所示网络由两个以太网组成，它们通过一个路由器连接，使得 IPv4 数据包可以在两个网络之间转发。客户端上运行一个浏览器程序，该程序与名为 `www` 的 Web 服务器建立连接，以获取一个 HTML 文档。

---

### Hinweise / 提示

**Deutsch:**  
Benutzen Sie beim Eintragen in die Tabelle von...

- MAC-Adressen nur das letzte Byte,
- IP-Adressen nur die letzten 2 Byte,
- Broadcast als die Angabe B-Cast.
- Die erste Tabellenzeile ist als Beispiel vorgegeben.

Der Client kennt:

- die IP-Adresse seines lokalen DNS-Servers,
- die URL des abzufragenden Web-Objekts und
- eine Default-Route über .1.1.

Der DNS-Server ist autoritativ für alle Teilnehmer in der Abbildung.

Caches, also ARP, DNS usw.:

- es existieren keine aktuellen Cache-Werte,
- empfangene aufgelöste Adressen werden aggressiv zwischengespeichert und müssen nicht wieder angefragt werden.

Eine PDU einer Schicht NNN passt immer in eine PDU der Schicht N−1N-1N−1.

Vernachlässigen Sie Übertragungsfehler, Verluste oder verworfene Nachrichten.

**中文：**  
填写表格时请使用：

- MAC 地址只写最后 1 个 Byte，
- IP 地址只写最后 2 个 Byte，
- 广播写作 B-Cast。
- 第一行作为示例已给出。

客户端已知：

- 本地 DNS 服务器的 IP 地址，
- 要请求的 Web 对象的 URL，
- 通过 .1.1 的默认路由。

DNS 服务器对图中所有参与者具有权威性。

缓存，例如 ARP、DNS 等：

- 当前不存在缓存值，
- 接收到的已解析地址会被积极缓存，不需要再次查询。

第 N 层的 PDU 总能放入第 N−1层的 PDU 中。

忽略传输错误、丢失或被丢弃的消息。

---

### Netzwerkteilnehmer / 网络参与者

| Gerät / 设备                     | IP-Adresse / IP 地址 | MAC-Adresse / MAC 地址 |
| ------------------------------ | ------------------ | -------------------- |
| Client                         | 192.168.1.2        | 00:30:05:79:55:C0    |
| Router, Client-Seite / 路由器客户端侧 | 192.168.1.1        | 00:30:05:79:55:A1    |
| Router, Server-Seite / 路由器服务器侧 | 192.168.2.1        | 00:30:05:79:55:A2    |
| DNS Server                     | 192.168.2.2        | 00:30:05:79:55:D0    |
| Web-Server `www`               | 192.168.2.3        | 00:30:05:79:55:E0    |

---

### Lösung / 解答表格

|Pkt|MAC-Adr von / MAC 源|MAC-Adr zu / MAC 目的|IP-Adr von / IP 源|IP-Adr zu / IP 目的|Port von / 源端口|Port zu / 目的端口|TCP Flags / TCP 标志|Payload / Erklärung / 负载或说明|
|--:|---|---|---|---|---|---|---|---|
|1|:C0|B-Cast|-|-|-|-|-|ARP: wer hat .1.1? / ARP：谁有 .1.1？|
|2|:A1|:C0|-|-|-|-|-|ARP: ich habe .1.1! / ARP：我有 .1.1！|
|3|:C0|:A1|.1.2|.2.2|X|dns(53)|-|DNS Query: www? / DNS 查询：www？|
|4|:A2|B-Cast|-|-|-|-|-|ARP: wer hat .2.2? / ARP：谁有 .2.2？|
|5|:D0|:A2|-|-|-|-|-|ARP: ich habe .2.2! / ARP：我有 .2.2！|
|6|:A2|:D0|.1.2|.2.2|12345|dns(53)|-|DNS Query: www? / DNS 查询：www？|
|7|:D0|:A2|.2.2|.1.2|dns(53)|12345|-|DNS Response: .2.3 / DNS 响应：.2.3|
|8|:A1|:C0|.2.2|.1.2|dns(53)|12345|-|DNS Response: .2.3 / DNS 响应：.2.3|
|9|:C0|:A1|.1.2|.2.3|4711|www(80)|SYN|Conn-Req / 连接请求|
|10|:A2|B-Cast|-|-|-|-|-|ARP: wer hat .2.3? / ARP：谁有 .2.3？|
|11|:E0|:A2|-|-|-|-|-|ARP: ich habe .2.3! / ARP：我有 .2.3！|
|12|:A2|:E0|.1.2|.2.3|4711|www(80)|SYN|Conn-Req / 连接请求|
|13|:E0|:A2|.2.3|.1.2|www(80)|4711|SYN,ACK|Conn-Req-Ack / 连接请求确认|
|14|:A1|:C0|.2.3|.1.2|www(80)|4711|SYN,ACK|Conn-Req-Ack / 连接请求确认|
|15|:C0|:A1|.1.2|.2.3|4711|www(80)|ACK|Conn-Est / 连接建立|
|16|:A2|:E0|.1.2|.2.3|4711|www(80)|ACK|Conn-Est / 连接建立|

---

## Begriffe / 术语对照

|Deutsch|中文|
|---|---|
|Schicht|层|
|Anwendungsschicht|应用层|
|Darstellungsschicht|表示层|
|Kommunikationssteuerungsschicht|会话层|
|Transportschicht|传输层|
|Vermittlungsschicht|网络层|
|Sicherungsschicht|数据链路层|
|Bitübertragungsschicht|物理层|
|Übertragungsrate|传输速率|
|Bandbreite|带宽|
|Rauschen|噪声|
|Signal-Rausch-Verhältnis|信噪比|
|Codierungsverfahren|编码方法|
|Manchester-Codierung|曼彻斯特编码|
|Baud-Rate|波特率|
|Rahmengröße|帧大小|
|Netzmaske|子网掩码|
|Fragmentierung|分片|
|Routing-Protokoll|路由协议|
|User Agent|用户代理|
|Mail Server|邮件服务器|
|Latenz|延迟|
|Jitter|抖动|
|Datenverlust|数据丢失|
|Broadcast|广播|
|ARP|地址解析协议|
|DNS|域名系统|
|TCP|传输控制协议|
|SYN|同步请求|
|SYN,ACK|同步确认|
|ACK|确认|