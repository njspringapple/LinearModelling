# Uebungsblatt 3 - 中德对照解答

## 1. Protokollschichtung / 协议分层

![Blatt 03 Seite 1: Internetmodell und OSI-Modell](pictures/blatt-03_pages-1-2-1.png)

**(a) OSI-Modell / OSI 七层**

| Nr. | Deutsch | 中文 |
|---:|---|---|
| 7 | Anwendung | 应用层 |
| 6 | Darstellung | 表示层 |
| 5 | Sitzung | 会话层 |
| 4 | Transport | 传输层 |
| 3 | Vermittlung | 网络/网际层 |
| 2 | Sicherung | 数据链路层 |
| 1 | Bituebertragung | 物理层 |

**(b) Aufgaben / 任务**

**DE:** Die Anwendungsschicht stellt konkrete Anwendungsprotokolle bereit. Die Darstellungsschicht behandelt Kodierung, Datenformat und Verschluesselung. Die Sitzungsschicht verwaltet Dialoge und Sitzungen. Die Transportschicht bietet Ende-zu-Ende-Kommunikation zwischen Prozessen. Die Vermittlungsschicht uebernimmt Adressierung und Routing. Die Sicherungsschicht uebertraegt Frames auf einem lokalen Link. Die Bituebertragungsschicht uebertraegt Signale.

**中文：** 应用层提供具体应用协议；表示层处理编码、格式、加密；会话层管理会话状态；传输层提供端到端进程通信；网络层负责寻址和路由；链路层负责本地链路上的帧传输、MAC、差错检测；物理层负责信号、电气/光学/无线传输。可靠性、分段、寻址等功能可能在多层出现，但作用范围不同。

**(c) Vorteile/Nachteile / 优缺点**

**DE:** Vorteile sind Modularitaet, Austauschbarkeit und Standardisierung; jede Schicht muss nur den Dienst der darunterliegenden Schicht kennen. Nachteile sind zusaetzliche Header, Verarbeitungsaufwand und manchmal doppelte Funktionen.

**中文：** 优点是模块化、可替换、易标准化；每层只依赖下层服务。缺点是会有额外头部和处理开销；严格分层有时导致信息重复或性能优化困难。

**(d) OSI vs Internet-Anwendungsschicht**

**DE:** Das Internetmodell fasst OSI-Schichten 5-7 meist zur Anwendungsschicht zusammen. Anwendungen muessen daher Funktionen wie Darstellung, Verschluesselung oder Sitzungsverwaltung oft selbst oder ueber Bibliotheken leisten.

**中文：** Internet 模型通常把 OSI 的第 5 到第 7 层合并为应用层。因此表示、加密、会话管理等功能往往由应用程序自身或应用库来完成。

## 2. Datenpakete in Python / Python 数据包

Gegeben:

```python
msg = bytes([0x48, 0x65, 0x6c, 0x6c, 0x6f])
recipient = ("192.168.1.135", 6243)
sendto(msg, recipient)
```

**(a)**  
**DE:** `msg` ist Nutzdaten/SDU der Schicht N; `recipient` ist Schnittstelleninformation/ICI, denn sie steuert, wohin gesendet wird.

**中文：** `msg` 是第 N 层要传送的有效载荷，也就是 `(N)-SDU`。`recipient` 包含目的 IP 和端口，用来告诉服务应该发往哪里，因此属于接口控制信息 `(N)-ICI`。

**(b)**  
**DE:** Die `(N)-PDU` besteht aus `(N)-PCI + (N)-SDU`. Aus dem Funktionsaufruf allein sieht man die fertige PDU nicht, weil Header/PCI intern von `sendto` erzeugt werden.

**中文：** `(N)-PDU` 由本层控制信息 `(N)-PCI` 加上本层服务数据单元 `(N)-SDU` 构成。函数调用中只能看到 SDU 和目的信息，看不到完整 PDU，因为首部通常在 `sendto` 内部生成。

**(c)**  
**DE:** Die `(N+1)-PDU` in `msg` wird auf Schicht N als `(N)-SDU` behandelt und durch Encapsulation mit `(N)-PCI` zu einer `(N)-PDU`.

**中文：** 上一层交下来的 `(N+1)-PDU` 到了第 N 层后，被第 N 层当作自己的 `(N)-SDU`。第 N 层再加上自己的 PCI 进行封装，得到 `(N)-PDU`。

**(d)** ASCII:

```text
0x48 0x65 0x6c 0x6c 0x6f = Hello
```

## 3. Bestandteile des Schichtenmodells / 分层模型的数据单元

![Blatt 03 Seite 2: Schnittbildung und Modellgrafiken](pictures/blatt-03_pages-1-2-2.png)

**(a)**  
**DE:** Eine PDU entsteht, indem eine Schicht zu einer SDU eigene Protokollkontrollinformationen (PCI, Header/Trailer) hinzufuegt.

**中文：** PDU 的形成过程就是封装：某一层把上层交来的 SDU 加上自己的控制信息 PCI，例如首部或尾部，形成该层的协议数据单元 PDU。

**(b)**

1. **DE:** Die Peer-Entity befindet sich auf derselben Schicht N im anderen System.  
   **中文：** 对等实体位于另一台系统的同一层 N。
2. **DE:** Die `(N)-PDU` wird fuer Schicht `N-1` zur `(N-1)-SDU`.  
   **中文：** 第 N 层的 PDU 交给下一层后，会被第 `N-1` 层当作 SDU。
3. **DE:** Beim Senden wandern Anwendungsdaten nach unten, jede Schicht kapselt ein. Beim Empfangen wandern Daten nach oben, jede Schicht entfernt und interpretiert ihren Header.  
   **中文：** 发送时数据自上而下，每层添加自己的头部；接收时数据自下而上，每层解释并去掉自己的头部。

**(c) Schnittarten**

**DE:** Die Abbildung mit LMU/DFN/MIT beschreibt Systemschnitte, weil reale Systeme oder Netze voneinander abgegrenzt werden. Die Abbildung HTTP/TCP/IP/Ethernet zeigt Dienstschnitte zwischen Schichten und Protokollschnitte horizontal zwischen Peer-Protokollen.

**中文：** LMU/DFN/MIT 那张图强调不同真实系统或网络之间的边界，因此对应系统切分。HTTP/TCP/IP/Ethernet 那张图中，垂直方向的层间边界是服务接口，水平方向同层实体之间的逻辑通信是协议接口。

## 4. ISO/OSI-Kritik / 对 OSI 的批评

**(a)**  
**DE:** Besonders umstritten sind Sitzung und Darstellung, weil viele reale Internetanwendungen diese Aufgaben selbst oder gar nicht getrennt implementieren.

**中文：** 最有争议的是会话层和表示层，因为许多真实的互联网应用并不会把这些功能独立成层，而是由应用自己、库或具体协议一起完成。

**(b)**  
**DE:** Kritikpunkte: zu komplex, spaet standardisiert, wenig implementierungsnah, teilweise schwerfaellige Protokollfamilie; TCP/IP war pragmatischer und frueher breit eingesetzt.

**中文：** 批评点包括：模型和协议族过于复杂，标准化较晚，离实际实现较远，工程上不够轻量。TCP/IP 更务实，也更早被广泛部署。

**(c)**  
**DE:** Der Wert des OSI-Modells liegt in seiner Rolle als didaktisches Referenzmodell mit klaren Begriffen fuer Kapselung, Schnittstellen und Verantwortlichkeiten.

**中文：** OSI 的价值主要在教学和概念整理：它清楚地区分封装、接口、层次职责，帮助理解网络系统如何分工。

**(d)**  
**DE:** Fuer reale Protokolle sollte das Internetmodell im Vordergrund stehen; das OSI-Modell bleibt ein begriffliches Werkzeug. Am besten werden beide gemeinsam gelehrt.

**中文：** 实际协议学习应以 Internet 模型为主，因为真实互联网主要按 TCP/IP 工作；OSI 模型适合作为概念工具。最好的教学方式是两者结合。

**Wissen / 知识点：** OSI 是解释分层的“地图”，TCP/IP 是互联网工程中真正主导的协议族。
