
# Rechnernetze und verteilte Systeme  
# 计算机网络与分布式系统

## Gedächtnisprotokoll Klausur SoSe 2013  
## 2013 年夏季学期考试回忆版

> **Deutsch:** Die Reihenfolge der Fragen ist vermutlich nicht ganz richtig und es fehlen eventuell einige Fragen.  
> **中文：** 题目顺序可能并不完全正确，并且可能缺少一些题目。  
>
> **Deutsch:** Rote Stellen bedeuten, dass die Angabe unsicher war.  
> **中文：** 红色部分表示原始资料中不确定的内容。

---

# 1. Offene Fragen  
# 1. 开放题

## 1.1 ISO/OSI-Schichten eines Routers  
## 1.1 路由器实现的 ISO/OSI 层

**Frage / 问题**

> Nennen Sie die ISO/OSI-Schicht(en), die ein Router implementiert!  
> 请说明路由器实现了哪些 ISO/OSI 层。

**Antwort / 解答**

Ein Router implementiert mindestens:

| OSI-Schicht | Deutsch                | 中文    |
| ----------: | ---------------------- | ----- |
|           1 | Bitübertragungsschicht | 物理层   |
|           2 | Sicherungsschicht      | 数据链路层 |
|           3 | Vermittlungsschicht    | 网络层   |

**Kurzantwort / 简答**

```text
Schicht 1 bis 3, hauptsächlich Schicht 3.
第 1 到第 3 层，核心是第 3 层网络层。
````

---

## 1.2 MTU-Schicht

## 1.2 MTU 所在层

**Frage / 问题**

> Auf welcher Schicht wird die MTU festgelegt?  
> MTU 是在哪一层确定的？

**Antwort / 解答**

```text
Schicht 2, Sicherungsschicht.
第 2 层，数据链路层。
```

**Erklärung / 说明**

Die MTU beschreibt die maximale Größe einer Nutzlast, die über ein bestimmtes Link-Layer-Netz übertragen werden kann.  
MTU 描述的是某条链路能够承载的最大网络层数据包大小，因此通常由数据链路层技术决定。

---

## 1.3 Koppelkomponente, die fehlerhafte Rahmen normalerweise nicht verwirft

## 1.3 通常不丢弃错误帧的连接设备

**Frage / 问题**

> Nennen Sie eine Koppelkomponente in verteilten Systemen, die normalerweise keine Pakete mit fehlerhaften Rahmen verwirft!  
> 请举出一种通常不会因为帧错误而丢弃包的连接组件。

**Antwort / 解答**

```text
Repeater oder Hub.
中继器或集线器。
```

**Erklärung / 说明**

Ein Repeater bzw. Hub arbeitet auf Schicht 1 und interpretiert keine vollständigen Rahmen. Daher prüft er normalerweise keine Prüfsummen wie FCS.  
中继器或集线器工作在第 1 层，不理解完整的数据链路层帧，因此通常不会检查帧校验序列 FCS，也就不会基于帧错误主动丢弃。

---

## 1.4 Nyquist

## 1.4 奈奎斯特定理

**Frage / 问题**

> Was ist laut dem Satz von Nyquist die beschränkende Größe für die Übertragung in einem rauschfreien Medium?  
> 根据奈奎斯特定理，在无噪声信道中限制传输的量是什么？

**Antwort / 解答**

```text
Die Bandbreite des Mediums.
介质的带宽。
```

**Erklärung / 说明**

Für ein rauschfreies Medium gilt nach Nyquist:

$R_{\max} = 2B \cdot \log_2(M)$

| Symbol       | Bedeutung Deutsch                      | 中文含义      |
| ------------ | -------------------------------------- | --------- |
| $⁡R_{\max}$​ | maximale Datenrate                     | 最大数据率     |
| B            | Bandbreite                             | 带宽        |
| M            | Anzahl unterscheidbarer Signalzustände | 可区分的信号状态数 |

---

## 1.5 Herausforderungen verteilter Systeme

## 1.5 分布式系统的挑战

**Frage / 问题**

> Nennen Sie zwei Herausforderungen von verteilten Systemen im Vergleich zu Einzelsystemen!  
> 相比单机系统，分布式系统有哪些挑战？请举两个。

**Antwort / 解答**

Mögliche Antworten:

| Deutsch              | 中文     |
| -------------------- | ------ |
| Nebenläufigkeit      | 并发性    |
| Teilausfälle         | 局部故障   |
| fehlende globale Uhr | 缺少全局时钟 |
| Konsistenzprobleme   | 一致性问题  |
| Netzwerklatenz       | 网络延迟   |
| Nachrichtenverlust   | 消息丢失   |
| Sicherheit           | 安全性    |

**Kurzantwort / 简答**

```text
Teilausfälle und fehlende globale Uhr.
局部故障和缺少全局时钟。
```

---

## 1.6 IPv4-Fragmentierung erkennen

## 1.6 如何识别 IPv4 分片

**Frage / 问题**

> Woran wird bei IPv4 erkannt, dass eine PDU fragmentiert ist?  
> 在 IPv4 中如何判断一个 PDU 被分片了？

**Antwort / 解答**

An folgenden Feldern im IPv4-Header:

| Feld                | Bedeutung Deutsch                              | 中文                |
| ------------------- | ---------------------------------------------- | ----------------- |
| Identification      | gleiche ID für alle Fragmente desselben Pakets | 同一原始包的所有分片具有相同 ID |
| Flags, besonders MF | More Fragments Flag                            | MF 标志表示后面还有分片     |
| Fragment Offset     | Position des Fragments im ursprünglichen Paket | 分片在原始包中的偏移        |

**Kurzantwort / 简答**

```text
Am MF-Flag und/oder am Fragment Offset.
通过 MF 标志和/或 Fragment Offset 字段识别。
```

---

## 1.7 Von SDU zu PDU

## 1.7 从 SDU 到 PDU

**Frage / 问题**

> Was muss in einer Schicht zu einer SDU hinzugefügt werden, um eine PDU zu erhalten?  
> 某一层中，SDU 加上什么会变成 PDU？

**Antwort / 解答**

```text
Protokollkontrollinformationen, also Header und gegebenenfalls Trailer.
协议控制信息，即头部，有时还包括尾部。
```

**Schema / 示意图**

```text
+----------------------+----------------------+----------------------+
| Header               | SDU                  | Trailer              |
| 头部                 | 服务数据单元          | 尾部                 |
+----------------------+----------------------+----------------------+

= PDU / 协议数据单元
```

---

# 2. ISO/OSI-Modell

# 2. ISO/OSI 模型

## 2.1 System- und Dienstschnitt

## 2.1 系统接口与服务接口

**Frage / 问题**

> Zeichnen Sie in folgende Grafik den System- und Dienstschnitt ein und beschriften Sie diese.  
> 在图中画出系统接口和服务接口并标注。

## Lösung / 解答

**Begriffe / 概念**

| Deutsch       | 中文   | Bedeutung                                               |
| ------------- | ---- | ------------------------------------------------------- |
| Systemschnitt | 系统接口 | zwischen zwei Systemen, also zwischen Client und Server |
| Dienstschnitt | 服务接口 | zwischen zwei Schichten innerhalb desselben Systems     |

**ASCII-Grafik / ASCII 图**

```text
                 Systemschnitt / 系统接口
        <---------------------------------------->
        Kommunikation gleicher Schichten / 同层通信

+-------------------+                  +-------------------+
|      Client       |                  |      Server       |
|                   |                  |                   |
|   Schicht i       |                  |   Schicht i       |
|   第 i 层         |                  |   第 i 层         |
|        |          |                  |          |        |
|========+==========|                  |==========+========|
| Dienstschnitt     |                  | Dienstschnitt     |
| 服务接口          |                  | 服务接口          |
|        |          |                  |          |        |
|   Schicht i+1     |                  |   Schicht i-1     |
|   第 i+1 层       |                  |   第 i-1 层       |
|                   |                  |                   |
+---------+---------+                  +---------+---------+
          |                                      |
          +--------------------------------------+
               physische/logische Verbindung
               物理或逻辑连接
```

**Hinweis / 说明**

- **Dienstschnitt / 服务接口**：vertikal innerhalb eines Systems zwischen benachbarten Schichten.  
    在同一个系统内部，相邻层之间的垂直接口。
- **Systemschnitt / 系统接口**：horizontal zwischen verschiedenen Systemen auf gleicher Ebene.  
    在不同系统之间，同层实体之间的水平接口。

---

# 3. IPv6

# 3. IPv6

## 3.1 IPv6-Adresse maximal verkürzen

## 3.1 IPv6 地址最大缩写

**Frage / 问题**

> Verkürzen Sie folgende IPv6-Adresse maximal:  
> 最大程度缩写以下 IPv6 地址：

```text
1337:0000:0000:0000:1000:0000:0000:0001
```

**Antwort / 解答**

```text
1337::1000:0:0:1
```

**Erklärung / 说明**

Die längste Nullsequenz ist:

```text
0000:0000:0000
```

Diese wird durch `::` ersetzt.  
最长的连续零段是三个 `0000`，用 `::` 替代。

Andere führende Nullen werden entfernt:

```text
0001 -> 1
0000 -> 0
```

---

## 3.2 IPv6-Subnetting

## 3.2 IPv6 子网划分

**Gegeben / 已知**

Aus dem Bild ist der Anfang der Adresse erkennbar als:

```text
DE94:DEF0:0000:0000:0000:0000:0000:0000/?
```

Die rote Markierung macht die genaue Präfixlänge unsicher. In solchen Aufgaben ist sehr wahrscheinlich ein Präfix wie `/28` oder `/32` gemeint.  
图中红色部分表示原始资料不确定，因此前缀长度不完全清楚。此类题通常给出类似 `/28` 或 `/32` 的前缀。

Da gefragt wird, wie viele Bits für 4 Subnetze benötigt werden:

$4 = 2^2$

Daher werden **2 zusätzliche Bits** benötigt.

**Antwort i / 解答 i**

```text
2 zusätzliche Bits für die Netz-ID.
网络 ID 需要额外 2 位。
```

---

### Fall A: Wenn das Ausgangsnetz `/28` ist

### 情况 A：如果原网络是 `/28`

Neue Präfixlänge:

28 + 2 = 30

```text
DE94:DEC0::/30
DE94:DF00::/30
DE94:DF40::/30
DE94:DF80::/30
```

**Achtung / 注意**

Bei `/28` ist `DE94:DEF0::/28` keine kanonische Netzadresse. Die kanonische Netzadresse wäre:

```text
DE94:DEC0::/28
```

---

### Fall B: Wenn das Ausgangsnetz `/32` ist

### 情况 B：如果原网络是 `/32`

Neue Präfixlänge:

32 + 2 = 34

```text
DE94:DEF0::/34
DE94:DEF0:4000::/34
DE94:DEF0:8000::/34
DE94:DEF0:C000::/34
```

---

**Kurzantwort / 简答**

```text
Für 4 Subnetze braucht man 2 zusätzliche Bits.
Bei Ausgangspräfix /p entstehen Subnetze mit /p+2.
分成 4 个子网需要额外 2 位。
如果原前缀是 /p，则新子网前缀为 /p+2。
```

---

# 4. Unbekanntes Protokoll

# 4. 未知协议

**Gegebene PDU / 给定 PDU**

```text
+---------------------------+------------------------------+
| Sequenznummer 4 Bit       | Quittungsnummer 4 Bit        |
| 序列号 4 位               | 确认号 4 位                  |
+---------------------------+------------------------------+
|                                                          |
|                       Nutzdaten                          |
|                       用户数据                           |
|                                                          |
+----------------------------------------------------------+
```

---

## 4.1 Warum Sequenznummer zur Nummerierung der PDUs?

## 4.1 为什么认为序列号用于给 PDU 编号？

**Frage / 问题**

> Warum ist es sinnvoll davon auszugehen, dass die Sequenznummer dazu benutzt wird, die PDUs zu nummerieren?  
> 为什么可以合理认为序列号用于给 PDU 编号？

**Antwort / 解答**

Weil die PDU ein explizites Feld namens **Sequenznummer** besitzt. Ein solches Feld dient typischerweise dazu, gesendete PDUs zu nummerieren, Reihenfolge zu erkennen und Duplikate zu identifizieren.  
因为 PDU 中有明确的 **Sequenznummer / 序列号** 字段。该字段通常用于给 PDU 编号、识别顺序以及检测重复包。

```text
Sequenznummer vorhanden -> Nummerierung der PDUs plausibel.
存在序列号字段 -> 合理推断用于 PDU 编号。
```

---

## 4.2 Netz mit 4 Teilnehmern?

## 4.2 能否用于 4 个参与者的网络？

**Frage / 问题**

> Kann man dieses Protokoll für ein Netz mit 4 Teilnehmern verwenden?  
> 该协议能用于有 4 个参与者的网络吗？

**Antwort / 解答**

```text
Nein, nicht sinnvoll bzw. nicht eindeutig.
不能，或者说不能唯一寻址。
```

**Begründung / 理由**

In der PDU gibt es keine Felder für:

```text
Quelladresse / 源地址
Zieladresse / 目的地址
```

Ohne Adressfelder kann bei 4 Teilnehmern nicht eindeutig angegeben werden, wer Sender und Empfänger ist.  
由于没有源地址和目的地址字段，无法在 4 个参与者中唯一标识通信双方。

---

## 4.3 Können Fehler in den Nutzdaten erkannt werden?

## 4.3 能否检测用户数据中的错误？

**Frage / 问题**

> Können Fehler in den Nutzdaten erkannt werden?  
> 能否检测 Nutzdaten 中的错误？

**Antwort / 解答**

```text
Nein.
不能。
```

**Begründung / 理由**

Es gibt kein Prüfsummenfeld, kein CRC-Feld und keine sonstige Fehlererkennungsinformation in der PDU.  
PDU 中没有校验和字段、CRC 字段或其他错误检测字段。

---

## 4.4 Maximale Sendefenstergröße

## 4.4 最大发送窗口大小

**Frage / 问题**

> Wie groß darf das Sendefenster maximal sein, damit Duplikate zuverlässig erkannt werden können?  
> 为了可靠识别重复包，发送窗口最大可以多大？

**Gegeben / 已知**

Sequenznummer:

```text
4 Bit
```

Anzahl verschiedener Sequenznummern:

$2^4 = 16$

Für Sliding Window mit zuverlässiger Duplikaterkennung gilt typischerweise:

$W_{\max} = \frac{2^k}{2}$​

Bei k=4:

$W_{\max} = \frac{16}{2} = 8$

**Antwort / 解答**

```text
8 PDUs
8 个 PDU
```

---

# 5. Generatorpolynom / CRC

# 5. 生成多项式 / CRC

## 5.1 Länge der CRC-Prüfsumme

## 5.1 CRC 校验和长度

**Frage / 问题**

> Es sei folgendes Generatorpolynom gegeben:  
> 给定生成多项式：

$G = x^3 + 1$

> Wie lang ist die CRC-Prüfsumme mit diesem Generatorpolynom?  
> 该生成多项式对应的 CRC 校验和长度是多少？

**Antwort / 解答**

Der Grad des Generatorpolynoms ist:

3

Daher ist die CRC-Prüfsumme:

```text
3 Bit
```

**Generator als Bitfolge / 生成多项式对应比特串**

$x^3 + 1 \Rightarrow 1001$

---

## 5.2 Fehlererkennung bei empfangener Bitfolge

## 5.2 接收端错误检测

**Gegeben / 已知**

Empfangene Bitfolge:

```text
0011001
```

Generator:

```text
1001
```

Die im Bild gezeigte Rechnung:

```text
0011000 : 1001 = ...
 1000
 1001
 ----
    1
```

**Frage / 问题**

> Kennzeichnen Sie die Stelle in der Rechnung, an der bereits erkannt werden konnte, dass ein Fehler aufgetreten ist.  
> 标出计算中已经可以看出发生错误的位置，并说明原因。

**Antwort / 解答**

Ein Fehler wird erkannt, sobald am Ende der CRC-Division ein von `000` verschiedener Rest übrig bleibt.

Hier bleibt als Rest:

```text
1
```

bzw. mit 3 Bit geschrieben:

```text
001
```

Da der Rest nicht Null ist:

```text
Rest != 000 -> Fehler erkannt
余数不为 000 -> 检测到错误
```

**Markierung / 标记**

```text
0011001 : 1001

...
Rest = 001  <--- Fehler hier erkennbar
              此处可识别错误
```

**Begründung / 理由**

Eine korrekte CRC-codierte Bitfolge muss durch das Generatorpolynom ohne Rest teilbar sein.  
正确的 CRC 编码比特串必须能被生成多项式整除，余数应为 0。

---

## 5.3 CRC-Prüfsumme von `110001`

## 5.3 计算 `110001` 的 CRC 校验和

**Gegeben / 已知**

Nachricht:

```text
110001
```

Generator:

```text
1001
```

Grad:

```text
3
```

Daher werden 3 Nullen angehängt:

```text
110001000
```

Modulo-2-Division:

```text
110001000 ÷ 1001
```

Rechnung:

```text
110001000
1001
----
010101000
 1001
 ----
000111000
   1001
   ----
000011100
    1001
    ----
000001010
      1001
      ----
000000011
```

Der Rest sind die letzten 3 Bit:

```text
011
```

**Antwort / 解答**

```text
CRC-Prüfsumme = 011
CRC 校验和 = 011
```

Zu sendende Bitfolge:

```text
110001011
```

**Stelle der Prüfsumme / 校验和位置**

```text
110001 011
Nachricht CRC
消息    校验和
```

---

# 6. PPP

# 6. PPP 协议

## 6.1 Beide Unterprotokolle

## 6.1 两个子协议

**Frage / 问题**

> Nennen Sie die beiden Unterprotokolle von PPP.  
> 请写出 PPP 的两个子协议。

**Antwort / 解答**

```text
LCP und NCP
```

|Abkürzung|Deutsch|中文|
|---|---|---|
|LCP|Link Control Protocol|链路控制协议|
|NCP|Network Control Protocol|网络控制协议|

---

## 6.2 Was macht LCP?

## 6.2 LCP 做什么？

**Antwort / 解答**

LCP ist zuständig für Aufbau, Konfiguration, Test und Abbau der PPP-Verbindung.  
LCP 负责 PPP 链路的建立、配置、测试和释放。

```text
LCP: Link aufbauen, konfigurieren, prüfen, abbauen.
LCP：建立、配置、检测、关闭链路。
```

---

## 6.3 Was macht NCP?

## 6.3 NCP 做什么？

**Antwort / 解答**

NCP konfiguriert die Protokolle der Vermittlungsschicht, z. B. IPv4 oder IPv6 über PPP.  
NCP 用于配置网络层协议，例如通过 PPP 承载 IPv4 或 IPv6。

```text
NCP: Netzwerkschicht-Protokolle konfigurieren.
NCP：配置网络层协议。
```

---

# 7. Modulation von analogen Daten

# 7. 模拟数据调制

Gegeben:

s(t)=A⋅sin⁡(2⋅π⋅f⋅t+ω)s(t) = A \cdot \sin(2 \cdot \pi \cdot f \cdot t + \omega)s(t)=A⋅sin(2⋅π⋅f⋅t+ω)

---

## 7.1 Drei veränderbare Größen

## 7.1 调制可改变的三个量

|Größe|Deutsch|中文|Bedeutung|
|---|---|---|---|
|AAA|Amplitude|振幅|Höhe / Stärke des Signals，信号强度|
|fff|Frequenz|频率|Schwingungen pro Sekunde，每秒振荡次数|
|ω\omegaω|Phase|相位|zeitliche Verschiebung der Schwingung，波形的时间偏移|

**Antwort / 解答**

```text
Größe 1: Amplitude A
Größe 2: Frequenz f
Größe 3: Phase ω
```

---

## 7.2 Zuordnung der Diagramme

## 7.2 图像对应的调制量

Die drei Bilder zeigen ungefähr:

```text
Bild 1: unterschiedliche Höhe der Ausschläge
Bild 2: unterschiedliche Breite / Periodendauer
Bild 3: zeitliche Verschiebung / Lage der Welle
```

**Zuordnung / 对应关系**

|Diagramm|Veränderter Term|Deutsch|中文|
|--:|---|---|---|
|1|AAA|Amplitude|振幅|
|2|fff|Frequenz|频率|
|3|ω\omegaω|Phase|相位|

**ASCII-Skizze / ASCII 示意**

```text
1) Amplitude A / 振幅变化

Signal
  ^
  |       /\          /\
  |      /  \        /  \
  | /\  /    \  /\  /    \
--+--------------------------> Zeit
        hohe Amplitude
        高振幅


2) Frequenz f / 频率变化

Signal
  ^
  | /\  /\       /--------\
  |/  \/  \     /          \
--+--------------------------> Zeit
   hohe f       niedrige f
   高频          低频


3) Phase ω / 相位变化

Signal
  ^
  |   /\    /\     /\
  |  /  \  /  \   /  \
--+--------------------------> Zeit
    verschobene Lage
    相位偏移
```

---

## 7.3 Bits pro Signalschritt bei 8 Signalzuständen

## 7.3 8 个信号状态每个信号步传几位？

**Frage / 问题**

> Wie viele Bits werden pro Signalschritt übertragen, wenn 8 Signalzustände unterscheidbar sind?  
> 如果有 8 个可区分信号状态，每个信号步传输多少 bit？

**Rechnung / 计算**

log⁡2(8)=3\log_2(8) = 3log2​(8)=3

**Antwort / 解答**

```text
3 Bit pro Signalschritt.
每个信号步 3 bit。
```

---

## 7.4 Bandbreite und maximale Übertragungsrate

## 7.4 带宽与最大传输速率

**Gegeben / 已知**

Grenzfrequenzen:

```text
f_min = 50 MHz
f_max = 100 MHz
```

---

### 7.4.1 Bandbreite

### 7.4.1 带宽

B=fmax⁡−fmin⁡B = f_{\max} - f_{\min}B=fmax​−fmin​

B=100 MHz−50 MHz=50 MHzB = 100\,\text{MHz} - 50\,\text{MHz} = 50\,\text{MHz}B=100MHz−50MHz=50MHz

**Antwort / 解答**

```text
50 MHz
```

---

### 7.4.2 Maximale Übertragungsrate mit Binärcodierung

### 7.4.2 二进制编码下最大传输速率

Für ein rauschfreies ideales Medium nach Nyquist:

$R_{\max} = 2B \cdot \log_2(M)$

Bei Binärcodierung:

M=2

$\log_2(2) = 1$

Also:

$R_{\max} = 2 \cdot 50\,\text{MHz} \cdot 1$

$R_{\max} = 100\,\text{Mbit/s}$

**Antwort / 解答**

```text
100 Mbit/s
```

---

# 8. TCP

# 8. TCP

## 8.1 TCP-Flags im Ablaufdiagramm

## 8.1 时序图中的 TCP 标志

**Gegeben / 已知**

```text
A initial sequence number = 1000
B initial sequence number = 4001
A sends 500 Byte to B.
A 的初始序列号 = 1000
B 的初始序列号 = 4001
A 向 B 发送 500 Byte 数据。
```

Zusätzlich / 另外：

```text
Verbindungsabbau wird von A initiiert.
连接关闭由 A 发起。
```

Das Diagramm zeigt von oben nach unten ungefähr:

1. Verbindungsaufbau
2. Nachricht von A an B
3. Quittung
4. Verbindungsabbau

**TCP-Sequenzdiagramm / TCP 时序图**

```text
A                                                        B
|                                                        |
|  SYN, Seq=1000 --------------------------------------> |
|                                                        |
|  <---------------------------- SYN,ACK, Seq=4001, Ack=1001 |
|                                                        |
|  ACK, Seq=1001, Ack=4002 ----------------------------> |
|                                                        |
|  PSH,ACK, Seq=1001, Ack=4002, Len=500 --------------> |
|                                                        |
|  <---------------------------- ACK, Seq=4002, Ack=1501 |
|                                                        |
|  FIN,ACK, Seq=1501, Ack=4002 ------------------------> |
|                                                        |
|  <---------------------------- ACK, Seq=4002, Ack=1502 |
|                                                        |
|  <---------------------------- FIN,ACK, Seq=4002, Ack=1502 |
|                                                        |
|  ACK, Seq=1502, Ack=4003 ----------------------------> |
|                                                        |
```

**Nur Flags für die Kästchen / 只填标志**

Von oben nach unten:

|Nr.|Richtung|Flags Deutsch|中文|
|--:|---|---|---|
|1|A → B|SYN|连接建立请求|
|2|B → A|SYN, ACK|同意建立并确认|
|3|A → B|ACK|确认|
|4|A → B|PSH, ACK|数据传输|
|5|B → A|ACK|确认数据|
|6|A → B|FIN, ACK|A 请求关闭|
|7|B → A|ACK|B 确认关闭请求|
|8|B → A|FIN, ACK|B 请求关闭|
|9|A → B|ACK|A 最终确认|

---

## 8.2 Sequenznummer beim vierten Pfeil von A aus

## 8.2 从 A 发出的第 4 个箭头的序列号

**Frage / 问题**

> Angenommen die Quittungsnummer beim zweiten Pfeil von B aus ist yyy.  
> Was ist dann die Sequenznummer beim vierten Pfeil von A aus?  
> 假设从 B 发出的第二个箭头中的确认号为 yyy，那么从 A 发出的第四个箭头的序列号是多少？

Interpretation:

- Zweiter Pfeil von B aus ist ACK auf die 500-Byte-Nachricht.
- Diese ACK-Nummer ist yyy.
- Damit bestätigt B alle Bytes bis y−1y-1y−1.
- A hat danach als nächste Sequenznummer yyy.
- Der vierte Pfeil von A aus ist dann der FIN von A.

**Antwort / 解答**

```text
Sequenznummer = y
序列号 = y
```

Wenn man mit den konkreten Zahlen rechnet:

```text
A startet mit Seq = 1000
SYN verbraucht 1 -> nächste Seq = 1001
Daten 500 Byte -> nächste Seq = 1501
```

Also:

```text
y = 1501
Sequenznummer des FIN von A = 1501
```

---

## 8.3 Um wie viel erhöht sich die Sequenznummer insgesamt bei B?

## 8.3 B 的序列号总共增加多少？

B sendet:

|Aktion|Erhöhung|
|---|--:|
|SYN|+1|
|reine ACKs|+0|
|FIN|+1|

Keine Nutzdaten von B angegeben.  
题目没有说明 B 发送用户数据。

1+1=21 + 1 = 21+1=2

**Antwort / 解答**

```text
2
```

B startet bei:

```text
4001
```

Ende:

```text
4003
```

---

## 8.4 Protokolle und Schichten

## 8.4 协议与层

**Frage / 问题**

> Benennen Sie zu jedem aufgeführten Protokoll den Schichtnamen.  
> 给出每个协议所在的层名。

|Protokoll|OSI-Schicht Deutsch|中文|
|---|---|---|
|HTTP|Anwendungsschicht|应用层|
|ICMP|Vermittlungsschicht|网络层|
|ARP|Sicherungsschicht / zwischen Schicht 2 und 3|数据链路层 / 介于 2、3 层之间|
|TCP|Transportschicht|传输层|

**Hinweis / 说明**

ARP 在教材中常被视为第 2 层协议，或介于第 2 层和第 3 层之间，因为它把 IP 地址解析为 MAC 地址。

---

# 9. Fragmentierung

# 9. 分片

## 9.1 Gegebenes Netz

## 9.1 给定网络

```text
+----+      MTU 1000       +----+      MTU 280       +----+      MTU 500       +----+
| E1 |---------------------| R1 |--------------------| R2 |--------------------| E2 |
+----+                     +----+                    +----+                    +----+
```

---

## 9.2 Wo wird fragmentiert?

## 9.2 在哪里分片？

**Frage / 问题**

> Nennen Sie die Komponente(n), an denen fragmentiert wird.  
> 请指出在哪些组件处发生分片。

**Antwort / 解答**

Bei IPv6 fragmentieren Router nicht. Fragmentierung erfolgt nur am Sender.  
IPv6 中路由器不进行分片，只有源主机进行分片。

```text
E1
```

**Erklärung / 说明**

Der kleinste MTU-Wert auf dem Pfad ist:

```text
280 Byte
```

E1 muss die Pakete also so fragmentieren, dass sie über diese MTU passen.  
路径最小 MTU 是 280 Byte，因此 E1 必须按 280 Byte 的限制分片。

---

## 9.3 Fragmentierung von 600 Byte Nutzdaten

## 9.3 600 Byte 用户数据的 IPv6 分片

**Gegeben / 已知**

```text
IPv6-Nutzdaten: 600 Byte
Original IPv6 Header: 40 Byte
IPv6 Header inkl. Extension Header und Fragmentation Header: 48 Byte
Kleinste MTU: 280 Byte
```

Originalpaket:

```text
+------------------+--------------------------------------+
| Header 40 Byte   | Nutzdaten 600 Byte                   |
| 头部 40 字节     | 用户数据 600 字节                    |
+------------------+--------------------------------------+
```

Für jedes Fragment gilt:

```text
maximale Fragmentgröße = 280 Byte
Header pro Fragment = 48 Byte
```

Maximale Nutzdaten pro Fragment:

280 - 48 = 232

Da bei IPv6 Fragmentdaten außer beim letzten Fragment ein Vielfaches von 8 Byte sein müssen:

232 / 8 = 29

232 ist gültig.

Aufteilung der 600 Byte:

```text
600 = 232 + 232 + 136
```

---

## 9.4 Fragmente bei E2

## 9.4 E2 收到的分片

|Fragmentnr.|Fragment Deutsch|中文|
|--:|---|---|
|1|Header 48 Byte + Nutzdaten 232 Byte = 280 Byte|头部 48 字节 + 数据 232 字节 = 280 字节|
|2|Header 48 Byte + Nutzdaten 232 Byte = 280 Byte|头部 48 字节 + 数据 232 字节 = 280 字节|
|3|Header 48 Byte + Nutzdaten 136 Byte = 184 Byte|头部 48 字节 + 数据 136 字节 = 184 字节|

**Grafische Darstellung / 图形表示**

```text
Fragment 1:
+----------------------+------------------------------------------------+
| Header 48 Byte       | Nutzdaten 232 Byte                            |
| 头部 48 字节         | 数据 232 字节                                 |
+----------------------+------------------------------------------------+
Gesamt: 280 Byte


Fragment 2:
+----------------------+------------------------------------------------+
| Header 48 Byte       | Nutzdaten 232 Byte                            |
| 头部 48 字节         | 数据 232 字节                                 |
+----------------------+------------------------------------------------+
Gesamt: 280 Byte


Fragment 3:
+----------------------+----------------------------+
| Header 48 Byte       | Nutzdaten 136 Byte        |
| 头部 48 字节         | 数据 136 字节             |
+----------------------+----------------------------+
Gesamt: 184 Byte
```

---

# 10. Routertabellen

# 10. 路由表

## 10.1 Gegebenes Routernetz

## 10.1 给定路由器网络

Aus der Abbildung ergibt sich ungefähr folgendes Netz:

```text
          B
          |
          |
A---------C---------E
|         |         |
|         |         |
+---------D---------+
```

Kanten / 连接关系：

```text
A-C
C-B
C-E
C-D
A-D
D-E
```

Alle Kanten werden als gleich gewichtet angenommen.  
假设所有链路代价相同。

---

## 10.2 Kürzester Quell-Senken-Baum für B

## 10.2 以 B 为源的最短路径树

Von B aus:

|Ziel|Kürzester Weg|Distanz|
|---|---|--:|
|C|B-C|1|
|A|B-C-A|2|
|D|B-C-D|2|
|E|B-C-E|2|

Ein korrekter kürzester Baum ist:

```text
          B
          |
          C
        / | \
       A  D  E
```

**Deutsch / 中文说明**

- B ist direkt mit C verbunden.  
    B 直接连接 C。
- Alle anderen Router sind von B aus über C erreichbar.  
    其他路由器均可通过 C 到达。

---

## 10.3 Netz mit Subnetzen

## 10.3 含子网的网络

Die Abbildung zeigt:

```text
                 Internet
                    |
                    |
                    B
                    |
                    |
Subnetz 2           C--------------E-----------Subnetz 3
   |                |              |
   |                |              |
   A----------------+              |
   |                               |
   |                               |
   +---------------D---------------+
                   |
                   |
                Subnetz 4

Subnetz 1 hängt an A.
Subnetz 2 hängt an A.
Subnetz 3 hängt an E.
Subnetz 4 hängt an D.
Internet hängt an B.
```

Kompakter:

```text
                 [Internet]
                     |
                     B
                     |
                     C
                   / | \
                  A  D  E
                 /|  |   \
        [S1] [S2] [S4] [S3]
```

---

## 10.4 Routingtabelle von C für fünf Subnetze

## 10.4 C 到五个子网的路由表

Gesucht:

> Geben Sie die vollständige Routingtabelle von C für alle fünf Subnetze an.  
> 给出 C 到所有五个子网的完整路由表。

Subnets:

```text
Subnetz 1
Subnetz 2
Subnetz 3
Subnetz 4
Internet
```

Aus Sicht von C:

|Subnetz|Erreichbar über Deutsch|中文|
|---|---|---|
|Subnetz 1|A|通过 A|
|Subnetz 2|A|通过 A|
|Subnetz 3|E|通过 E|
|Subnetz 4|D|通过 D|
|Internet|B|通过 B|

**Antwort / 解答**

```text
+-----------+----------------+
| Subnetz   | erreichbar über|
+-----------+----------------+
| Subnetz 1 | A              |
| Subnetz 2 | A              |
| Subnetz 3 | E              |
| Subnetz 4 | D              |
| Internet  | B              |
+-----------+----------------+
```

---

# 11. DNS-Anfragen

# 11. DNS 查询

## 11.1 Gegebener dig-Auszug

## 11.1 给定 dig 输出

Der Ausschnitt zeigt ungefähr:

```text
dig +trace mail.nm.ifi.lmu.de

.                    NS   root-servers
de.                  NS   C.DE.NET, ...
lmu.de.              NS   dns1.lrz-muenchen.de, ...
mail.nm.ifi.lmu.de.  CNAME pcheger0.nm.ifi.lmu.de.
pcheger0.nm.ifi.lmu.de. A 141.84.218.30
nm.ifi.lmu.de.       NS ...
```

---

## 11.2 Wie viele DNS-Anfragen wurden ausgeführt?

## 11.2 执行了多少次 DNS 查询？

Im `dig +trace` sieht man Antworten von:

1. lokalem Resolver / Root-Hinweis bzw. Root-Server-Liste
2. Root-Server für `.de`
3. `.de`-Server für `lmu.de`
4. `lmu.de`-Server für `nm.ifi.lmu.de` und Zielhost

Im gezeigten Auszug stehen vier `Received ... from ...`-Zeilen.

**Antwort / 解答**

```text
4 DNS-Anfragen
4 次 DNS 查询
```

---

## 11.3 IP von `nm.ifi.lmu.de`

## 11.3 `nm.ifi.lmu.de` 的 IP

Im Auszug ist sichtbar:

```text
pcheger0.nm.ifi.lmu.de.  A  141.84.218.30
mail.nm.ifi.lmu.de.      CNAME pcheger0.nm.ifi.lmu.de.
```

Daher hat `mail.nm.ifi.lmu.de` indirekt die IP:

```text
141.84.218.30
```

Falls wirklich nach `nm.ifi.lmu.de` gefragt ist, zeigt der Ausschnitt vor allem Nameserver-Einträge für diese Zone, aber keinen eindeutigen A-Record für `nm.ifi.lmu.de` selbst.  
如果严格问的是 `nm.ifi.lmu.de` 本身，图中主要显示其 NS 记录，不一定显示它自己的 A 记录。

**Wahrscheinliche Prüfungsantwort / 可能的考试答案**

```text
141.84.218.30
```

---

## 11.4 Iterativ oder rekursiv?

## 11.4 迭代还是递归？

**Frage / 问题**

> War die Anfrage iterativ oder rekursiv? Begründen Sie!  
> 查询是迭代还是递归？请说明。

**Antwort / 解答**

```text
Iterativ.
迭代查询。
```

**Begründung / 理由**

Bei `dig +trace` werden die zuständigen Nameserver schrittweise abgefragt:

```text
Root -> de. -> lmu.de. -> nm.ifi.lmu.de.
```

Jeder Server liefert entweder eine Antwort oder einen Verweis auf den nächsten zuständigen Nameserver.  
`dig +trace` 会一步步询问 root、TLD、权威服务器，每一步返回下一层服务器信息，因此是迭代解析。

---

## 11.5 Kann man mit diesen Informationen eine E-Mail an `post@ifi.lmu.de` schicken?

## 11.5 能否仅凭这些信息给 `post@ifi.lmu.de` 发邮件？

**Antwort / 解答**

```text
Nein, nicht sicher.
不能确定，通常不够。
```

**Begründung / 理由**

Für das Ausliefern einer E-Mail an:

```text
post@ifi.lmu.de
```

braucht man den MX-Record der Domain:

```text
ifi.lmu.de
```

Im gegebenen Ausschnitt geht es aber um:

```text
mail.nm.ifi.lmu.de
```

und es sind keine MX-Records für `ifi.lmu.de` angegeben.

**中文说明**

要投递 `post@ifi.lmu.de`，需要查询 `ifi.lmu.de` 的 MX 记录。图中信息主要是 `mail.nm.ifi.lmu.de` 的 CNAME 和 A 记录，因此不能仅凭这些信息确认能否投递。

---

## 11.6 DNS-Funktionstypen

## 11.6 DNS 记录类型的函数形式

Gegeben:

```text
A: f_A(Hostname) = IPv4-Adresse
A 记录：主机名 -> IPv4 地址
```

### i. NS

```text
NS: f_NS(Domainname) = zuständiger Nameserver
NS：域名 -> 权威名称服务器
```

Beispiel:

```text
f_NS(lmu.de) = dns1.lrz-muenchen.de
```

---

### ii. CNAME

```text
CNAME: f_CNAME(Aliasname) = kanonischer Hostname
CNAME：别名 -> 规范主机名
```

Beispiel:

```text
f_CNAME(mail.nm.ifi.lmu.de) = pcheger0.nm.ifi.lmu.de
```

---

### iii. MX

```text
MX: f_MX(Domainname) = Mailserver der Domain
MX：域名 -> 邮件服务器
```

Beispiel:

```text
f_MX(ifi.lmu.de) = mailserver.ifi.lmu.de
```

---

# 12. HTTP-Anfragen

# 12. HTTP 请求

## 12.1 Szenario

## 12.1 场景

Gegeben ist ein Netz mit Client, HTTP-Server, DNS-Server und einer Komponente X.

```text
       +-------------+
       | HTTP Server |
       +-------------+
              |
              |
              o-------------+-------------+----------+
              |             |             |          |
       +-------------+   +-----+      +--------+     |
       | DNS Server  |   |  X  |------| Client |     |
       +-------------+   +-----+      +--------+
```

Etwas übersichtlicher:

```text
+-------------+        +-----+        +--------+
| HTTP Server |--------|  X  |--------| Client |
+-------------+        +-----+        +--------+
        |
        |
+-------------+
| DNS Server  |
+-------------+
```

---

## 12.2 Fall 1: X ist eine Bridge

## 12.2 情况 1：X 是桥接器

### 12.2.1 MAC-Adresse für DNS-Anfrage

### 12.2.1 DNS 查询发往哪个 MAC 地址？

**Frage / 问题**

> An welche MAC-Adresse schickt der Client seine DNS-Anfrage?  
> 客户端把 DNS 查询发送到哪个 MAC 地址？

**Antwort / 解答**

Da eine Bridge auf Schicht 2 arbeitet und alle Geräte im selben IP-Netz liegen, sendet der Client den Ethernet-Frame direkt an die MAC-Adresse des DNS-Servers.

```text
An die MAC-Adresse des DNS-Servers.
发送到 DNS 服务器的 MAC 地址。
```

**Hinweis / 说明**

Falls der Client die MAC-Adresse noch nicht kennt, benutzt er zuerst ARP.  
如果客户端还不知道 DNS 服务器的 MAC 地址，会先使用 ARP 查询。

---

### 12.2.2 IP-Adresse in der HTTP-Anfrage

### 12.2.2 HTTP 请求中的 IP 地址

**Frage / 问题**

> Welche IP-Adresse steht in der HTTP-Anfrage?  
> HTTP 请求中写的是哪个 IP 地址？

**Antwort / 解答**

Auf IP-Ebene ist die Zieladresse die IP-Adresse des HTTP-Servers.

```text
Ziel-IP = IP-Adresse des HTTP-Servers.
目的 IP = HTTP 服务器的 IP 地址。
```

**Erklärung / 说明**

Eine Bridge verändert IP-Adressen nicht.  
桥接器不会修改 IP 地址。

---

## 12.3 Fall 2: X ist ein Router

## 12.3 情况 2：X 是路由器

Gegeben:

```text
Router-Port 0: 92.16.0.1
Router-Port 1: 92.17.0.1
```

Das bedeutet: Client und Server liegen in unterschiedlichen IP-Netzen.  
这意味着客户端和服务器位于不同 IP 网络中。

---

### 12.3.1 MAC-Adresse für HTTP-Anfrage

### 12.3.1 HTTP 请求发往哪个 MAC 地址？

**Frage / 问题**

> An welche MAC-Adresse schickt der Client seine HTTP-Anfrage?  
> 客户端把 HTTP 请求发送到哪个 MAC 地址？

**Antwort / 解答**

Der Client sendet den Ethernet-Frame an die MAC-Adresse seines Default Gateways, also an den Router-Port im lokalen Netz des Clients.

```text
An die MAC-Adresse des Router-Ports auf der Client-Seite.
发送到客户端所在网络中路由器接口的 MAC 地址。
```

Wenn der Client im Netz von Port 1 liegt:

```text
MAC-Adresse von Router-Port 1.
```

Wenn der Client im Netz von Port 0 liegt:

```text
MAC-Adresse von Router-Port 0.
```

Die IP-Zieladresse bleibt trotzdem die IP-Adresse des HTTP-Servers.  
但 IP 目的地址仍然是 HTTP 服务器的 IP。

---

### 12.3.2 IP-Adresse, an die der Router die DNS-Anfrage schickt

### 12.3.2 路由器把 DNS 查询发往哪个 IP 地址？

**Frage / 问题**

> An welche IP-Adresse schickt der Router die DNS-Anfrage?  
> 路由器把 DNS 查询发往哪个 IP 地址？

**Antwort / 解答**

Ein Router leitet IP-Pakete anhand der Ziel-IP-Adresse weiter und verändert die Ziel-IP normalerweise nicht.

```text
An die IP-Adresse des DNS-Servers.
发往 DNS 服务器的 IP 地址。
```

**Erklärung / 说明**

Der Router ändert beim normalen Routing nur die Link-Layer-Adressen pro Hop, nicht aber die IP-Zieladresse.  
普通路由转发中，路由器会更改每一跳的 MAC 地址，但不会改变 IP 目的地址。

---

## 12.4 Warum braucht ein Router IP-Adressen an seinen Ports?

## 12.4 为什么路由器接口需要 IP 地址，而桥接器不需要？

**Frage / 问题**

> Warum braucht ein Router im Gegensatz zu einer Bridge IP-Adressen an seinen Ports?  
> 为什么路由器的端口需要 IP 地址，而桥接器的端口不需要？

**Antwort / 解答**

Ein Router arbeitet auf Schicht 3 und verbindet unterschiedliche IP-Netze. Jeder Router-Port gehört zu einem eigenen IP-Subnetz und muss dort als nächster Hop bzw. Gateway erreichbar sein. Deshalb braucht jeder Router-Port eine IP-Adresse.

Eine Bridge arbeitet dagegen auf Schicht 2. Sie verbindet Segmente desselben IP-Netzes und leitet Frames anhand von MAC-Adressen weiter. Daher braucht sie für die reine Weiterleitung keine IP-Adresse an jedem Port.

**中文解答**

路由器工作在第 3 层，用来连接不同的 IP 子网。每个路由器端口属于一个独立的 IP 子网，并且通常作为该子网中的默认网关，因此每个端口需要 IP 地址。

桥接器工作在第 2 层，只根据 MAC 地址转发帧，连接的是同一 IP 网络中的不同链路段，因此单纯转发时不需要每个端口都有 IP 地址。

**Kurzantwort / 简答**

```text
Router: Schicht 3, verbindet IP-Netze, Ports brauchen IP-Adressen.
Bridge: Schicht 2, verbindet LAN-Segmente, forwarding über MAC-Adressen.
```