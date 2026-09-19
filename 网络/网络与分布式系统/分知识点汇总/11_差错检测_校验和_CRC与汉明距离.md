# 差错检测：校验和、CRC 与汉明距离

## 知识点总结

- 奇偶校验、BCC、Internet Checksum、CRC 都是错误检测机制。
- Hamming 距离决定检测/纠错能力。
- CRC 使用模 2 除法，生成多项式最高次数决定校验位数。

## 完整题目与解答汇总

### 题目 1: 3. Fehlererkennung und -korrektur (H)

**类型：** 作业  

**来源说明：** Uebungsblatt 10, Aufgabe 3  


#### 题目中文翻译 / 中文题意

对消息 RNVS 做 7-bit/8-bit ASCII 编码，计算 Hamming 距离、二维奇偶校验、Internet Checksum 和 CRC。

#### 德文原题

```text
3. Fehlererkennung und -korrektur (H)
Ein typischer Fehler bei der Übertragung von Daten ist die Verfälschung, bei der Stellen des Bitstroms
invertiert werden. In dieser Aufgabe soll die Nachricht RNVS übertragen werden.
(a) Überführen Sie die Nachricht in eine Binärdarstellung, indem Sie die Buchstaben in 7-Bit ASCII
codieren. Hinweis: Beachten Sie die Großschreibung der Buchstaben bei der ASCII Codierung.
(b) Berechnen Sie die Hamming-Distanz aller Codewörter für die Symbole R,N,V sowie S.
(c) Die Nachricht soll in Form einer Paritätsmatrix versendet werden.
i. Schreiben Sie hierzu die 4 Codewörter untereinander und berechnen Sie für jede Zeile / Spalte
entsprechende Paritätsbits.
ii. Zeigen Sie beispielhaft, dass anhand der Paritätsmatrix 1-bit Fehler sowohl erkannt als auch
korrigiert werden können.
iii. Zeigen Sie beispielhaft, dass 2-bit Fehler erkannt, allerdings nicht korrigiert werden können.
iv. In welchen Szenarien kann ein Fehler weder erkannt noch korrigiert werden?
(d) Die Nachricht soll über das UDP Transportprotokol, das die 16-bit Internet Checksumme verwendet,
versandt werden.
i. Teilen Sie die Nachricht dazu in zwei 16-bit Segmente auf und berechnen Sie anschließend die
Internet-Checksumme nach Vorschrift. Hinweis: Führende Nullen sind notwendig, damit 8 Bit
pro Buchstaben zur Codierung verwendet werden.
ii. Zeigen Sie, dass auf Empfänger-Seite die Nachricht korrekt angekommen ist.
(e) Gegeben sei das Generatorpolynom G = x16 + x14 + x11 + x7 + x6 + x5 + 1.
Berechnen Sie die CRC-Prüfsumme über die 7-Bit ASCII codierte Darstellung der aus vier Buch-
staben bestehenden Nachricht RNVS (insgesamt 28 Bit Nutzdaten). Kennzeichnen Sie in Ihrer
Rechnung die Prüfsumme deutlich.
```

#### 解答

**3. Fehlererkennung und -korrektur / 差错检测与纠正**

**(a) ASCII / ASCII 编码**

**中文说明：** 题目 (a) 要求用 7-bit ASCII；题目 (d) 的 Internet checksum 提示需要每个字符补成 8 bit，所以表中同时列出 7-bit 和 8-bit。

| Zeichen | 7-bit ASCII | 8-bit ASCII |
|---|---|---|
| R | `1010010` | `01010010` |
| N | `1001110` | `01001110` |
| V | `1010110` | `01010110` |
| S | `1010011` | `01010011` |

**(b) Hamming-Distanzen, 7-bit / 7 位码字的汉明距离**

**中文说明：** 汉明距离就是两个等长比特串中不同位置的个数。例如 R=`1010010`，V=`1010110`，只有一个位置不同，所以距离为 1。

| Paar | Distanz |
|---|---:|
| R-N | 3 |
| R-V | 1 |
| R-S | 1 |
| N-V | 2 |
| N-S | 4 |
| V-S | 2 |

**(c) Paritaetsmatrix mit gerader Paritaet**

| Zeichen | Bits | Zeilenparitaet |
|---|---|---:|
| R | `1010010` | 1 |
| N | `1001110` | 0 |
| V | `1010110` | 0 |
| S | `1010011` | 0 |
| Spaltenparitaet | `0011001` | 1 |

**DE:** Ein 1-Bit-Fehler erzeugt genau eine falsche Zeile und eine falsche Spalte; der Schnittpunkt ist korrigierbar. Ein 2-Bit-Fehler erzeugt meist mehrere Paritaetsverletzungen und ist erkennbar, aber nicht eindeutig korrigierbar. Nicht erkannt werden koennen bestimmte Muster mit gerader Fehlerzahl pro betroffener Zeile und Spalte, z.B. vier Ecken eines Rechtecks.

**中文：** 单比特错误会导致恰好一行和一列的奇偶校验错误，交点就是出错位置，因此可以纠正。两个比特错误通常能检测出来，但无法唯一定位。若错误模式让每个受影响行和列都有偶数个错误，例如矩形四个角同时翻转，则可能检测不出来。

**(d) Internet-Checksumme, 8-bit ASCII / Internet 校验和**

**中文说明：** Internet checksum 使用 16-bit 字相加后取反。`RNVS` 用 8-bit ASCII 分成两个 16-bit 字：`RN` 和 `VS`。

16-bit-Worte:

```text
RN = 01010010 01001110 = 0x524E
VS = 01010110 01010011 = 0x5653
Summe = 0xA8A1
Checksumme = Einerkomplement = 0x575E
```

**DE:** Der Empfaenger addiert Datenworte plus Checksumme:

**中文：** 接收端把两个数据字和校验和一起按一补码加法相加：

```text
0x524E + 0x5653 + 0x575E = 0xFFFF
```

**DE:** Damit ist die Nachricht nach der Internet-Checksumme korrekt.

**中文：** 结果为全 1，即 `0xFFFF`，说明按 Internet checksum 检查时消息正确。

**(e) CRC mit G = x^16 + x^14 + x^11 + x^7 + x^6 + x^5 + 1**

**DE:** Generatorbits:

**中文：** 生成多项式对应的二进制除数为：

```text
10100100011100001
```

**DE:** Ueber `RNVS` in 7-bit ASCII (`1010010100111010101101010011`) ergibt die Division modulo 2 den Rest:

**中文：** 对 7-bit ASCII 的 `RNVS` 比特串 `1010010100111010101101010011` 进行模 2 除法，得到 16 位 CRC 余数：

```text
0100100100011011 = 0x491B
```

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 2: 4. CRC (H)

**类型：** 作业  

**来源说明：** Uebungsblatt 10, Aufgabe 4  


#### 题目中文翻译 / 中文题意

给定生成多项式 `G=x^3+1`，计算消息的 CRC 校验位，并验证接收比特串是否正确。

#### 德文原题

```text
4. CRC (H)
(a) Gegeben sei das Generatorpolynom G = x3 + 1.
i. Durch wie viele Bits wird G bei CRC repräsentiert?
ii. Es soll die Nachricht 11 00 11 CRC-geschützt übertragen werden. Berechnen Sie die zu über-
tragende Bitfolge (inkl. CRC-Prüfsumme) unter Verwendung des Generatorpolynoms G.
iii. Nehmen Sie an, dass Sie die CRC-geschützte Bitfolge 10 01 10 01 empfangen haben. Zeigen
Sie, dass die empfangene Bitfolge unter Verwendung des Generatorpolynoms G korrekt ist (inkl.
Rechnung). Markieren Sie in Ihrer Rechnung die Stelle, an der die Korrektheit sichtbar wird.
```

#### 解答

**4. CRC**

**(a i)**  
**DE:** `G = x^3 + 1` hat Grad 3 und wird durch 4 Bits dargestellt:

**中文：** `G = x^3 + 1` 的最高次数是 3，所以对应 4 位二进制多项式：

```text
1001
```

**(a ii)**  
**DE:** An die Nachricht `110011` werden drei Nullen angehaengt und dann wird durch `1001` geteilt. Rest:

**中文：** 原消息是 `110011`，因为生成多项式次数为 3，所以先在末尾补 3 个 0，再用 `1001` 做模 2 除法。余数为：

```text
101
```

**DE:** Zu uebertragende Bitfolge:

**中文：** 最终发送的比特串是原消息加 CRC 余数：

```text
110011101
```

**(a iii)**  
**DE:** Die empfangene Folge `10011001` geteilt durch `1001` ergibt Rest `000`; daran sieht man die Korrektheit.

**中文：** 接收到的 `10011001` 用同一个生成多项式 `1001` 去除，如果余数为 `000`，说明 CRC 检查通过。

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 3: Paritätsbit für 1000001 (ungerade Parität)

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Paritätsbit für 1000001 (ungerade Parität)  


#### 题目中文翻译 / 中文题意

为1000001计算奇校验位

#### 德文原题

```text
### Paritätsbit für 1000001 (ungerade Parität)

**为1000001计算奇校验位**
```

#### 解答

**Lösung / 答案：**

- 数据：1000001
- 1的个数：2（偶数）
- 奇校验需要1的总数为奇数
- **Paritätsbit = 1**
- 完整数据：**10000011**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 4: Block Check Character für 1100101|00-11101（偶校验）

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Block Check Character für 1100101|00-11101（偶校验）  


#### 题目中文翻译 / 中文题意

消息的块校验字符

#### 德文原题

```text
### Block Check Character für 1100101|00-11101（偶校验）

**消息的块校验字符**
```

#### 解答

**Lösung / 答案：**

BCC是对每一列进行偶校验：

```
1 1 0 0 1 0 1
0 0 1 1 1 0 1
--------------
1 1 1 1 0 0 0  (BCC)
```

每列1的个数加上BCC位后应为偶数。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 5: Hamming-Distanz H(C) = min({H(x,y) | x≠y ∀ x,y ∈ C})

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Hamming-Distanz H(C) = min({H(x,y) | x≠y ∀ x,y ∈ C})  


#### 题目中文翻译 / 中文题意

汉明距离公式
定义： 码字集合C的汉明距离是任意两个不同码字之间比特差异的最小值。

#### 德文原题

```text
### Hamming-Distanz H(C) = min({H(x,y) | x≠y ∀ x,y ∈ C})

**汉明距离公式**

**定义：** 码字集合C的汉明距离是任意两个不同码字之间比特差异的最小值。

---
```

#### 解答

**Hamming-Distanz H(C) = min({H(x,y) | x≠y ∀ x,y ∈ C})**

**汉明距离公式**

**定义：** 码字集合C的汉明距离是任意两个不同码字之间比特差异的最小值。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 6: (a) Was ist Hamming-Distanz H von C? Wie viele verfälschte Bits können erkannt werden?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt (a) Was ist Hamming-Distanz H von C? Wie viele verfälschte Bits können erkannt werden?  


#### 题目中文翻译 / 中文题意

C的汉明距离是多少？能检测多少位错误？
|Code C|Hamming-Distanz H(C)|可检测的位错误数|
|Paritätsprüfung / 奇偶校验|2|1|
公式：
检测d位错误需要：H(C) ≥ d + 1
纠正d位错误需要：H(C) ≥ 2d + 1

#### 德文原题

```text
### (a) Was ist Hamming-Distanz H von C? Wie viele verfälschte Bits können erkannt werden?

**C的汉明距离是多少？能检测多少位错误？**

|Code C|Hamming-Distanz H(C)|可检测的位错误数|
|---|---|---|
|Paritätsprüfung / 奇偶校验|2|1|
|BCC|2|1|

**公式：**

- 检测d位错误需要：H(C) ≥ d + 1
- 纠正d位错误需要：H(C) ≥ 2d + 1

---
```

#### 解答

**(a) Was ist Hamming-Distanz H von C? Wie viele verfälschte Bits können erkannt werden?**

**C的汉明距离是多少？能检测多少位错误？**

|Code C|Hamming-Distanz H(C)|可检测的位错误数|
|---|---|---|
|Paritätsprüfung / 奇偶校验|2|1|
|BCC|2|1|

**公式：**

- 检测d位错误需要：H(C) ≥ d + 1
- 纠正d位错误需要：H(C) ≥ 2d + 1

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 7: (b) Gibt es Übertr.Fehler die zuv. korrigiert werden?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt (b) Gibt es Übertr.Fehler die zuv. korrigiert werden?  


#### 题目中文翻译 / 中文题意

有可以可靠纠正的传输错误吗？
对于简单奇偶校验（H=2）：不能纠正任何错误，只能检测奇数位错误。
示例：
发送: 1100101|0
接收: 1100101|0
标记错误位。

#### 德文原题

```text
### (b) Gibt es Übertr.Fehler die zuv. korrigiert werden?

**有可以可靠纠正的传输错误吗？**

对于简单奇偶校验（H=2）：**不能纠正任何错误**，只能检测奇数位错误。

**示例：**

```
发送: 1100101|0
接收: 1100101|0
BCC:  00000001
```

标记错误位。

---

## VII. Routing & IPv4-Multicasting

## 路由与IPv4多播

---
```

#### 解答

**(b) Gibt es Übertr.Fehler die zuv. korrigiert werden?**

**有可以可靠纠正的传输错误吗？**

对于简单奇偶校验（H=2）：**不能纠正任何错误**，只能检测奇数位错误。

**示例：**

```
发送: 1100101|0
接收: 1100101|0
BCC:  00000001
```

标记错误位。

---

**VII. Routing & IPv4-Multicasting**

**路由与IPv4多播**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 8: Frage 32 / 第32题

**类型：** 考卷  

**来源说明：** Klausur 2017, Frage 32  


#### 题目中文翻译 / 中文题意

给定生成多项式 G = x³ + 1。
G在CRC中用多少位表示？

#### 德文原题

```text
### Frage 32 / 第32题

**Gegeben sei das Generatorpolynom G = x³ + 1.**  
**给定生成多项式 G = x³ + 1。**

**(a) Durch wie viele Bits wird G bei CRC repräsentiert?**  
**G在CRC中用多少位表示？**
```

#### 解答

**参考答案 / Lösung:** **4**

G = x³ + 1 = 1001（二进制），需要4位。

**(b) Es soll die Nachricht 11 00 11 CRC-geschützt übertragen werden. Berechnen Sie die zu übertragende Bitfolge (inkl. CRC-Prüfsumme!) unter der Verwendung des Generatorpolynoms G.**  
**消息11 00 11要进行CRC保护传输。使用生成多项式G计算要传输的位序列（包括CRC校验和）。**

**参考答案 / Lösung:**

**要传输的位序列：11 00 11 101**

**计算过程：**

1. 原始消息：110011
2. 添加(n-1)=3个零：110011000
3. 除以G=1001：

```
110011000 ÷ 1001 = ...
110011000
1001
----
 1011
 1001
 ----
  0101
  0000
  ----
   1010
   1001
   ----
    0110
    0000
    ----
     1100
     1001
     ----
      101  ← 余数（CRC校验和）
```

4. 传输：110011 + 101 = **110011101**

**(c) Nehmen Sie an, dass Sie die CRC-geschützte Bitfolge 10 01 10 01 empfangen haben. Zeigen Sie, dass die empfangene Bitfolge unter Verwendung des Generatorpolynoms G korrekt ist (inkl. Rechnung). Markieren Sie in Ihrer Rechnung die Stelle, an der der Empfänger die Korrektheit ablesen kann.**  
**假设您收到了CRC保护的位序列10 01 10 01。证明使用生成多项式G接收的位序列是正确的（包括计算）。在计算中标记接收方可以判断正确性的位置。**

**参考答案 / Lösung:**

```
10011001 ÷ 1001 = ...
10011001
1001
----
 0001
 0000
 ----
  0010
  0000
  ----
   0100
   0000
   ----
    1001
    1001
    ----
       0  ← 余数为0，表示正确！
```

**余数 = 0 → 校验和正确，传输无错误**

---

**总结 / Zusammenfassung**

本试卷涵盖了计算机网络的核心主题：

1. **OSI模型**：七层架构及各层功能
2. **协议基础**：TCP、IP、UDP、ICMP等
3. **寻址**：IPv4、IPv6、子网划分、CIDR
4. **DNS**：域名解析过程
5. **路由**：OSPF、路由表
6. **分片**：IP数据包分片机制
7. **TCP**：流量控制、拥塞控制、滑动窗口
8. **CRC**：循环冗余校验计算

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 9: Frage 4 / 第4题

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 4  


#### 题目中文翻译 / 中文题意

以下关于块校验字符(BCC)的哪些陈述是正确的？
这是一种错误检测方法。
它属于自纠错码家族。
它传输额外的奇偶校验位。
它能检测所有连续错误。

#### 德文原题

```text
### Frage 4 / 第4题

**Welche Aussagen treffen auf Block Check Character zu?**  
**以下关于块校验字符(BCC)的哪些陈述是正确的？**

- ☒ Es ist ein Verfahren zur Fehlererkennung.
    - 这是一种错误检测方法。
- ○ Es gehört zur Familie selbstkorrigierender Codes.
    - 它属于自纠错码家族。
- ☒ Es überträgt zusätzliche Paritätsbits.
    - 它传输额外的奇偶校验位。
- ○ Es erkennt alle zusammenhängenden Fehler.
    - 它能检测所有连续错误。
```

#### 解答

**解析：**

- ✓ 第一项正确：BCC是错误检测方法
- ✗ 第二项错误：BCC只能检测错误，不能纠正
- ✓ 第三项正确：BCC通过添加奇偶校验位实现
- ✗ 第四项错误：BCC不能检测所有连续错误（如偶数个位翻转可能检测不到）

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 10: Frage 5 / 第5题

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 5  


#### 题目中文翻译 / 中文题意

CRC方法使用生成多项式 G = x⁵ + x³ + x + 1。该多项式的正确位序列表示是什么？

#### 德文原题

```text
### Frage 5 / 第5题

**Beim CRC-Verfahren wird ein Generatorpolynom G = x⁵ + x³ + x + 1 verwendet. Welche ist die richtige Darstellung dieses Polynoms als Bitfolge?**  
**CRC方法使用生成多项式 G = x⁵ + x³ + x + 1。该多项式的正确位序列表示是什么？**

- ○ 110 101
- ○ 110 10
- ☒ 101 011
- ○ 110 011
```

#### 解答

**解析：**  
G = x⁵ + x³ + x + 1

将各项系数对应到位置：

- x⁵ → 位5 = 1
- x⁴ → 位4 = 0
- x³ → 位3 = 1
- x² → 位2 = 0
- x¹ → 位1 = 1
- x⁰ → 位0 = 1

从高位到低位：**101011**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 11: Frage 33 / 第33题

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 33  


#### 题目中文翻译 / 中文题意

给定生成多项式 G = x³ + x + 1。
G在CRC中用多少位表示？

#### 德文原题

```text
### Frage 33 / 第33题

**Gegeben sei das Generatorpolynom G = x³ + x + 1.**  
**给定生成多项式 G = x³ + x + 1。**

**(a) Durch wie viele Bits wird G bei CRC repräsentiert?**  
**G在CRC中用多少位表示？**
```

#### 解答

**Lösung / 答案：** **4**

G = x³ + x + 1 = 1011（二进制），需要4位。

**(b) 消息11 00 11要进行CRC保护传输。计算要传输的位序列。**

**Lösung / 答案：**

**计算过程：**

1. 原始消息：110011
2. 添加(n-1)=3个零：110011000
3. 除以G=1011：

```
110011000 ÷ 1011
110011000
1011
----
 01111
 1011
 ----
  1001
  1011
  ----
   0100
   0000
   ----
    1000
    1011
    ----
     011  ← 余数（CRC校验和）
```

**zu übertragende Bitfolge / 要传输的位序列：** **110011011**

**(c) 验证接收到的位序列10 01 11 00是否正确**

**Lösung / 答案：**

```
10011100 ÷ 1011
10011100
1011
----
 0101
 0000
 ----
  1011
  1011
  ----
   0001
   0000
   ----
    0010
    0000
    ----
     0100
     0000
     ----
      100  ← 余数不为0！
```

**结果：余数 = 100 ≠ 0**

**这说明接收的位序列包含错误！** 如果余数为0才表示正确。

---

**总结 / Zusammenfassung**

本试卷涵盖了计算机网络的核心主题：

| 章节   | 主题                     | 分值  |
| ---- | ---------------------- | --- |
| I    | 选择题（OSI模型、协议、NAT、TCP等） | 10分 |
| II   | OSI层模型（PDU/SDU、接口划分）   | 4分  |
| III  | DNS（查询过程、传输协议）         | 6分  |
| IV   | 协议协作（MAC寻址、路由）         | 6分  |
| V    | IP分片                   | 5分  |
| VI   | 寻址（IPv4子网、IPv6）        | 10分 |
| VII  | TCP（拥塞控制、流量控制）         | 10分 |
| VIII | 以太网、CSMA/CD            | 9分  |
| IX   | 信号调制、采样                | 5分  |
| X    | CRC校验                  | 7分  |

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 12: Frage 25 / 第25题

**类型：** 考卷  

**来源说明：** Klausur 2019, Frage 25  


#### 题目中文翻译 / 中文题意

(a) 将消息分成两个16-bit段并计算Internet校验和 (3分)

#### 德文原题

```text
### Frage 25 / 第25题

**(a) 将消息分成两个16-bit段并计算Internet校验和 (3分)**
```

#### 解答

**Lösung / 答案：**

**步骤：**

1. 分成两个16位段（每个字符8位，需要补前导零）：
    
    - 段1：R + N = 01010010 01001110
    - 段2：V + S = 01010110 01000011
    
2. 二进制加法：
    

```
  01010010 01001110  (RN)
+ 01010110 01000011  (VS)
-------------------
  10101000 10010001
```

3. 取反得到校验和：

```
  10101000 10010001 → 取反 → 01010111 01101110
```

**Checksumme / 校验和：** **0101 0111 0110 1110** 或十六进制 **576E**

**(b) Welche Schritte muss der Empfänger ausführen, um die Nachricht als korrekt zu verifizieren? (2分)**  
**接收方需要执行哪些步骤来验证消息是否正确？**

**Lösung / 答案：**

1. 将接收到的数据（包括校验和）分成16位段
2. 对所有段（包括校验和）进行二进制加法
3. 处理进位（回卷加到最低位）
4. 对结果取反
5. 如果结果全为0，则消息正确；否则有错误

或者简单地说：将所有16位字（包括校验和）相加，结果应该是全1（0xFFFF）。

---

**总结 / Zusammenfassung**

| 章节     | 主题                | 分值      |
| ------ | ----------------- | ------- |
| I      | 一般知识（延迟、强度、子网、路由） | 12分     |
| II     | Wireshark分析       | 8分      |
| III    | DNS               | 6分      |
| IV     | 协议协作              | 7分      |
| V      | IP分片              | 7分      |
| VI     | IP和路由             | 9分      |
| VII    | TCP（拥塞控制、流量控制）    | 17分     |
| VIII   | CSMA/CD和以太网       | 9分      |
| IX     | UDP校验和            | 5分      |
| **总计** |                   | **80分** |

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---
