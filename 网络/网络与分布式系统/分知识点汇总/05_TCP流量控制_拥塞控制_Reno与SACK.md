# TCP 流量控制、拥塞控制、Reno 与 SACK

## 知识点总结

- 流量控制保护接收方，拥塞控制保护网络。
- Slow Start 指数增长，Congestion Avoidance 线性增长。
- Reno 用 Fast Retransmit/Fast Recovery 处理三重复 ACK。

## 完整题目与解答汇总

### 题目 1: 2. Selektive Quittungen (H)

**类型：** 作业  

**来源说明：** Uebungsblatt 06, Aufgabe 2  


#### 题目中文翻译 / 中文题意

在 TCP 中使用选择确认 SACK：给定 8 个段和丢失模式，写出累计 ACK 与 SACK block。

#### 德文原题

```text
2. Selektive Quittungen (H)
Der Verlust einzelner Segmente kann zur unnötigen Wiederholung großer Datenmengen führen, insbe-
sondere bei Pfaden mit hoher Netzverzögerung. Durch die Einführung selektiver Quittungen (SACK,
siehe RFC 2018) kann dieses Problem gemildert werden.
1 501 1001 1501 2001
ACK SACK
Statt wie bei den „normalen” kumulativen TCP-Quittungen den bis dahin korrekt empfangenen zusam-
menhängenden Byte-Strom zu quittieren, kann ein Empfänger mit selektiven Quittungen zusätzliche
Segmente (oder zusammenhängende Folgen von Segmenten, sogenannte Blöcke) im Options-Feld des
TCP-Headers als empfangen notieren. Hierzu wird ein Bereich des Byte-Stroms mit Anfangs- und End-
Byte notiert. Auf Grundlage der obigen Abbildung würde etwa eine Quittung mit AckNr=501 und
SACK-Block=(1001,2001) gesendet; das verlorene Segment mit Bytes 501–1000 (schraffiert) wird so
ausgespart. Es können mehrere solche Blöcke in den Optionen angegeben werden.
Gehen Sie von einem Sender mit aktueller SeqNr=5000 aus. Der Sender sendet 8 Segmente von jeweils
500 Byte Länge. Wie werden kumulative und selektive Quittungen benutzt, wenn:
(a) die ersten vier Segmente empfangen werden, die letzten vier aber verloren gehen?
(b) das zweite, vierte und sechste und achte Segment verloren gehen?
Geben Sie ihre Lösung in Form einer Tabelle oder als Sequenzdiagramm an.
```

#### 解答

**2. Selektive Quittungen / SACK**

**DE:** Sender startet mit SeqNr 5000 und sendet 8 Segmente a 500 Byte.

**中文：** 发送方当前序列号为 5000，连续发送 8 个段，每段 500 B。因此每个段对应的字节范围如下：

| Segment | Bytebereich |
|---:|---|
| 1 | 5000-5499 |
| 2 | 5500-5999 |
| 3 | 6000-6499 |
| 4 | 6500-6999 |
| 5 | 7000-7499 |
| 6 | 7500-7999 |
| 7 | 8000-8499 |
| 8 | 8500-8999 |

**(a)**  
**DE:** Erste vier Segmente werden empfangen, die letzten vier gehen verloren.

**中文：** 前四个段都收到，后四个段丢失。因此连续字节流已经收到到 6999，下一个期望字节是 7000：

```text
Kumulatives ACK = 7000
SACK = leer, weil kein Block nach der Luecke empfangen wurde
```

**(b)**  
**DE:** Die Segmente 2, 4, 6 und 8 gehen verloren.

**中文：** 第 2、4、6、8 个段丢失。第 1 个段收到后，累计 ACK 只能到 5500；后面收到的第 3、5、7 段是不连续块，需要用 SACK 标出：

Nach Segment 1: ACK 5500. Danach kommen 3,5,7 ausserhalb der Luecken an:

```text
ACK = 5500
SACK-Bloecke = (6000,6500), (7000,7500), (8000,8500)
```

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 2: 3. TCP Reno (H)

**类型：** 作业  

**来源说明：** Uebungsblatt 06, Aufgabe 3  


#### 题目中文翻译 / 中文题意

根据 TCP Reno 的拥塞窗口图，识别 Slow Start、Congestion Avoidance、丢包轮次、threshold、Fast Recovery 和 Tahoe/Reno 差异。

#### 德文原题

```text
3. TCP Reno (H)
45
40
35
30
25
20
15
10
5
0
1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
Runde
)SSM(
niWgnoC
Abbildung 2: Verhalten von TCP Reno
In Abbildung 2 ist das Verhalten von TCP-Reno nach RFC 2581 zu sehen. TCP-Reno verhält sich ähnlich
wie TCP-Tahoe, jedoch mit Unterstützung für Fast Recovery:
Empfängt der Sender 3 ACK-Duplikate, geht er – anstelle von Slow Start – in den Fast Recovery Zustand
über. In dieser Phase wird das verlorene Paket noch vor dem Timeout wiederholt. Anschließend wird
zwar auch Threshold auf CongWin gesetzt, allerdings startet auch das neue CongWin mit diesem Wert.
2
TCP-Reno fährt also direkt mit der linearen Phase (Congestion Avoidance) fort.
In dieser Aufgabe nehmen wir an, dass während Fast Recovery keine weiteren Duplikate auftreten. Das
heißt, dass das erneut übertragene Segment — die Ursache für die ACK-Duplikate — erfolgreich quittiert
wird.
(a) Identifizieren Sie die Intervalle, in denen TCP Slow Start aktiv ist.
(b) Identifizieren Sie die Intervalle, in denen TCP Congestion Avoidance aktiv ist.
(c) In welchen Runden trat ein Paketverlust auf? Wurde dieser durch duplizierte ACKs oder durch die
Überschreitung des Timeouts erkannt?
(d) Welchen Wert hat Threshold zu Beginn (in der 1. Runde), in der 18. und in der 24. Runde?
(e) In welcher Runde wird das 70. Segment gesendet?
(f) Angenommen in der 26. Runde wird ein Paketverlust durch ein (dreifaches) ACK-Duplikat festge-
stellt. Wie werden die Werte von CongWin und Threshold anschließend sein?
(g) Angenommen es würde TCP Tahoe statt Reno genutzt, wie würde sich das Verhalten nach dem
ersten Paketverlust ändern? Welche Werte haben Threshold und CongWin in der 19. Runde?
```

#### 解答

**3. TCP Reno**

![Blatt 06 Seite 2: TCP Reno Congestion Window](pictures/blatt-06_pages-1-2-2.png)

**DE:** Aus der Grafik liest man ungefaehr folgende CongWin-Werte ab:

```text
Runde:   1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
CongWin: 1  2  4  8 16 32 33 34 35 36 37 38 39 40 41 42 21 22 23 24 25 26  1  2  4  8
```

**中文：** 按图中的点读数，答案如下：

| Teil | Antwort |
|---|---|
| (a) Slow Start | Runde `1-6`，以及 timeout 后的 Runde `23-26` |
| (b) Congestion Avoidance | Runde `6-16`，以及 Fast Recovery 后的 Runde `17-22` |
| (c) Paketverlust | Runde `16` 后窗口从 `42` 降到 `21`，是 `3 Duplicate ACKs / Fast Recovery`；Runde `22` 后窗口从 `26` 降到 `1`，是 `Timeout` |
| (d) Threshold | Runde 1: `32`；Runde 18: `21`；Runde 24: `13` |
| (e) 70. Segment | 到 Runde 6 共发送 `1+2+4+8+16+32 = 63` 个段，所以第 70 个段在 Runde `7` |
| (f) Verlust in Runde 26 durch 3 DupACK | Runde 26 时 `CongWin=8`，所以 `Threshold=4`，Fast Recovery 结束后的新 `CongWin=4` |
| (g) TCP Tahoe statt Reno | 第一次丢包后 Tahoe 不做 Fast Recovery，而是 `Threshold=21`、`CongWin=1` 并重新 Slow Start；因此 Runde 19 时 `Threshold=21`、`CongWin=4` |

**Wissen / 知识点：** Reno 区分 Timeout 和 Triple-Duplicate-ACK。Timeout 更严重，`CongWin` 回到 1；Triple-Duplicate-ACK 触发 Fast Retransmit/Fast Recovery，阈值减半后继续线性增长。

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 3: 4. Approximierung von TCP Durchsatz und Verlustrate

**类型：** 作业  

**来源说明：** Uebungsblatt 06, Aufgabe 4  


#### 题目中文翻译 / 中文题意

推导 TCP 拥塞避免阶段的丢包率公式，并由丢包率近似 TCP 平均吞吐率。

#### 德文原题

```text
4. Approximierung von TCP Durchsatz und Verlustrate
Aus den Vorlesungsfolien geht hervor, dass der TCP Throughput während einer bestehenden Verbindung
immer im Intervall zwischen W bis W pendelt. Während jeder Congestion-Avoidance-Phase geht
2∗RT D RT D
1 Paket verloren.
(a) Zeigen Sie, dass die Verlustrate L (Anteil verloren gegangener Segmente) dem folgenden Term
entspricht:
1
L =
3 W 2 + 3 W
8 4
(b) Zeigen Sie, dass mit gegebener Verlustrate L während einer TCP-Verbindung die durchschnittliche
Übertragungsrate1 mit dem folgenden Term approximiert werden kann:
1.22 · MSS
R ≈ √
RTD · L
1erfolgreich übertragener Segmente
```

#### 解答

**4. TCP-Durchsatz und Verlustrate**

**DE:** In Congestion Avoidance pendelt das Fenster zwischen `W/2` und `W`. Pro Zyklus werden ungefaehr:

**中文：** 在拥塞避免阶段，窗口在 `W/2` 到 `W` 之间线性增长。一个周期中发送的段数约为：

```text
(W/2 + W) / 2 * (W/2) = 3W^2/8
```

Segmente gesendet; genauer mit diskreter Summation kommt:

```text
3W^2/8 + 3W/4
```

**DE:** Da pro Zyklus ein Segment verloren geht:

**中文：** 因为假设每个周期丢 1 个段，所以丢包率等于“1 除以一个周期内发送的段数”：

```text
L = 1 / (3W^2/8 + 3W/4)
```

**DE:** Nach Umstellen gilt naeherungsweise:

**中文：** 将丢包率公式近似变形，可以得到窗口大小和吞吐率的近似关系：

```text
W approx sqrt(8/(3L))
R approx (0.75 * W * MSS) / RTD
R approx 1.22 * MSS / (RTD * sqrt(L))
```

**Wissen / 知识点：** TCP 吞吐量与丢包率的平方根成反比：丢包率稍微上升，吞吐量会明显下降。

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 4: 4. Staukontrolle bei TCP (H)

**类型：** 作业  

**来源说明：** Uebungsblatt 07, Aufgabe 4  


#### 题目中文翻译 / 中文题意

在 TCP Slow Start 中给定 RTT、数据量、段大小和速率，画传输时序并比较有无 Slow Start 的总传输时间。

#### 德文原题

```text
4. Staukontrolle bei TCP (H)
Bei TCP kommt Slow-Start als Mechanismus zur Staukontrolle zum Einsatz. Es sollen über eine TCP-
Verbindung mit Netzverzögerung von 100 ms (d.h. RTD = 200 ms) 7500 B Nutzdaten in 15 Segmenten
gleicher Größe S = 500 B vom Server zum Client übertragen werden.
Hinweis: Nehmen Sie einen nicht erreichbaren Wert für Threshold an, und vernachlässigen Sie Ver-
luste, Empfangspuffergröße, Verarbeitungsverzögerung, die Übertragungszeit für Protokollheader sowie
den Verbindungsabbau. Gehen Sie also davon aus, dass die beiden TCP-Instanzen sofort nach dem
Verbindungsaufbau die Übertragung beginnen und danach die Slow-Start-Phase nicht verlassen.
(a) Sei die Übertragungsrate R = 20 kB/s. Erstellen Sie ein Sequenzdiagramm für die Übertragung und
tragen Sie die Größe des jeweils aktuellen Überlastfensters (CongWin) in das Diagramm ein!
(b) Bestimmen Sie die Übertragungsdauer gemessen vom Absenden des SYN des Clients bis alle Nutz-
daten empfangen wurden:
i. mit Slow-Start.
ii. ohne Slow-Start mit fester Fenstergröße von 20.
(c) Wie lange würde die Übertragung jeweils mit und ohne Slow-Start für R = 500 kB/s dauern?
```

#### 解答

**4. Staukontrolle bei TCP / TCP Slow Start**

Gegeben: RTD `200 ms`, Einweg `100 ms`, 15 Segmente a `500 B`, Threshold unerreichbar.

**(a)** Slow Start sendet pro RTT: `1, 2, 4, 8` Segmente, zusammen `15`. CongWin entwickelt sich also `1 -> 2 -> 4 -> 8 -> 16`.

**中文：** Slow Start 每个 RTT 让拥塞窗口大约翻倍，所以四轮可以分别发送 `1, 2, 4, 8` 个段，正好覆盖 15 个段。拥塞窗口变化为 `1 -> 2 -> 4 -> 8 -> 16`。

**(b)**  
**DE:** Bei `R = 20 kB/s` dauert ein Segment `500/20000 = 25 ms`. Einschliesslich 3-Way-Handshake bis Server senden kann grob `300 ms`, danach vier Slow-Start-Runden. Die letzte Runde enthaelt 8 Segmente und braucht wegen Serialisierung `8*25 ms = 200 ms` plus Wegzeit. Ohne Slow Start mit Fenster 20 koennen alle 15 Segmente sofort in einer Sendefolge gesendet werden: `15*25 ms + 100 ms` nach Verbindungsaufbau.

**中文：** 当速率 `R = 20 kB/s` 时，一个 500 B 的段发送时间是 `500 / 20000 = 25 ms`。从客户端 SYN 开始算，三次握手到服务器可以开始发送，大约需要 `300 ms`。之后 Slow Start 需要四轮发送：第 1 轮发 1 个，第 2 轮发 2 个，第 3 轮发 4 个，第 4 轮发 8 个。最后一轮因为 8 个段要排队串行发出，所以仅发送就要 `8 * 25 ms = 200 ms`，再加上最后一个段到客户端的传播时间。若没有 Slow Start，固定窗口为 20，则 15 个段可以连续发出，发送时间为 `15 * 25 ms = 375 ms`，再加上最后一个段传播到客户端的 `100 ms`。

Kurzform:

| Fall | Dauer ab SYN, grob |
|---|---:|
| mit Slow Start | `300 ms + 3*200 ms + 200 ms + 100 ms = 1200 ms` |
| Fenster 20 | `300 ms + 375 ms + 100 ms = 775 ms` |

**(c)** Bei `R = 500 kB/s` dauert ein Segment `1 ms`; dann dominiert RTT:

| Fall | Dauer ab SYN, grob |
|---|---:|
| mit Slow Start | `300 ms + 3*200 ms + 8 ms + 100 ms = 1008 ms` |
| Fenster 20 | `300 ms + 15 ms + 100 ms = 415 ms` |

**中文：** 当 `R = 500 kB/s` 时，一个段只需要 `500 / 500000 = 1 ms` 发送，发送时间已经很小，主要耗时变成 RTT 和握手等待。Slow Start 仍然要等多轮 RTT，所以约 `1008 ms`；固定窗口 20 可以一次连续发送完 15 个段，所以约 `415 ms`。

**Wissen / 知识点：** 带宽大、RTT 大时，窗口机制非常关键；Slow Start 的指数增长避免一开始压垮网络，但会增加短连接时延。

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 5: Flusssteuerung / 流量控制

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Flusssteuerung / 流量控制  


#### 题目中文翻译 / 中文题意

定义： 防止发送方发送速度超过接收方处理速度的机制。
TCP使用滑动窗口实现流量控制，接收方通过Window字段告知发送方可用缓冲区大小。

#### 德文原题

```text
### Flusssteuerung / 流量控制

**定义：** 防止发送方发送速度超过接收方处理速度的机制。

TCP使用**滑动窗口**实现流量控制，接收方通过Window字段告知发送方可用缓冲区大小。

---
```

#### 解答

**Flusssteuerung / 流量控制**

**定义：** 防止发送方发送速度超过接收方处理速度的机制。

TCP使用**滑动窗口**实现流量控制，接收方通过Window字段告知发送方可用缓冲区大小。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 6: Staukontrolle / 拥塞控制

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Staukontrolle / 拥塞控制  


#### 题目中文翻译 / 中文题意

定义： 防止网络过载的机制。
TCP拥塞控制算法：Slow Start, Congestion Avoidance, Fast Retransmit, Fast Recovery

#### 德文原题

```text
### Staukontrolle / 拥塞控制

**定义：** 防止网络过载的机制。

TCP拥塞控制算法：Slow Start, Congestion Avoidance, Fast Retransmit, Fast Recovery

---
```

#### 解答

**Staukontrolle / 拥塞控制**

**定义：** 防止网络过载的机制。

TCP拥塞控制算法：Slow Start, Congestion Avoidance, Fast Retransmit, Fast Recovery

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 7: (b) Phase a?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt (b) Phase a?  


#### 题目中文翻译 / 中文题意

阶段a是什么？

#### 德文原题

```text
### (b) Phase a?

**阶段a是什么？**

**Lösung / 答案：** **Slow Start / 慢启动**

特点：窗口指数增长

---
```

#### 解答

**Lösung / 答案：** **Slow Start / 慢启动**

特点：窗口指数增长

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 8: (c) Phase b?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt (c) Phase b?  


#### 题目中文翻译 / 中文题意

阶段b是什么？

#### 德文原题

```text
### (c) Phase b?

**阶段b是什么？**

**Lösung / 答案：** **Congestion Avoidance / 拥塞避免**

特点：窗口线性增长（每RTT增加1个MSS）

---
```

#### 解答

**Lösung / 答案：** **Congestion Avoidance / 拥塞避免**

特点：窗口线性增长（每RTT增加1个MSS）

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 9: (d) Vervollst. Phase c

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt (d) Vervollst. Phase c  


#### 题目中文翻译 / 中文题意

完成阶段c的图
阶段c是timeout后的恢复：
cwnd重置为1
ssthresh设为原cwnd的一半
重新进入Slow Start

#### 德文原题

```text
### (d) Vervollst. Phase c

**完成阶段c的图**

阶段c是**timeout后的恢复**：

- cwnd重置为1
- ssthresh设为原cwnd的一半
- 重新进入Slow Start

---
```

#### 解答

**(d) Vervollst. Phase c**

**完成阶段c的图**

阶段c是**timeout后的恢复**：

- cwnd重置为1
- ssthresh设为原cwnd的一半
- 重新进入Slow Start

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 10: (f) Verfahren zur Optimierung d. Überlastkontrollalgorithmus nach timeout?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt (f) Verfahren zur Optimierung d. Überlastkontrollalgorithmus nach timeout?  


#### 题目中文翻译 / 中文题意

超时后优化拥塞控制算法的方法？

#### 德文原题

```text
### (f) Verfahren zur Optimierung d. Überlastkontrollalgorithmus nach timeout?

**超时后优化拥塞控制算法的方法？**
```

#### 解答

**Lösung / 答案：** 如果题目严格说的是 **timeout 已经发生之后**，则应当：

- `ssthresh = cwnd/2`
- `cwnd = 1`
- 重新进入 `Slow Start`

**补充：** `Fast Retransmit / Fast Recovery` 是为了在出现 `3 Duplicate ACKs` 时尽早恢复，避免等到 timeout；它不是 timeout 已经发生后的处理方式。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 11: Flusssteuerung / 流量控制

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Flusssteuerung / 流量控制  


#### 题目中文翻译 / 中文题意

场景： 发送方发送数据，接收方有4096B缓冲区

#### 德文原题

```text
## Flusssteuerung / 流量控制

**场景：** 发送方发送数据，接收方有4096B缓冲区

---
```

#### 解答

**Flusssteuerung / 流量控制**

**场景：** 发送方发送数据，接收方有4096B缓冲区

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 12: Frage 8 / 第8题

**类型：** 考卷  

**来源说明：** Klausur 2017, Frage 8  


#### 题目中文翻译 / 中文题意

以下哪些参数与TCP慢启动相关？
传输阈值的高度
接收窗口的大小
拥塞避免缓冲区的大小
拥塞窗口的大小

#### 德文原题

```text
### Frage 8 / 第8题

**Welche/r der folgenden Parameter sind/ist relevant bei TCP Slow Start?**  
**以下哪些参数与TCP慢启动相关？**

- ☒ Höhe der Übertragungsschwelle (threshold)
    - 传输阈值的高度
- ☒ Größe des Empfangsfensters (receive window)
    - 接收窗口的大小
- ○ Größe des Stauvermeidungspuffers (congestion avoidance buffer)
    - 拥塞避免缓冲区的大小
- ☒ Größe des Überlastungsfensters (congestion window)
    - 拥塞窗口的大小
```

#### 解答

**参考答案 / Lösung:**

- ✓ Threshold：决定何时从慢启动切换到拥塞避免
- ✓ Receive Window：限制发送方的发送量
- ✗ Congestion avoidance buffer：不是标准TCP参数
- ✓ Congestion Window：慢启动期间指数增长的核心参数

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 13: Frage 9 / 第9题

**类型：** 考卷  

**来源说明：** Klausur 2017, Frage 9  


#### 题目中文翻译 / 中文题意

以下关于TCP流量控制的哪些陈述是正确的？
流量控制减轻发送方的负担。
流量控制减轻接收方的负担。
流量控制限制传输网络中的连接数量。
TCP-Tahoe和TCP-Reno是流量控制方法。
TCP使用滑动窗口协议进行流量控制。

#### 德文原题

```text
### Frage 9 / 第9题

**Welche Aussagen über TCP-Flusssteuerung treffen zu?**  
**以下关于TCP流量控制的哪些陈述是正确的？**

- ○ Die Flusssteuerung entlastet den Sender.
    - 流量控制减轻发送方的负担。
- ☒ Die Flusssteuerung entlastet den Empfänger.
    - 流量控制减轻接收方的负担。
- ○ Die Flusssteuerung begrenzt die Anzahl der Verbindungen in Transitnetzen.
    - 流量控制限制传输网络中的连接数量。
- ○ TCP-Tahoe und TCP-Reno sind Verfahren der Flusssteuerung.
    - TCP-Tahoe和TCP-Reno是流量控制方法。
- ☒ TCP verwendet das Sliding-Window-Protokoll zur Flusssteuerung.
    - TCP使用滑动窗口协议进行流量控制。
```

#### 解答

**参考答案 / Lösung:**

- ✗ 第一项错误：流量控制保护接收方，不是发送方
- ✓ 第二项正确：防止接收方缓冲区溢出
- ✗ 第三项错误：这是拥塞控制的功能，不是流量控制
- ✗ 第四项错误：Tahoe和Reno是拥塞控制算法
- ✓ 第五项正确：滑动窗口是流量控制的核心机制

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 14: Frage 25: Überlastkontrolle / 拥塞控制

**类型：** 考卷  

**来源说明：** Klausur 2017, Frage 25  


#### 题目中文翻译 / 中文题意

下图显示了TCP-Tahoe拥塞控制中一个关键参数的时间演变。过程分为a、b、c三个阶段。
y轴的标签是什么？

#### 德文原题

```text
### Frage 25: Überlastkontrolle / 拥塞控制

**Die folgende Abbildung zeigt die Zeitliche Entwicklung einer kritischen Kenngröße bei der TCP-Tahoe Überlastkontrolle. Der Ablauf ist in die Phasen a, b und c eingeteilt.**  
**下图显示了TCP-Tahoe拥塞控制中一个关键参数的时间演变。过程分为a、b、c三个阶段。**

**(a) Wie lautet die Beschriftung für die y-Achse?**  
**y轴的标签是什么？**
```

#### 解答

**参考答案 / Lösung:** **Congestion Window (cwnd) / 拥塞窗口**

**(b) Benennen Sie Phase a:**  
**命名阶段a：**

**参考答案 / Lösung:** **Exponentielle Phase (Slow Start) / 指数阶段（慢启动）**

窗口指数增长，每个RTT翻倍。

**(c) Benennen Sie Phase b:**  
**命名阶段b：**

**参考答案 / Lösung:** **Lineare Phase (Congestion Avoidance) / 线性阶段（拥塞避免）**

达到阈值后，窗口线性增长。

**(d) Was ist ein frühzeitiger Indikator für Paketverluste, noch bevor ein timeout für ein Paket auftritt?**  
**在数据包超时之前，什么是数据包丢失的早期指标？**

**参考答案 / Lösung:** **Duplicate ACKs / 重复确认**

收到3个重复ACK表示可能发生了丢包。

**(e) Nennen Sie ein Verfahren zur Optimierung des Überlastkontrollalgorithmus nach Auftreten eines timeouts.**  
**列举一种在超时发生后优化拥塞控制算法的方法。**

**参考答案 / Lösung:** **Threshold halbieren und neuer Slow Start / 阈值减半并重新慢启动**

TCP Tahoe：超时后cwnd=1，ssthresh=cwnd/2，重新开始慢启动。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 15: Frage 26: Flusssteuerung / 流量控制

**类型：** 考卷  

**来源说明：** Klausur 2017, Frage 26  


#### 题目中文翻译 / 中文题意

下图显示了使用TCP的数据传输（忽略连接建立和拆除）。
用序列号(Seq, 1x)、确认号(ACK, 2x)和窗口大小(Win, 2x)补充图中的五个缺失值。

#### 德文原题

```text
### Frage 26: Flusssteuerung / 流量控制

**Die Folgende Abbildung zeigt eine Datenübertragung unter Verwendung von TCP unter Vernachlässigung des Verbindungsauf- und -abbaus.**  
**下图显示了使用TCP的数据传输（忽略连接建立和拆除）。**

**(a) Vervollständigen Sie die fünf fehlenden Angaben in der Abbildung mit Sequenznummern (Seq, 1x), Beätätigungsnummern (ACK, 2x) und Fenstergrößen (Win, 2x).**  
**用序列号(Seq, 1x)、确认号(ACK, 2x)和窗口大小(Win, 2x)补充图中的五个缺失值。**
```

#### 解答

**参考答案 / Lösung:**

根据图示：

- 第一个数据段：1024 Bytes, Seq = 0
- 第一个ACK：ACK = **1024**, Win = **3072**
- 第二个数据段：2048 Bytes, Seq = **1024**
- 第二个ACK：ACK = **3072**, Win = **1024**

**解释：**

- 发送1024字节后，接收方期望下一个字节是1024，所以ACK=1024
- 接收方缓冲区4096字节，已用1024，剩余3072，所以Win=3072
- 第二个数据段从1024开始
- 收到2048字节后，共收到3072字节，ACK=3072，缓冲区剩余1024

**(b) Markieren Sie den belegten Pufferspeicher des Empfängers in oben stehender Abbildung.**  
**在上图中标记接收方已占用的缓冲区。**

**参考答案 / Lösung:**

- 初始状态：缓冲区4096字节，全空
- 收到第一个1024字节后：1个方格被标记（×）
- 收到第二个2048字节后：3个方格被标记（×××）

**(c) Durch welches Ereignis wird der Pufferspeicher auf Empfängerseite wieder frei?**  
**什么事件使接收方的缓冲区重新释放？**

**参考答案 / Lösung:** **Einlesen/Verarbeiten der Daten auf Empfängerseite / 接收方读取/处理数据**

当应用程序从TCP缓冲区读取数据时，缓冲区空间被释放。

---

**8 Cyclic Redundancy Check**

**8 循环冗余校验**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 16: Frage 9 / 第9题

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 9  


#### 题目中文翻译 / 中文题意

以下关于TCP流量控制的哪些陈述是正确的？
流量控制减轻发送方的负担。
流量控制减轻接收方的负担。
流量控制用于限制传输网络中的连接数量。
Tahoe算法是流量控制方法。
TCP使用滑动窗口协议进行流量控制。

#### 德文原题

```text
### Frage 9 / 第9题

**Welche Aussagen über TCP-Flusssteuerung treffen zu?**  
**以下关于TCP流量控制的哪些陈述是正确的？**

- ○ Die Flusssteuerung entlastet den Sender.
    - 流量控制减轻发送方的负担。
- ☒ Die Flusssteuerung entlastet den Empfänger.
    - 流量控制减轻接收方的负担。
- ○ Die Flusssteuerung dient zur Begrenzung der Anzahl von Verbindungen in Transitnetzen.
    - 流量控制用于限制传输网络中的连接数量。
- ○ Der Tahoe-Algorithmus ist ein Verfahren zur Flusssteuerung.
    - Tahoe算法是流量控制方法。
- ☒ TCP verwendet das Sliding-Window-Protokoll zur Flusssteuerung.
    - TCP使用滑动窗口协议进行流量控制。
```

#### 解答

**解析：**

- ✗ 第一项错误：流量控制保护接收方，不是发送方
- ✓ 第二项正确：防止接收方缓冲区溢出
- ✗ 第三项错误：这是拥塞控制的功能
- ✗ 第四项错误：Tahoe是拥塞控制算法，不是流量控制
- ✓ 第五项正确：滑动窗口是TCP流量控制的核心机制

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 17: Frage 26: Überlastkontrolle / 拥塞控制

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 26  


#### 题目中文翻译 / 中文题意

y轴的正确标签是什么？

#### 德文原题

```text
### Frage 26: Überlastkontrolle / 拥塞控制

**(a) Wie lautet die korrekte Beschriftung für die y-Achse?**  
**y轴的正确标签是什么？**
```

#### 解答

**Lösung / 答案：** **Congestion Window (cwnd) / 拥塞窗口** 或 **Segmente / 段数**

**(b) Benennen Sie Phase a:**  
**命名阶段a：**

**Lösung / 答案：** **Slow Start / 慢启动**（指数增长阶段）

**(c) Benennen Sie Phase b:**  
**命名阶段b：**

**Lösung / 答案：** **Congestion Avoidance / 拥塞避免**（线性增长阶段）

**(d) Was ist ein frühzeitiger Indikator für Segmentverlust?**  
**什么是数据段丢失的早期指标？**

**Lösung / 答案：** **Duplicate ACKs / 重复确认**（3个重复ACK）

**(e) Nennen Sie ein Verfahren zur Optimierung des Überlastkontrollalgorithmus nach Auftreten eines timeouts.**  
**列举一种超时后优化拥塞控制算法的方法。**

**Lösung / 答案：**

- **ssthresh = cwnd/2**（阈值减半）
- **cwnd = 1**（拥塞窗口重置为1）
- **重新开始慢启动**

**注意：** `Fast Recovery` 是 `3 Duplicate ACKs` 时的 Reno 机制，用来避免等到 timeout；如果 timeout 已经发生，处理就是回到 `cwnd=1` 并慢启动。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 18: Frage 27: Flusssteuerung / 流量控制

**类型：** 考卷  

**来源说明：** Klausur 2018, Frage 27  


#### 题目中文翻译 / 中文题意

(a) 补充图中的五个缺失值

#### 德文原题

```text
### Frage 27: Flusssteuerung / 流量控制

**(a) 补充图中的五个缺失值**
```

#### 解答

**Lösung / 答案：**

|参数|值|
|---|---|
|第一个ACK|ACK = **1024**|
|第一个Win|Win = **3072**|
|第二个Seq|Seq = **1024**|
|第二个ACK|ACK = **3072**|
|第二个Win|Win = **1024**|

**解释：**

- 发送1024字节（Seq=0）后，接收方期望下一个字节是1024，所以ACK=1024
- 缓冲区4096字节，已用1024，剩余Win=3072
- 第二个数据段从Seq=1024开始，发送2048字节
- 收到后，共收到3072字节，ACK=3072
- 缓冲区剩余4096-3072=1024，Win=1024

**(b) 标记接收方已占用的缓冲区**

- 第一次接收后：1/4格被占用（1024/4096）
- 第二次接收后：3/4格被占用（3072/4096）

**(c) 什么事件使接收方的缓冲区重新释放？**

**Lösung / 答案：** **应用程序从TCP缓冲区读取数据** / Anwendung liest Daten aus dem Puffer

---

**VIII. Ethernet, CSMA**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 19: Frage 19: Staukontrolle / 拥塞控制

**类型：** 考卷  

**来源说明：** Klausur 2019, Frage 19  


#### 题目中文翻译 / 中文题意

条件：
初始 CongWin = 1
第6轮和第9轮发生超时导致的段丢失
(a) 在图中画出11轮的窗口大小变化 (4分)

#### 德文原题

```text
### Frage 19: Staukontrolle / 拥塞控制

**条件：**

- 初始 CongWin = 1
- Threshold = 8
- 第6轮和第9轮发生超时导致的段丢失

**(a) 在图中画出11轮的窗口大小变化 (4分)**
```

#### 解答

**Lösung / 答案：**

| Runde | CongWin | Phase                | 说明                    |
| ----- | ------- | -------------------- | --------------------- |
| 1     | 1       | Slow Start           | 初始                    |
| 2     | 2       | Slow Start           | 指数增长                  |
| 3     | 4       | Slow Start           | 指数增长                  |
| 4     | 8       | Slow Start           | 达到Threshold           |
| 5     | 9       | Congestion Avoidance | 线性增长                  |
| 6     | 10      | **Timeout!**         | 丢包                    |
| 7     | 1       | Slow Start           | cwnd重置为1，新threshold=5 |
| 8     | 2       | Slow Start           | 指数增长                  |
| 9     | 4       | **Timeout!**         | 丢包                    |
| 10    | 1       | Slow Start           | cwnd=1，新threshold=2   |
| 11    | 2       | Slow Start/CA        | 达到threshold           |

**(b) In welcher Einheit wird die Y-Achse gemessen? (1分)**  
**Y轴的单位是什么？**

**Lösung / 答案：** **Segmente / 段** 或 **MSS (Maximum Segment Size)**

**(c) Benennen Sie die Phase von Runde 1 bis Runde 4. (1分)**  
**命名第1-4轮的阶段。**

**Lösung / 答案：** **Slow Start / 慢启动**

**(d) Benennen Sie die Phase von Runde 4 bis Runde 6. (1分)**  
**命名第4-6轮的阶段。**

**Lösung / 答案：** **Congestion Avoidance / 拥塞避免**

**(e) Welchen Wert hat Threshold in Runde 7? (1分)**  
**第7轮的Threshold值是多少？**

**Lösung / 答案：** **5**

第6轮超时时cwnd=10，新threshold = cwnd/2 = 10/2 = 5

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 20: Frage 20: TCP Reno (Fast Recovery)

**类型：** 考卷  

**来源说明：** Klausur 2019, Frage 20  


#### 题目中文翻译 / 中文题意

TCP Reno在什么事件时进入快速恢复阶段？

#### 德文原题

```text
### Frage 20: TCP Reno (Fast Recovery)

**(a) Bei welchem Ereignis wechselt TCP Reno in die Fast Recovery Phase? (1分)**  
**TCP Reno在什么事件时进入快速恢复阶段？**
```

#### 解答

**Lösung / 答案：** **3 Duplicate ACKs / 3个重复确认**

**(b) Welchen Wert hat CongWin, wenn in Runde 6 ein Segmentverlust durch 3 Quittungsduplikate festgestellt wird? (1分)**  
**如果在第6轮通过3个重复ACK检测到段丢失，CongWin的值是多少？**

**Lösung / 答案：** 进入 Fast Recovery 当下为 **8**；Fast Recovery 结束后回到 **5**。

TCP Reno在3个重复ACK时：

- ssthresh = cwnd/2 = 10/2 = 5
- cwnd = ssthresh + 3 = 8（进入 Fast Recovery 时）
- 当重传段被新的 ACK 确认后，cwnd 设置为 ssthresh，即 5

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 21: Frage 21: Flusssteuerung / 流量控制

**类型：** 考卷  

**来源说明：** Klausur 2019, Frage 21  


#### 题目中文翻译 / 中文题意

场景：TCP连接结束阶段，C0和S0交叉
什么决定了TCP序列号的增长？

#### 德文原题

```text
### Frage 21: Flusssteuerung / 流量控制

**场景：TCP连接结束阶段，C0和S0交叉**

**(a) Was ist ausschlaggebend für das Fortzählen von TCP-Sequenznummern? (1分)**  
**什么决定了TCP序列号的增长？**
```

#### 解答

**Lösung / 答案：** **发送的字节数 / Anzahl der gesendeten Bytes**

每发送一个字节，序列号增加1。SYN和FIN也各占一个序列号。

**(b) Welche gesetzten Flags im Header eines TCP-Segments werden immer durch ein ACK quittiert? (1分)**  
**TCP段头部中哪些设置的标志总是需要ACK确认？**

**Lösung / 答案：** **SYN, FIN**

SYN和FIN标志都需要确认，并且各消耗一个序列号。

**(c) Welche Flags im TCP-Header sind in dem Segment C1 gesetzt? (1分)**  
**段C1中设置了哪些TCP标志？**

**Lösung / 答案：** **ACK**

C1是对S1的确认。

**(d) Der Client begann den Verbindungsaufbau. Wieviele Segmente ohne Nutzdaten hat der Client im Verlauf der gesamten Kommunikation mindestens an den Server geschickt? (1分)**  
**客户端发起连接。客户端在整个通信过程中至少发送了多少个不含数据的段给服务器？**

**Lösung / 答案：** **3**

- 三次握手的SYN（1个）
- 三次握手的ACK（1个）
- 四次挥手的最后ACK（1个）

**(e) 补充表格中C0和S1的序列号和确认号 (3分)**

**Lösung / 答案：**

已知：S0 发送1000字节，Seq=1051, Ack=6000

|Segment|Seq-nummer|Ack-Nummer|
|---|---|---|
|S0|1051|6000|
|C0|**6000**|**2051**|
|S1|**2051**|**6001**|

**解释：**

- S0：服务器发送1000字节数据（Seq=1051），期望客户端下一个字节是6000（Ack=6000）
- C0：客户端FIN段，Seq=6000（接着上次），Ack=1051+1000=2051（确认收到S0的数据）
- S1：服务器FIN段，Seq=2051，Ack=6000+1=6001（FIN占一个序列号）

---

**VIII. CSMA und Ethernet (9分)**

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---
