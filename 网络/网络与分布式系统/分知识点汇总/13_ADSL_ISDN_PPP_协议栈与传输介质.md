# ADSL、ISDN、PPP、协议栈与传输介质

## 知识点总结

- ADSL/ISDN/PPP/ATM 题常考协议栈封装和单位换算。
- PPP 的 LCP 管链路，NCP 管网络层协议配置。

## 完整题目与解答汇总

### 题目 1: 第1题

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt 第1题  


#### 题目中文翻译 / 中文题意

ISO/OSI哪一层涉及介质的物理特性？
选项： Bitübertragung（比特传输层）, Anwendung（应用层）, Kommunikation（通信层）, Vermittlung（网络层）

#### 德文原题

```text
### 第1题

**Welche Schichten d. ISO/OSI betr. physikal. Eigensch. von Medien (Bitübertr., Anw., Komm., Vermittl.)?**  
**ISO/OSI哪一层涉及介质的物理特性？**

**选项：** Bitübertragung（比特传输层）, Anwendung（应用层）, Kommunikation（通信层）, Vermittlung（网络层）
```

#### 解答

**Lösung / 答案：** **Bitübertragungsschicht / 物理层 (Schicht 1)**

物理层定义电气、机械、功能和过程特性，包括电压、频率、电缆类型等。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 2: Störeffekte, welche bei elektr. aber nicht Lichtwellenleiter?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Störeffekte, welche bei elektr. aber nicht Lichtwellenleiter?  


#### 题目中文翻译 / 中文题意

电导体有但光纤没有的干扰效应？

#### 德文原题

```text
### Störeffekte, welche bei elektr. aber nicht Lichtwellenleiter?

**电导体有但光纤没有的干扰效应？**
```

#### 解答

**Lösung / 答案：**

- **Elektromagnetische Interferenz (EMI) / 电磁干扰**
- **Übersprechen (Crosstalk) / 串扰**
- **Induktion / 感应**

光纤使用光信号，不受电磁干扰影响。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 3: Lichtwellenleiter: 2 Klassen (Kerndurchmesser) nennen

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Lichtwellenleiter: 2 Klassen (Kerndurchmesser) nennen  


#### 题目中文翻译 / 中文题意

光纤：列举两类（按纤芯直径）

#### 德文原题

```text
### Lichtwellenleiter: 2 Klassen (Kerndurchmesser) nennen

**光纤：列举两类（按纤芯直径）**
```

#### 解答

**Lösung / 答案：**

|类型|Kerndurchmesser|特性|
|---|---|---|
|**Multimode / 多模**|50-62.5 μm|短距离，较便宜|
|**Singlemode / 单模**|8-10 μm|长距离，高带宽|

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 4: ADSL频段分配

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt ADSL频段分配  


#### 题目中文翻译 / 中文题意

|频段|范围(kHz)|用途|
|0-4|POTS（普通电话）||
|4-138|上行数据||
|138-276|保护带||
|276-1041|下行数据||

#### 德文原题

```text
### ADSL频段分配

|频段|范围(kHz)|用途|
|---|---|---|
|0-4|POTS（普通电话）||
|4-138|上行数据||
|138-276|保护带||
|276-1041|下行数据||

---
```

#### 解答

**ADSL频段分配**

|频段|范围(kHz)|用途|
|---|---|---|
|0-4|POTS（普通电话）||
|4-138|上行数据||
|138-276|保护带||
|276-1041|下行数据||

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 5: Mit welchem M.Verfahren arbeitet Splitter?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Mit welchem M.Verfahren arbeitet Splitter?  


#### 题目中文翻译 / 中文题意

分离器使用什么复用方法？

#### 德文原题

```text
### Mit welchem M.Verfahren arbeitet Splitter?

**分离器使用什么复用方法？**
```

#### 解答

**Lösung / 答案：** **FDM (Frequency Division Multiplexing) / 频分复用**

Splitter将语音信号（低频）与数据信号（高频）分离。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 6: Wie groß ist Bandbreite für ADSL?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Wie groß ist Bandbreite für ADSL?  


#### 题目中文翻译 / 中文题意

ADSL的带宽是多少？

#### 德文原题

```text
### Wie groß ist Bandbreite für ADSL?

**ADSL的带宽是多少？**
```

#### 解答

**Lösung / 答案：**

- **下行：** 约1-8 Mbit/s（可达24 Mbit/s for ADSL2+）
- **上行：** 约0.5-1 Mbit/s

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 7: Bei ISDN-Tels → Sprachdatenrate 8kHz; Abtastwerte anhand Tabelle ITU-T Spec. G.711 auf 8Bit. Was ist min Übertragrate?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt Bei ISDN-Tels → Sprachdatenrate 8kHz; Abtastwerte anhand Tabelle ITU-T Spec. G.711 auf 8Bit. Was ist min Übertragrate?  


#### 题目中文翻译 / 中文题意

ISDN电话 → 语音采样率8kHz，根据ITU-T G.711每样本8位。最小传输率是多少？

#### 德文原题

```text
### Bei ISDN-Tels → Sprachdatenrate 8kHz; Abtastwerte anhand Tabelle ITU-T Spec. G.711 auf 8Bit. Was ist min Übertragrate?

**ISDN电话 → 语音采样率8kHz，根据ITU-T G.711每样本8位。最小传输率是多少？**
```

#### 解答

**Lösung / 答案：**  
8000 samples/s×8 bit/sample=64 kbit/s8000 \text{ samples/s} \times 8 \text{ bit/sample} = \textbf{64 kbit/s}8000 samples/s×8 bit/sample=64 kbit/s

这就是ISDN的一个B信道（B-Kanal）的速率。

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 8: PPP: 2 Unterprotokolle + Funktionen?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt PPP: 2 Unterprotokolle + Funktionen?  


#### 题目中文翻译 / 中文题意

PPP的两个子协议及功能？

#### 德文原题

```text
### PPP: 2 Unterprotokolle + Funktionen?

**PPP的两个子协议及功能？**
```

#### 解答

**Lösung / 答案：**

|子协议|功能|
|---|---|
|**LCP (Link Control Protocol)**|建立、配置、测试数据链路连接|
|**NCP (Network Control Protocol)**|配置不同的网络层协议（如IPCP for IP）|

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 9: ADSL协议栈计算题

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt ADSL协议栈计算题  


#### 题目中文翻译 / 中文题意

场景： 通过ATM模拟以太网的ADSL，使用PPPoE隧道建立TCP/IP连接
协议栈：
TCP/UDP (20B/8B) | 用户数据
ATM-H (5B) | ATM负载 (48B) ...

#### 德文原题

```text
### ADSL协议栈计算题

**场景：** 通过ATM模拟以太网的ADSL，使用PPPoE隧道建立TCP/IP连接

**协议栈：**

```
TCP/UDP (20B/8B) | 用户数据
----------------
IP-H (20B)
----------------
PPP-H (5B) | FCS (4B)
----------------
ETH-H (14B) | FCS
----------------
ATM-H (5B) | ATM负载 (48B) ...
```

---
```

#### 解答

**ADSL协议栈计算题**

**场景：** 通过ATM模拟以太网的ADSL，使用PPPoE隧道建立TCP/IP连接

**协议栈：**

```
TCP/UDP (20B/8B) | 用户数据
----------------
IP-H (20B)
----------------
PPP-H (5B) | FCS (4B)
----------------
ETH-H (14B) | FCS
----------------
ATM-H (5B) | ATM负载 (48B) ...
```

---

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---

### 题目 10: TCP 438B UD versenden → Wie viel Byte an Daten von ATM-Schicht sind an ADSL-Schicht zu übergeben?

**类型：** 考卷  

**来源说明：** Klausur 2015, Abschnitt TCP 438B UD versenden → Wie viel Byte an Daten von ATM-Schicht sind an ADSL-Schicht zu übergeben?  


#### 题目中文翻译 / 中文题意

发送TCP 438B用户数据 → ATM层向ADSL层传递多少字节？

#### 德文原题

```text
### TCP 438B UD versenden → Wie viel Byte an Daten von ATM-Schicht sind an ADSL-Schicht zu übergeben?

**发送TCP 438B用户数据 → ATM层向ADSL层传递多少字节？**
```

#### 解答

**Lösung / 答案：**

**计算过程：**

1. TCP层：438B + 20B(TCP头) = 458B
    
2. IP层：458B + 20B(IP头) = 478B
    
3. PPP层：478B + 5B(PPP头) = 483B
    
4. 以太网层：483B + 14B(ETH头) + 4B(FCS) = 501B
    
5. ATM分片：
    
    - ATM负载 = 48B
    - 需要的ATM信元 = ⌈501/48⌉ = 11个信元
    - 每个信元 = 5B(头) + 48B(负载) = 53B
6. **总计 = 11 × 53 = 583B**
    

---

**总结 / Zusammenfassung**

这份手写笔记涵盖了计算机网络的核心主题：

|主题|重点内容|
|---|---|
|OSI/TCP-IP模型|各层功能对比|
|物理层|电缆类型、光纤、信号传播|
|数据链路层|CSMA/CD、错误检测|
|网络层|IP、路由、多播|
|传输层|TCP拥塞/流量控制|
|应用层|DNS、HTTP等协议|
|接入技术|ADSL、ISDN、协议栈|

**备注：** 本题归入本知识点；题目中文翻译/题意、德文原题、解答和来源已保留。


---
