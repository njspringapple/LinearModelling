# Uebungsblatt 1 - 中德对照解答

## 1. Grundlagen / 分布式系统基础

**DE Aufgabenidee:** Zwei Gelehrte spielen Schach per Brieftauben. Gesucht sind Stoerfaelle und Gegenmassnahmen.

**中文题意：** 用“信鸽下棋”类比分布式系统通信，列出问题、后果和避免方法。

**中文解题思路：** 这道题不是要谈真实的鸽子，而是要把“远程下棋”抽象成两个分布式节点通过不可靠信道交换消息。只要消息可能丢失、延迟、重复、乱序或损坏，就会出现和分布式系统类似的问题。

| Stoerfall / 问题 | Folge / 后果 | Vermeidung / 缓解 |
|---|---|---|
| Brief geht verloren / 消息丢失 | Zug kommt nie an / 对方不知道走法 | Quittung, Timeout, erneutes Senden / ACK、超时重传 |
| Brief kommt doppelt an / 重复消息 | Zug wird evtl. zweimal ausgefuehrt / 可能重复执行 | Sequenznummern / 序列号 |
| Briefe kommen vertauscht an / 乱序 | falsche Spielstellung / 状态错误 | Nummerierte Zuege, nur erwartete Nummer akzeptieren / 按序号接收 |
| Brief wird verfaelscht / 内容损坏 | falscher Zug / 错误执行 | Pruefsumme, Signatur, Plausibilitaetscheck / 校验和或签名 |
| Sehr lange Laufzeit / 延迟大 | beide warten, unklarer Zustand / 双方等待 | Fristen, Statusanfragen, erneute Synchronisation / 超时与状态同步 |

**Wissen / 知识点：** 分布式系统的核心困难是“不共享内存、不共享时钟、通信不可靠”。因此协议需要确认、重传、编号、校验和状态同步。

## 2. Von Netzen und Baeumen / 图和树

**(a) Vollvermaschter Graph mit 4 Knoten / 4 个点全互联**

**DE:** Ein vollvermaschter Graph mit 4 Knoten ist der vollstaendige Graph `K4`.

**中文：** 4 个节点的完全图 `K4` 有 `4*3/2 = 6` 条边，每个节点度数为 3。

**(b) Baeume mit 4 Knoten / 4 个节点树的形状**

**DE:** Die vier Formen kann man nach der Hoehe und der Verteilung der Kinder beschreiben.

**中文：** 四种有根树形态可按高度描述：

1. Hoehe 3: 一条链 `root - a - b - c`。
2. Hoehe 2: 根有一个孩子，该孩子有两个孩子。
3. Hoehe 2: 根有两个孩子，其中一个孩子再有一个孩子。
4. Hoehe 1: 根直接连接三个叶子，星形结构。

**(c) Binaerbaeume / 二叉树**

**DE:** Ein Binaerbaum hat pro Knoten hoechstens zwei Kinder.

**中文：** 二叉树中每个节点最多有 2 个孩子，因此最少节点数对应“一条链”，最多节点数对应“每层都满”。

1. Hoehe 5 mindestens: 最少是一条长度为 5 的链，所以有 `5 + 1 = 6` 个节点。
2. Hoehe 2 hoechstens: 满二叉树节点数 `1 + 2 + 4 = 7`。

**Wissen / 知识点：** 高度为 `h` 的二叉树最少 `h+1` 个节点，最多 `2^(h+1)-1` 个节点；高度为 `h` 的满二叉树有 `2^h` 个叶子。

## 3. Adressierung in Baeumen / 树中的寻址

![Blatt 01 Seite 2: Binaerbaeume und Praefixbaum](pictures/blatt-01_page-2-2.png)

**(a)** 图中两个二叉树高度为 3，因此叶子数为 `2^3 = 8`。

**(b)** Abstand 3 bedeutet: 从根到该节点经过 3 条边。  
Bei benannten Knoten / 若节点命名：路径是经过的节点名序列。  
Bei benannten Kanten / 若边命名：路径是 3 位二进制串：

```text
000, 001, 010, 011, 100, 101, 110, 111
```

**(c)**  
**DE:** Aus dem Pfad `1` allein kann man nicht allgemein erkennen, ob es ein innerer Knoten oder ein Blatt ist. Das haengt von der Hoehe des Baums ab.  
**中文：** 只看路径 `1` 不能判断这个节点是不是叶子，必须知道整棵树的高度。如果树高度为 1，它是叶子；如果树更高，它可能还是内部节点。

**(d)** Ein beliebiger Knoten mit Abstand 8 kann z.B. so adressiert werden:

```text
10110010
```

**(e)**  
**DE:** Alle Nachfahren eines Knotens mit Pfad `p` haben `p` als Praefix.  
**中文：** 路径为 `p` 的节点，其所有后代地址都以 `p` 开头。这就是“前缀表示一个子树”的含义。

**(f)** Pfad:

```text
1100000010101000
```

| Format | Wert |
|---|---|
| Dezimal | 49320 |
| Hexadezimal | C0A8 |
| Tupel zweier 8-Bit-Haelften | `(192, 168)` |

**Wissen / 知识点：** 树路径与二进制前缀完全对应，这正是 IP 前缀、路由聚合、CIDR 的直观基础。

## 4. Zahlen sind auch nur Baeume / 数字也是树

高度 4 的完整二叉树有 `2^4 = 16` 个叶子。固定长度为 `l` 的前缀后，还剩 `4-l` 位可变。

**DE:** Ein Praefix legt die ersten Bits des Pfads fest; alle verbleibenden Bits koennen frei gewaehlt werden.

**中文：** 前缀固定了路径最前面的若干位，后面的位可以任意变化，所以叶子数量按 `2^(总高度 - 前缀长度)` 计算。

**(a)** Praefix `0`: `2^(4-1) = 8` Pfade.

**(b)** Bei Hoehe `n`: `2^(n-1)` Pfade.

**(c)** Fuer Hoehe 4 gilt:

| Teilbaum | Praefixlaenge | Blaetter |
|---|---:|---:|
| `1111/1` | 1 | 8 |
| `0000/2` | 2 | 4 |
| `1010/3` | 3 | 2 |
| `0101/4` | 4 | 1 |
| `0110/0` | 0 | 16 |

**Wissen / 知识点：** 前缀越长，子树越小；高度固定为 `n` 时，前缀长度 `l` 对应 `2^(n-l)` 个叶子。
