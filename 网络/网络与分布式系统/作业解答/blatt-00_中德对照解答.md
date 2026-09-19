# Uebungsblatt 0 - 中德对照解答

## 1. Anforderungen des Internets / 互联网的连接需求

**DE Aufgabenidee:** In einem vollvermaschten Netz ist jedes Endgeraet direkt mit jedem anderen verbunden. Danach wird eine hierarchische Struktur betrachtet, in der jeweils hoechstens 5 Geraete/Knoten an einen Knoten angeschlossen sind.

**中文题意：** 比较全互联网络和 5 叉层次结构所需连接数。

### Loesung / 解答

**(a) Vollvermaschung / 全互联**

**DE:** Bei `N` Teilnehmern braucht man fuer jedes ungeordnete Paar genau eine Verbindung:

**中文：** 有 `N` 个参与者时，每两个不同参与者之间需要一条连接，也就是从 `N` 个点中任选 2 个点：

```text
K_N = N * (N - 1) / 2
```

| Teilnehmer | Verbindungen |
|---:|---:|
| 8 | 28 |
| 300 | 44850 |
| N | `N(N-1)/2` |

**(b) Hierarchie mit maximal 5 Kindern / 每层最多 5 个子节点**

**DE:** Wenn `N = 5^a` gilt, hat ein vollstaendiger 5-aerer Baum mit `N` Blaettern:

**中文：** 如果 `N = 5^a`，可以看成一棵满 5 叉树。最底层有 `N` 个终端，上一层有 `N/5` 个节点，再上一层有 `N/25` 个节点，以此类推。树中每个非根节点对应一条到父节点的连接：

```text
N + N/5 + N/25 + ... + 5 = 5(N - 1)/4
```

**DE:** Das Ergebnis ist die Zahl der Verbindungen.

**中文：** 所以上式就是所需连接总数。对于不是精确 5 的幂的情况，例如 8 和 300，可以按每层向上取整估算。

| Teilnehmer | Verbindungen |
|---:|---:|
| 8 | `8 + ceil(8/5) = 10` |
| 300 | `300 + 60 + 12 + 3 = 375` |
| `N = 5^a` | `5(N-1)/4` |

**Wissen / 知识点：** 全互联的边数是二次增长 `O(N^2)`；树形/层次结构近似线性增长 `O(N)`，所以互联网必须分层、聚合和路由，而不可能让所有终端物理直连。

## 2. Das Stellenwertsystem / 位值计数系统

**(a) Hexadezimalziffern / 十六进制数字**

```text
0 1 2 3 4 5 6 7 8 9 A B C D E F
```

**(b) 2, 4, 8, 10 in verschiedenen Basen / 不同进制表示**

| Dezimal | Hex | Oktal | Binaer |
|---:|---:|---:|---:|
| 2 | 2 | 2 | 10 |
| 4 | 4 | 4 | 100 |
| 8 | 8 | 10 | 1000 |
| 10 | A | 12 | 1010 |

**(c) Umrechnung / 转换**

| Dezimal | Binaer | Hex |
|---:|---:|---:|
| 16 | 10000 | 10 |
| 127 | 1111111 | 7F |
| 168 | 10101000 | A8 |
| 172 | 10101100 | AC |
| 192 | 11000000 | C0 |
| 255 | 11111111 | FF |

**(d) `2^32 - 1`**

**DE:** `2^32 - 1` ist binaer 32-mal die Ziffer `1`.

**中文：** `2^32` 在二进制中是 `1` 后面 32 个 `0`，减去 1 后就变成 32 个连续的 `1`：

```text
11111111111111111111111111111111
```

Also: 32 Stellen, 32 Einsen, 0 Nullen.

**Wissen / 知识点：** 十六进制每位对应 4 个二进制位；IPv4 地址常用 8 位分组，因此 `255 = 11111111 = FF` 很常见。

## 3. Rechnen in unterschiedlichen Zahlensystemen / 不同进制计算

### (a) Multiplikationstabelle / 乘法表

**中文说明：** 表中“Zahl dez.”是原数的十进制值，结果分别写成二进制和十六进制。比如十进制 16 乘以 8 得到 128，二进制是 `10000000`，十六进制是 `80`。

| Zahl dez. | Faktor | Binaer-Ergebnis | Hex-Ergebnis |
|---:|---:|---:|---:|
| 1 | 2 | 10 | 2 |
| 2 | 2 | 100 | 4 |
| 8 | 2 | 10000 | 10 |
| 10 | 2 | 10100 | 14 |
| 16 | 2 | 100000 | 20 |
| 1 | 8 | 1000 | 8 |
| 2 | 8 | 10000 | 10 |
| 8 | 8 | 1000000 | 40 |
| 10 | 8 | 1010000 | 50 |
| 16 | 8 | 10000000 | 80 |

### (b) Potenzen / 幂

**中文说明：** 在某个基数 `b` 的进制中，`b^k` 的表示就是 `1` 后面跟 `k` 个 `0`。因此四行看起来形式相同，只是每一行所在的进制不同。

| Ausdruck | Ergebnis im geforderten System |
|---|---|
| `2^2 ... 2^7` binaer | `100, 1000, 10000, 100000, 1000000, 10000000` |
| `8^2 ... 8^7` oktal | `100, 1000, 10000, 100000, 1000000, 10000000` |
| `10^2 ... 10^7` dezimal | `100, 1000, 10000, 100000, 1000000, 10000000` |
| `16^2 ... 16^7` hex | `100, 1000, 10000, 100000, 1000000, 10000000` |

**Wissen / 知识点：** 在基数为 `b` 的系统中，`b^k` 写成 `1` 后面跟 `k` 个零。
