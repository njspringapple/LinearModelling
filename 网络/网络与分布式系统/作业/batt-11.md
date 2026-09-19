## Aufgabe 1: Interpretation einer DNS-Antwort (H)

## 第1题：DNS 响应的解读（家庭作业）

Ein nützliches Diagnosewerkzeug für den DNS ist das Programm `dig`, das auf vielen Unix-Derivaten vorhanden ist.

DNS 的一个有用的诊断工具是 `dig` 程序，它存在于许多 Unix 衍生系统上。

### (a) Skizze des DNS-Verkehrs | DNS 流量示意图

**Aufgabe | 题目：** Zeichnen Sie eine Skizze, die den DNS-Verkehr zur Anfrage darstellt.

绘制一个展示 DNS 查询流量的示意图。

**Lösung | 解答：**

```
┌─────────────────────────────────────────────────────────────────────────────────┐
│                                                                                 │
│   Anfragender Host                     Lokaler DNS-Server                       │
│   (查询主机)                            (本地 DNS 服务器)                         │
│   IP: unbekannt                        IP: 192.168.218.30                       │
│                                                                                 │
│         │                                      │                                │
│         │────── Query: mail.nm.ifi.lmu.de ────►│                                │
│         │       (查询)                          │                                │
│         │                                      │                                │
│         │                                      │                                │
│         │                              ┌───────┴───────┐                        │
│         │                              │               │                        │
│         │                              ▼               │                        │
│         │                   d.root-servers.net         │                        │
│         │                   IP: 128.8.10.90            │                        │
│         │                   (根服务器)                  │                        │
│         │                              │               │                        │
│         │                   Query: NS für .de?         │                        │
│         │                   (查询 .de 的 NS)            │                        │
│         │                              │               │                        │
│         │                   Response: NS Records       │                        │
│         │                   (响应: NS 记录)             │                        │
│         │                              │               │                        │
│         │                              ▼               │                        │
│         │                        C.DE.NET              │                        │
│         │                   IP: 208.48.81.43           │                        │
│         │                   (.de 域名服务器)            │                        │
│         │                              │               │                        │
│         │                   Query: NS für lmu.de?      │                        │
│         │                   (查询 lmu.de 的 NS)         │                        │
│         │                              │               │                        │
│         │                   Response: NS Records       │                        │
│         │                              │               │                        │
│         │                              ▼               │                        │
│         │                dns3.lrz-muenchen.de          │                        │
│         │                   IP: 129.187.5.2            │                        │
│         │                   (lmu.de 域名服务器)         │                        │
│         │                              │               │                        │
│         │                   Query: mail.nm.ifi.lmu.de? │                        │
│         │                              │               │                        │
│         │                   Response: CNAME + A Record │                        │
│         │                   (响应: CNAME + A 记录)      │                        │
│         │                              │               │                        │
│         │◄───── Response: 141.84.218.30 ──────────────┘                        │
│         │       (响应)                                                          │
│                                                                                 │
└─────────────────────────────────────────────────────────────────────────────────┘
```

**Beteiligte Server | 涉及的服务器：**

|Server|Hostname|IP-Adresse|Rolle|
|---|---|---|---|
|本地 DNS 服务器|-|192.168.218.30|Lokaler DNS-Server|
|根服务器|d.root-servers.net|128.8.10.90|Root Nameserver|
|.de 服务器|C.DE.NET|208.48.81.43|TLD Nameserver|
|lmu.de 服务器|dns3.lrz-muenchen.de|129.187.5.2|Autoritativer NS für lmu.de|

---

### (b) Rekursiv oder iterativ? | 递归还是迭代？

**Aufgabe | 题目：** Ist die Anfrage rekursiv oder iterativ?

这个查询是递归的还是迭代的？

**Lösung | 解答：**

Die Anfrage ist **iterativ**.

这个查询是**迭代的**。

**Begründung | 理由：**

Der lokale DNS-Server (192.168.218.30) fragt nacheinander verschiedene Server an und erhält jeweils nur Verweise auf den nächsten zuständigen Server:

本地 DNS 服务器 (192.168.218.30) 依次查询不同的服务器，每次只获得指向下一个负责服务器的引用：

1. Root-Server → liefert NS-Records für `.de` | 根服务器 → 返回 `.de` 的 NS 记录
2. `.de`-Server → liefert NS-Records für `lmu.de` | `.de` 服务器 → 返回 `lmu.de` 的 NS 记录
3. `lmu.de`-Server → liefert die endgültige Antwort | `lmu.de` 服务器 → 返回最终答案

Die Option `+trace` zeigt genau diesen iterativen Ablauf.

`+trace` 选项正是显示了这种迭代过程。

---

### (c) Anfrage an den Root-Server | 对根服务器的查询

**Aufgabe | 题目：** Die Ausgabe enthält eine Anfrage an einen der DNS-Root-Server. Wonach wird er gefragt?

输出包含对 DNS 根服务器的查询。查询的是什么？

**Lösung | 解答：**

Der Root-Server wird nach den **autoritativen Nameservern für die Top-Level-Domain `.de`** gefragt.

根服务器被查询的是**顶级域 `.de` 的权威名称服务器**。

Die Antwort (Zeilen 18-23) enthält NS-Records für die `.de`-Domain:

响应（第 18-23 行）包含 `.de` 域的 NS 记录：

- C.DE.NET
- L.DE.NET
- F.NIC.de
- S.DE.NET
- A.NIC.de
- Z.NIC.de

---

### (d) Alias-Auflösung | 别名解析

**Aufgabe | 题目：** Der gesuchte Rechnername `mail.nm.ifi.lmu.de` ist ein Alias.

查询的计算机名 `mail.nm.ifi.lmu.de` 是一个别名。

#### i. Wie heißt die Maschine wirklich? | 这台机器的真实名称是什么？

**Lösung | 解答：**

Die Maschine heißt wirklich **`pcheger0.nm.ifi.lmu.de`** (Zeile 31).

这台机器的真实名称是 **`pcheger0.nm.ifi.lmu.de`**（第 31 行）。

```
mail.nm.ifi.lmu.de. 86400 IN CNAME pcheger0.nm.ifi.lmu.de.
```

Der CNAME-Record (Canonical Name) zeigt, dass `mail.nm.ifi.lmu.de` ein Alias für `pcheger0.nm.ifi.lmu.de` ist.

CNAME 记录（规范名称）表明 `mail.nm.ifi.lmu.de` 是 `pcheger0.nm.ifi.lmu.de` 的别名。

#### ii. Welche IP-Adresse hat sie? | 它的 IP 地址是什么？

**Lösung | 解答：**

Die IP-Adresse ist **141.84.218.30** (Zeile 32).

IP 地址是 **141.84.218.30**（第 32 行）。

```
pcheger0.nm.ifi.lmu.de. 86400 IN A 141.84.218.30
```

---

### (e) Weitere Aussagen zu DNS-Servern | 关于 DNS 服务器的更多信息

#### i. Wer betreibt die DNS-Server für lmu.de? | 谁运营 lmu.de 的 DNS 服务器？

**Lösung | 解答：**

Die DNS-Server für die Domäne `lmu.de` werden vom **LRZ (Leibniz-Rechenzentrum)** betrieben.

`lmu.de` 域的 DNS 服务器由 **LRZ（莱布尼茨计算中心）** 运营。

Zeilen 26-28 zeigen:

第 26-28 行显示：

- dns3.lrz-muenchen.de
- dns1.lrz-muenchen.de
- dns2.lrz-muenchen.de

---

#### ii. Welche DNS-Server für die Domäne der gesuchten Maschine? | 哪些 DNS 服务器可以响应所查询机器所在域的查询？

**Lösung | 解答：**

Für die Domäne `nm.ifi.lmu.de` sind folgende DNS-Server zuständig (Zeilen 33-38):

以下 DNS 服务器负责 `nm.ifi.lmu.de` 域（第 33-38 行）：

|Server|说明|
|---|---|
|acheron.ifi.lmu.de|IFI-Server|
|dns3.lrz-muenchen.de|LRZ-Server|
|dns1.nm.ifi.lmu.de|NM-Team Server|
|dns1.lrz-muenchen.de|LRZ-Server|
|dns2.lrz-muenchen.de|LRZ-Server|
|dns0.nm.ifi.lmu.de|NM-Team Server|

---

#### iii. Autoritative Antwort? | 权威响应？

**Aufgabe | 题目：** Wurde die gesuchte IP-Adresse von einem autoritativen Server geliefert?

所查询的 IP 地址是否由权威服务器提供？

**Lösung | 解答：**

**Ja**, die IP-Adresse wurde von einem autoritativen Server geliefert.

**是的**，IP 地址是由权威服务器提供的。

**Begründung | 理由：**

Der Server `dns3.lrz-muenchen.de` (129.187.5.2), der die endgültige Antwort geliefert hat (Zeile 39), ist in den NS-Records für `nm.ifi.lmu.de` aufgeführt (Zeile 34). Er ist daher ein autoritativer Nameserver für diese Zone.

提供最终答案的服务器 `dns3.lrz-muenchen.de` (129.187.5.2)（第 39 行）被列在 `nm.ifi.lmu.de` 的 NS 记录中（第 34 行）。因此，它是该区域的权威名称服务器。

---

### (f) DNS-Cache zur Analyse von Nutzerverhalten | 利用 DNS 缓存分析用户行为

**Aufgabe | 题目：** Gibt es eine Möglichkeit, die von Nutzern meist besuchten Web-Server ausfindig zu machen?

是否有可能找出用户最常访问的 Web 服务器？

**Lösung | 解答：**

**Ja**, es gibt mehrere Möglichkeiten:

**是的**，有几种方法：

1. **Analyse der Cache-Einträge | 缓存条目分析：** Häufig aufgerufene Domains erscheinen mit aktuellen TTL-Werten im Cache. Je öfter ein Eintrag erneuert wird, desto beliebter ist die Website.
    
    经常访问的域名会以当前的 TTL 值出现在缓存中。条目刷新越频繁，网站越受欢迎。
    
2. **Statistik über Cache-Hits | 缓存命中统计：** Durch Logging der DNS-Anfragen kann ermittelt werden, welche Domains am häufigsten angefragt werden.
    
    通过记录 DNS 查询，可以确定哪些域名被查询得最频繁。
    
3. **TTL-Analyse | TTL 分析：** Einträge, die kurz vor Ablauf der TTL bereits wieder im Cache sind, wurden vermutlich erneut angefragt.
    
    在 TTL 到期前就再次出现在缓存中的条目，可能已被再次查询。
    

**Einschränkung | 限制：** Dies ermöglicht nur die Identifikation von Domains, nicht von spezifischen Webseiten oder Inhalten.

这只能识别域名，而不是具体的网页或内容。

---

## Aufgabe 2: HTTP Requests und Response (H)

## 第2题：HTTP 请求和响应（家庭作业）

### HTTP GET-Request | HTTP GET 请求

```http
GET /gnu/gnu.html HTTP/1.1
Host: www.gnu.org
User-Agent: Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:67.0) Gecko/20100101 Firefox/67.0
Accept: text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8
Accept-Language: de-DE,en-US;q=0.7,en;q=0.3
Accept-Encoding: gzip, deflate, br
Connection: keep-alive
```

### (a) URL des angefragten Dokuments | 所请求文档的 URL

**Lösung | 解答：**

Die URL lautet: **`http://www.gnu.org/gnu/gnu.html`**

URL 是：**`http://www.gnu.org/gnu/gnu.html`**

**Zusammensetzung | 组成：**

- Protokoll: `http://` (implizit durch HTTP/1.1) | 协议
- Host: `www.gnu.org` (aus dem Host-Header) | 主机
- Pfad: `/gnu/gnu.html` (aus der Request-Zeile) | 路径

---

### (b) HTTP-Version | HTTP 版本

**Lösung | 解答：**

Der Browser nutzt **HTTP/1.1**.

浏览器使用 **HTTP/1.1**。

Dies ist in der ersten Zeile des Requests zu sehen: `GET /gnu/gnu.html HTTP/1.1`

这在请求的第一行可以看到：`GET /gnu/gnu.html HTTP/1.1`

---

### (c) Persistente oder nicht-persistente Verbindung? | 持久连接还是非持久连接？

**Lösung | 解答：**

Der Browser fragt eine **persistente Verbindung** an.

浏览器请求**持久连接**。

**Begründung | 理由：**

Der Header `Connection: keep-alive` zeigt, dass der Browser die TCP-Verbindung für weitere Anfragen offen halten möchte.

`Connection: keep-alive` 头部表明浏览器希望保持 TCP 连接打开以进行更多请求。

---

### (d) IP-Adresse des Hosts | 主机的 IP 地址

**Lösung | 解答：**

Die IP-Adresse des Hosts, auf dem der Browser ausgeführt wird, ist **aus dem HTTP-Request nicht ersichtlich**.

运行浏览器的主机的 IP 地址**无法从 HTTP 请求中看出**。

**Begründung | 理由：**

HTTP ist ein Anwendungsschichtprotokoll und enthält keine IP-Adressen. Die IP-Adresse ist Teil des IP-Headers in der Vermittlungsschicht, nicht des HTTP-Requests.

HTTP 是应用层协议，不包含 IP 地址。IP 地址是网络层 IP 头部的一部分，而不是 HTTP 请求的一部分。

Um die IP zu ermitteln, müsste man den IP-Header der darunterliegenden Schicht analysieren (z.B. mit Wireshark).

要确定 IP 地址，需要分析底层的 IP 头部（例如使用 Wireshark）。

---

### (e) Browser-Typ und Zweck | 浏览器类型及用途

**Lösung | 解答：**

**Browser-Typ | 浏览器类型：**

- **Firefox 67.0** auf Ubuntu Linux (x86_64)
- Firefox 67.0，运行在 Ubuntu Linux (x86_64) 上

Der User-Agent-String lautet:

```
Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:67.0) Gecko/20100101 Firefox/67.0
```

**Zweck der Übermittlung | 传输目的：**

1. **Content Negotiation | 内容协商：** Der Server kann je nach Browser-Fähigkeiten unterschiedliche Inhalte ausliefern (z.B. optimierte Seiten für mobile Browser).
    
    服务器可以根据浏览器能力提供不同的内容（例如，为移动浏览器优化的页面）。
    
2. **Statistiken | 统计：** Server-Betreiber können analysieren, welche Browser ihre Nutzer verwenden.
    
    服务器运营商可以分析用户使用哪些浏览器。
    
3. **Kompatibilität | 兼容性：** Ältere Webseiten können browser-spezifischen Code ausliefern.
    
    旧网站可以提供特定于浏览器的代码。
    

**Notwendigkeit | 必要性：**

**Nein**, die Übermittlung ist technisch nicht notwendig für die HTTP-Kommunikation. Sie ist optional, aber in der Praxis fast immer vorhanden.

**不是的**，传输在技术上对于 HTTP 通信来说不是必需的。它是可选的，但在实践中几乎总是存在。

---

### HTTP Response | HTTP 响应

```http
HTTP/1.1 200 OK
Date: Thu, 23 May 2019 08:27:34 GMT
Server: Apache/2.4.7
Content-Location: gnu.html
Accept-Ranges: bytes
Content-Encoding: gzip
Content-Length: 5751
Keep-Alive: timeout=3, max=98
Connection: Keep-Alive
Content-Type: text/html
Content-Language: en

<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN
...
```

### (f) Dokument gefunden? Zeit? | 文档找到了吗？时间？

**Lösung | 解答：**

**Ja**, der Server konnte das Dokument erfolgreich finden.

**是的**，服务器成功找到了文档。

**Begründung | 理由：** Der Statuscode **200 OK** zeigt, dass die Anfrage erfolgreich war.

状态码 **200 OK** 表明请求成功。

**Zeit der Generierung | 生成时间：**

Die Antwort wurde generiert am: **Donnerstag, 23. Mai 2019 um 08:27:34 GMT**

响应生成于：**2019年5月23日星期四 08:27:34 GMT**

(aus dem `Date`-Header | 来自 `Date` 头部)

---

### (g) Sprache der Antwort | 响应的语言

**Lösung | 解答：**

Die Antwort ist in **Englisch (en)** formuliert.

响应是用**英语 (en)** 编写的。

Dies zeigt der Header: `Content-Language: en`

这由头部显示：`Content-Language: en`

---

### (h) Größe des Dokuments | 文档大小

**Lösung | 解答：**

Das zurückgegebene Dokument enthält **5751 Bytes**.

返回的文档包含 **5751 字节**。

Dies zeigt der Header: `Content-Length: 5751`

这由头部显示：`Content-Length: 5751`

**Hinweis | 注意：** Dies ist die Größe der komprimierten Daten (`Content-Encoding: gzip`), nicht die unkomprimierte Größe.

这是压缩数据的大小（`Content-Encoding: gzip`），而不是未压缩的大小。

---

### (i) Erste 5 Bytes und persistente Verbindung | 前 5 个字节和持久连接

**Lösung | 解答：**

**Erste 5 Bytes des Dokuments | 文档的前 5 个字节：**

Die ersten 5 Bytes sind: **`<!DOC`**

前 5 个字节是：**`<!DOC`**

(Beginn von `<!DOCTYPE html...>`)

**Persistente Verbindung bestätigt? | 持久连接已确认？**

**Ja**, der Server hat die Anfrage nach einer persistenten Verbindung bestätigt.

**是的**，服务器已确认持久连接请求。

**Begründung | 理由：**

- `Connection: Keep-Alive` bestätigt die persistente Verbindung | 确认持久连接
- `Keep-Alive: timeout=3, max=98` gibt die Parameter an: | 给出参数：
    - Timeout von 3 Sekunden | 超时 3 秒
    - Maximal 98 weitere Anfragen auf dieser Verbindung | 此连接上最多还有 98 个请求

---

## Aufgabe 3: Email (H)

## 第3题：电子邮件（家庭作业）

```
┌─────────────┐                                      ┌─────────────┐
│   Sender    │                                      │  Empfänger  │
│  (发送者)    │                                      │   (接收者)   │
├─────────────┤                                      ├─────────────┤
│ User Agent  │                                      │ User Agent  │
└──────┬──────┘                                      └──────▲──────┘
       │                                                    │
       │ 1.                                            3.   │
       ▼                                                    │
┌─────────────┐              2.               ┌─────────────┐
│ Mail Server │ ─────────────────────────────►│ Mail Server │
│ des Senders │                               │des Empfängers│
│ (发送方邮件  │                               │ (接收方邮件  │
│   服务器)    │                               │   服务器)    │
└─────────────┘                               └─────────────┘
```

### (a) Protokolle auf den drei Übertragungswegen | 三条传输路径上的协议

**Aufgabe | 题目：** Welche Protokolle der Anwendungsschicht können auf den drei eingezeichneten Übertragungswegen eingesetzt werden?

在三条标注的传输路径上可以使用哪些应用层协议？

**Lösung | 解答：**

| Weg                                | Protokolle             | 协议          | Erklärung                                           | 说明                  |
| ---------------------------------- | ---------------------- | ----------- | --------------------------------------------------- | ------------------- |
| **1** (UA → Mail Server Sender)    | **SMTP**               | SMTP        | User Agent sendet E-Mail an den lokalen Mail Server | 用户代理将电子邮件发送到本地邮件服务器 |
| **2** (Mail Server → Mail Server)  | **SMTP**               | SMTP        | Mail Transfer Agents kommunizieren untereinander    | 邮件传输代理之间相互通信        |
| **3** (Mail Server → UA Empfänger) | **POP3** oder **IMAP** | POP3 或 IMAP | User Agent ruft E-Mails vom Server ab               | 用户代理从服务器获取电子邮件      |

---

### (b) Webbasierter E-Mail Account | 基于 Web 的电子邮件账户

**Aufgabe | 题目：** Welche zusätzlichen Protokolle sind involviert, wenn der Sender einen webbasierten E-Mail Account (z.B. GMail, GMX) verwendet?

如果发送者使用基于 Web 的电子邮件账户（如 GMail、GMX），会涉及哪些额外的协议？

**Lösung | 解答：**

**Zusätzliche Protokolle | 额外的协议：**

|Weg|Klassisch|Mit Webmail|
|---|---|---|
|**1**|SMTP|**HTTP/HTTPS**|
|**2**|SMTP|SMTP (unverändert|
|**3**|POP3/IMAP|**HTTP/HTTPS**|

**Erklärung | 说明：**

Bei webbasierten E-Mail-Diensten:

对于基于 Web 的电子邮件服务：

1. **Weg 1 (Sender → Mail Server):** Der Browser des Senders kommuniziert über **HTTP/HTTPS** mit dem Webserver des E-Mail-Anbieters. Der Webserver reicht die E-Mail dann intern an den Mail Server weiter (via SMTP oder interne APIs).
    
    **路径 1（发送者 → 邮件服务器）：** 发送者的浏览器通过 **HTTP/HTTPS** 与电子邮件提供商的 Web 服务器通信。然后 Web 服务器在内部将电子邮件转发到邮件服务器（通过 SMTP 或内部 API）。
    
2. **Weg 2:** Bleibt SMTP zwischen den Mail-Servern.
    
    **路径 2：** 邮件服务器之间保持使用 SMTP。
    
3. **Weg 3 (Mail Server → Empfänger):** Der Empfänger greift über **HTTP/HTTPS** auf sein Webmail-Interface zu, um E-Mails zu lesen.
    
    **路径 3（邮件服务器 → 接收者）：** 接收者通过 **HTTP/HTTPS** 访问其 Webmail 界面来阅读电子邮件。
    

**Zusätzlich | 另外：**

- **DNS** wird für die MX-Record-Auflösung benötigt
- **TLS/SSL** für verschlüsselte Verbindungen

---

### (c) Message Transfer System | 消息传输系统

**Aufgabe | 题目：** Welche dargestellten Systeme sind Teil des Message Transfer Systems?

图中哪些系统属于消息传输系统？

**Lösung | 解答：**

Zum **Message Transfer System (MTS)** gehören:

属于**消息传输系统 (MTS)** 的有：

- **Mail Server des Senders** | 发送方的邮件服务器
- **Mail Server des Empfängers** | 接收方的邮件服务器

**Nicht zum MTS gehören | 不属于 MTS 的有：**

- User Agent des Senders | 发送者的用户代理
- User Agent des Empfängers | 接收者的用户代理

**Begründung | 理由：**

Das MTS besteht aus den **Message Transfer Agents (MTAs)**, die für den Transport der E-Mail zwischen den Servern verantwortlich sind. Die User Agents (MUAs) sind Teil des **Message Handling Systems (MHS)**, aber nicht des MTS.

MTS 由**消息传输代理 (MTA)** 组成，它们负责在服务器之间传输电子邮件。用户代理 (MUA) 是**消息处理系统 (MHS)** 的一部分，但不是 MTS 的一部分。

---

### (d) Dienstgüteparameter | 服务质量参数

**Aufgabe | 题目：** E-Mail ist empfindlich gegen Datenverlust. Gibt es Dienstgüteparameter, gegen die E-Mail unempfindlich ist?

电子邮件对数据丢失敏感。是否有电子邮件不敏感的服务质量参数？

**Lösung | 解答：**

**Ja**, E-Mail ist **unempfindlich** gegen folgende Dienstgüteparameter:

**是的**，电子邮件对以下服务质量参数**不敏感**：

| Parameter                        | E-Mail Sensitivität | Begründung                                                                                                                       |
| -------------------------------- | ------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| **Verzögerung (Delay)**          | Tolerant / 容忍       | E-Mails müssen nicht in Echtzeit zugestellt werden. Verzögerungen von Sekunden bis Minuten (oder sogar Stunden) sind akzeptabel. |
| **Übertragungsrate (Bandwidth)** | Elastisch / 弹性      | E-Mails haben keine festen Bandbreitenanforderungen. Sie passen sich an die verfügbare Kapazität an.                             |
| **Jitter (Schwankungen)**        | Tolerant / 容忍       | Schwankungen in der Übertragungszeit sind irrelevant, da E-Mail nicht zeitkritisch ist.                                          |

**Begründung | 理由：**

E-Mail ist eine **asynchrone** Kommunikationsform. Im Gegensatz zu Echtzeitanwendungen wie VoIP oder Videostreaming gibt es keine zeitlichen Anforderungen. Solange die Nachricht vollständig und korrekt ankommt (keine Datenverluste), ist der Zeitpunkt der Zustellung zweitrangig.

电子邮件是一种**异步**通信形式。与 VoIP 或视频流等实时应用不同，它没有时间要求。只要消息完整正确地到达（没有数据丢失），传递时间就是次要的。

---

## Aufgabe 4: Rollenwechsel

## 第4题：角色转换

**Aufgabe | 题目：** Entwerfen Sie für jedes Kapitel der Vorlesung je eine Quizfrage mit einer richtigen und mindestens zwei falschen Antwortalternativen.

为讲座的每一章设计一道测验题，包含一个正确答案和至少两个错误选项。

**Lösung | 解答：**

### Kapitel 1: Einführung und Grundlagen | 第1章：引言和基础

**Frage | 问题：** Welches Schichtenmodell beschreibt die Kommunikation im Internet am genauesten?

哪个层次模型最准确地描述了互联网通信？

- A) Das ISO/OSI-Modell mit 7 Schichten | ISO/OSI 7 层模型
- B) **Das TCP/IP-Modell mit 4 Schichten** ✓ | TCP/IP 4 层模型
- C) Das IEEE 802-Modell mit 3 Schichten | IEEE 802 3 层模型
- D) Das HTTP-Modell mit 5 Schichten | HTTP 5 层模型

---

### Kapitel 2: Bitübertragungsschicht | 第2章：物理层

**Frage | 问题：** Was beschreibt das Nyquist-Theorem?

奈奎斯特定理描述了什么？

- A) Die maximale Reichweite eines Funksignals | 无线信号的最大范围
- B) **Die maximale Symbolrate über einen bandbreitenbegrenzten Kanal** ✓ | 带宽受限信道上的最大符号率
- C) Die minimale Fehlerrate bei digitaler Übertragung | 数字传输中的最小错误率
- D) Die optimale Kabeldicke für Ethernet | 以太网的最佳电缆厚度

---

### Kapitel 3: Sicherungsschicht | 第3章：数据链路层

**Frage | 问题：** Welches Verfahren wird bei CSMA/CD verwendet, wenn eine Kollision erkannt wird?

CSMA/CD 中检测到冲突时使用什么方法？

- A) Die Übertragung wird sofort wiederholt | 立即重新传输
- B) Die Übertragung wird abgebrochen und nie wiederholt | 中止传输且永不重试
- C) **Die Übertragung wird abgebrochen und nach einer zufälligen Wartezeit wiederholt** ✓ | 中止传输并在随机等待时间后重试
- D) Die Kollision wird ignoriert und die Übertragung fortgesetzt | 忽略冲突并继续传输

---

### Kapitel 4: Vermittlungsschicht | 第4章：网络层

**Frage | 问题：** Wie viele Bits umfasst eine IPv4-Adresse?

IPv4 地址包含多少位？

- A) 16 Bits
- B) **32 Bits** ✓
- C) 64 Bits
- D) 128 Bits

---

### Kapitel 5: Transportschicht | 第5章：传输层

**Frage | 问题：** Welches Transportprotokoll garantiert eine zuverlässige, geordnete Zustellung von Daten?

哪个传输协议保证数据的可靠、有序传递？

- A) UDP
- B) **TCP** ✓
- C) IP
- D) ICMP

---

### Kapitel 6: Anwendungsschicht | 第6章：应用层

**Frage | 问题：** Welcher DNS-Record-Typ ordnet einem Hostnamen eine IPv4-Adresse zu?

哪种 DNS 记录类型将主机名映射到 IPv4 地址？

- A) MX
- B) CNAME
- C) NS
- D) **A** ✓

---

### Kapitel 7: Verteilte Systeme | 第7章：分布式系统

**Frage | 问题：** Was ist das Hauptmerkmal einer Peer-to-Peer (P2P) Architektur?

点对点 (P2P) 架构的主要特征是什么？

- A) Ein zentraler Server verwaltet alle Ressourcen | 中央服务器管理所有资源
- B) **Alle Teilnehmer sind gleichberechtigt und können sowohl als Client als auch als Server agieren** ✓ | 所有参与者平等，可以同时作为客户端和服务器
- C) Nur der Server kann Daten senden | 只有服务器可以发送数据
- D) Die Kommunikation erfolgt ausschließlich über HTTP | 通信仅通过 HTTP 进行

---

### Zusätzliche Quizfragen | 额外的测验题

**HTTP-Frage | HTTP 问题：**

Welcher HTTP-Statuscode zeigt an, dass eine Ressource dauerhaft an eine neue URL verschoben wurde?

哪个 HTTP 状态码表示资源已永久移动到新 URL？

- A) 200 OK
- B) 404 Not Found
- C) **301 Moved Permanently** ✓
- D) 500 Internal Server Error

---

**E-Mail-Frage | 电子邮件问题：**

Welches Protokoll wird typischerweise verwendet, um E-Mails zwischen Mail-Servern zu übertragen?

哪个协议通常用于在邮件服务器之间传输电子邮件？

- A) POP3
- B) IMAP
- C) **SMTP** ✓
- D) FTP