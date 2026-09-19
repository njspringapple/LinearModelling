
## 一、传输层任务与端到端思想

Internet 分层模型应用回顾：应用层调用 `send_http_request(host, message)`，下面依次是传输层、网络层和网络接入层。

Internet 分层模型应用：应用层向传输层交付 **IDU**，其中 **ICI** 为 `{host, port}`，**ID** 为 `{request}`。

IDU = ICI + ID

**发送端数据流**：传输层负责消息传输，包括连接 **建立/释放**、**是否保证交付**、**不可靠信道上的传输和拥塞管理**；网络层负责**路径**；底层负责**比特传输**。

ISO/OSI 参考模型中的传输层：从端系统视角看，传输层提供虚拟、可靠的端到端连接。

**端到端视角**：传输层实体位于端系统 A 和端系统 B 中，中间网络由多个路由器和子网组成。

**对传输层的要求**：上层应用层期望传输层在端系统之间传输消息，**可带或不带可靠性保证**。下层网络层提供**经过中转系统的消息传输**，但通常不提供可靠性保证，而是 **Best Effort**；例如过载时网络层可能丢弃分组。

不可靠网络中消息投递可能出现的问题：**消息丢失、消息重复、内容被篡改或损坏、多条消息顺序交换、传输延迟**。

**端到端论证中的可靠性保证**：包括**交付保证、抑制消息重复、安全数据传输、保证投递顺序**。可靠性保证与通信性能之间存在目标冲突，需要做成本收益权衡。

ISO/OSI 端到端视角回顾：**端系统之间的传输层提供虚拟可靠连接**。

**传输服务实现方式**：

- 无连接服务没有实体之间的连接，**实体状态无关，错误检测/纠正能力低，但同步和处理开销低、传输性能高**。

- 面向连接服务依赖实体之间的上下文，提供**更完整错误处理和可靠性，但处理开销更高**。

**无连接服务**：实体之间不建立连接，没有接收确认，发送方不知道消息是否到达目标。ICI 包含所有目的地信息，消息自行穿过网络。管理开销低，因此性能高。

## 二、面向连接服务与握手

**面向连接服务**：**建立逻辑连接或会话**。错误处理策略包括不丢失消息、不重复消息、不篡改内容、不打乱多条消息顺序。连接管理还通过流量控制和拥塞控制处理消息传输延迟。

**面向连接服务技术**：**连接建立和释放**，序列号与确认，流量控制和拥塞控制。讲解步骤是先理论思考，再技术实现，最后看 TCP/UDP 的具体实现。

**连接建立**：**通信开始时需要建立逻辑连接**。问题包括如何告诉通信伙伴想建立新连接，以及接收方如何区分新的连接请求和旧的连接请求。解决思路是**握手协议**。

握手协议准备：为连接管理定义特殊消息类型，并把它们作为 PCI 控制信息的一部分。SYN 表示请求建立连接或同步；CLS 表示请求关闭连接；ACK 表示正确认；NAK 表示否认。这些都是协议描述的一部分。

**二次握手的问题**：旧请求或延迟消息可能导致新的建立请求被误解或忽略，因为接收方可能认为第一次尝试仍然有效。

**二次握手问题总结**：**通信错误有多种可能阻止连接建立**。练习是找另一个二次握手无法成功建立连接的例子。核心问题是如何解决这种歧义。

**握手协议扩展**：为了可靠建立连接，需要让双方都确认对方的同步信息，从而减少旧消息和延迟消息造成的混淆。

**三次握手**：3-Way Handshake 用三条控制消息完成连接建立，使双方都知道连接建立请求和确认是当前有效的。

**连接释放**：通信结束时，**也需要让通信伙伴知道逻辑连接要被关闭**。

**握手协议扩展到连接释放**：连接释放同样使用控制消息，例如 CLS 和 ACK，确保双方都知道连接方向是否关闭。

## 三、可靠通信：序列号与确认

**可靠通信回顾**：不可靠网络中可能发生**丢失、重复、损坏、乱序和延迟**，传输层需要机制来处理这些问题。

**序列号**：用于消息传输的计数器。每条消息带一个序列号，使接收方能够识别顺序、重复和缺失。

**确认**：接收方对收到的消息进行确认，发送方据此知道哪些消息已成功到达。

序列号记号：`PCI[N] = {i→j, SeqNr, AckNr, ...}` 表示控制信息中包含连接方向、序列号、确认号等。

**序列号示例**：消息依次编号，接收方根据序列号判断下一条期望消息。

序列号示例继续：**如果发生丢失、重复或乱序，序列号和确认号帮助检测并修复。**

序列号 Trace：通过一组消息交换展示 `SeqNr` 和 `AckNr` 如何随发送和确认更新。

**提高性能的优化**：**示例中每个字符都单独发送并确认，效率低。改进包括在 UD 中放更多数据，同时发送多个 PDU 并等待确认，使用累计确认或选择性确认。**

**进一步问题**：序列号到底有多少个，是否有限？这涉及序列号空间大小。消息在网络中能存活多久，也涉及消息最大生命周期。

- **序列号空间大小一**：序列号空间是有限的，通常由协议固定。若只有 10 个序列号，消息 11 会重复使用 0。若旧的延迟消息 0 仍在网络中，可能在新消息 11 前到达，接收方会把旧消息误认为新消息并丢弃真正的新消息。

- **序列号空间大小二**：需要知道消息最大延迟时间，并选择足够大的序列号空间，使最大延迟期间不会重复使用同一序列号。一种估计是序列号空间大于发送速率乘以 RTD。例子中 RTD 为 15 ms、发送速率 200/ms，需要 12 bit 编码序列号。

**消息最大生命周期**：可给所有消息**加时间戳和寿命**，**过旧消息直接丢弃**。这需要复杂估计或同步时钟，关联 Lamport 关于时间和事件顺序的经典论文。另一种方式是计数 hop。

**面向连接服务回顾**：逻辑连接、无丢失、无重复、顺序保证等已通过相关机制部分解决，仍需考虑内容损坏、流控和拥塞控制。

**总结**：面向连接协议常用序列号和确认实现。已发送消息被编号；丢失消息重新请求；重复消息被丢弃；乱序消息重新排序。

## 四、使用传输服务：端口、Socket、客户端服务器

使用传输服务：**编程基础**包括通信端点寻址、接口描述和编程模型。

**通信端点寻址**：传输层从端系统视角提供虚拟可靠端到端连接。问题是接收端 B 上是否有多个服务；如果有，如何选择正确服务。

图片过渡页，引出“**端口**”概念。

端口用于寻址通信端点。端口也称 **TSAP**，即 **Transport Service Access Point**。它是在给定主机上，传输层协议实例的通信端点地址。**目标是识别服务或协议，并区分同一端点对之间的多条连接**。
**端口特性**：端口空间是 **16 bit** 的扁平地址空间，**范围 0 到 65535**，由 IANA 管理。Well Known Ports 为 0-1023，Registered Ports 为 1024-49151，Dynamic/Private Ports 为 49152-65535。例子包括 **21 FTP、22 SSH、23 telnet、25 SMTP、80 HTTP、143 IMAP。**

传输服务接口：Socket 是平台无关、标准化的 API，用于访问网络全局通信端点和传输服务。连接发生在两个 Socket 之间。

Socket 图示：TCP 和 UDP 都有对应的 Socket 使用方式。

Socket 特性继续：**Socket 是双向的，可发送也可接收**。Internet Socket 包括 Datagram Socket 和 Stream Socket。连接端点是 IP 地址和端口的元组，连接总是在两个 Socket 之间。

**Client-Server 范式**：Server Socket 等待客户端连接请求；普通 Socket 在连接建立后由客户端和服务器进程用于通信。Server/Daemon 通常在建立连接后创建新线程或进程处理通信。

使用传输服务总结：端口用于端点寻址，**Socket 用于接口编程，Client-Server 是常见编程模型**。

发送端数据流回顾：连接信息在传输层体现为 Socket。

总结：端口用于通信端点寻址，Well Known/Registered 用于固定服务，Dynamic/Private 可自由使用。Socket 是编程接口，像文件一样 open/read/write/close。Client-Server 是传输服务的常见编程模型。

## 五、UDP 与 TCP 基础

传输层协议：无连接服务对应 UDP，特点是不可靠但快；面向连接服务对应 TCP，特点是可靠但慢。

UDP，User Datagram Protocol。本节介绍用户数据报协议。

**UDP 基础**：UDP 是无连接、不可靠的 Internet 传输协议，**见 RFC 768**。它没有连接建立、释放或状态；没有序列号，因此不保证顺序，也不检测重复或丢失；没有流量管理，也不调节发送速率。

**UDP 头部 PCI**：**源端口和目的端口各 16 bit，用于寻址。Datagram Length 为 16 bit**，IPv4 中最大长度为 65515 byte。**Checksum 为 16 bit**，是对 IP 伪头部和 UDP 数据报的可选校验和。

**UDP 特性**：头部为 64 bit，即 8 byte，比 TCP 至少 20 byte 更小。UDP 端口支持应用复用，UDP Socket 用于寻址应用。使用 UDP 的协议包括 TFTP、DNS、RPC、SNMP。

**UDP 使用**：**通过 UDP Socket 编程**，发送和接收独立数据报。通常每个数据报发送时都要指定接收方 IP 地址和端口。**顺序保证、丢失和重复处理、发送速率控制都由应用程序负责。**

TCP，Transmission Control Protocol。本节介绍传输控制协议。

**TCP 基础**：TCP 是**面向连接**、**可靠**的 Internet 传输协议，见 RFC 793。它支持端到端传输连接、全双工通信、丢失和重复等错误处理、流量控制和拥塞控制。连接建立使用三次握手。TCP 面向字节流，序列号和确认号针对 byte。

**TCP 头部字段二**：**源端口和目的端口各 16 bit，表示连接端点。序列号 32 bit**，是该段第一个用户数据 byte 的编号。**确认号 32 bit**，是下一期望用户数据 byte 的编号，只确认连续无错接收的数据。Header Length 为 4 bit，表示头部中 32 bit 字的数量。

**TCP 头部字段三**：NB 6 bit 未使用。Flags 包括 URG、ACK、PSH、RST、SYN、FIN。窗口大小 16 bit，表示从最后确认位置起还允许发送的 byte 数。

**TCP 头部字段四**：Checksum 是伪头部和 TCP 段所有 16 bit 字之和的一补码。Urgent Pointer 指向紧急数据之后的第一个 byte。Options 为可选属性，如最大段大小和 Timestamp。

**用 TCP 实现握手协议**：TCP Flags 可表达前面抽象控制消息。SYN 用于建立，ACK 用于确认，RST 可拒绝或重置，FIN 用于关闭；重复旧确认可视作确认重复。

**带序列号的连接建立/释放**：问题是 TCP 如何用序列号实现三次握手，尤其是在新连接开始时如何告诉通信伙伴初始序列号。接收方还必须区分延迟的旧连接请求和当前请求。

**连接建立**：三次握手。初始序列号由计时器选择，以便频繁建立/释放连接时仍保持唯一性。SYN 和 FIN 会消耗一个序列号，以保护连接管理。

**带序列号的三次握手**：第一条 SYN 包含 Pi 到 Pj 的初始序列号；SYN+ACK 确认该序列号并提供反向初始序列号；ACK 确认反向序列号，也可能携带第一批数据。延迟旧 SYN 可通过序列号识别。

**序列号记号**：`snn[i→j]` 是下一个要发送 byte 的序列号，`sne[i→j]` 是下一个期望 byte 的序列号。PCI 中的 `SeqNr` 是消息第一个 byte 的序列号，`AckNr` 是下一期望消息的序列号。

**三次握手后的数据传输**：**SYN 和 FIN 各使序列号加 1**；带用户数据的段按每个 byte 增加序列号；空段如纯确认则使用旧序列号。

**双向连接释放**：TCP 连接两个方向可分别关闭。一个方向发送 FIN 并得到 ACK 后，该方向关闭；另一方向仍可继续发送数据，直到也发送 FIN 并得到 ACK。

**TCP 特性**：TCP-PDU 称为 Segment。TCP 端口允许应用复用。TCP Socket 用于应用寻址。服务数据会被缓存。使用 TCP 的协议包括 HTTP、SMTP、IMAP、FTP。

TCP Socket 程序：客户端创建 Socket 并连接到主机 B 的端口 p；服务器创建 ServerSocket，等待客户端连接，accept 后用连接 Socket 读写数据，最后关闭连接。

客户端典型 TCP 状态。本页为状态图。

服务器典型 TCP 状态。本页为状态图。

**连接 Socket 示例**：服务器在 well-known port 上用 Server Socket 接收连接请求；连接建立后形成具体的 Socket 对，例如浏览器本地临时端口到 Web 服务器 80 端口。

Socket 使用：端到端连接绑定到服务或应用。Socket 是通信端点，由 IP 地址和端口号组成的端点元组表示，并由某个进程创建或拥有。连接由一对 Socket 给出。

`**netstat**` 示例：显示服务器监听的 TCP/UDP 端口和非服务器连接。**LISTEN** 表示等待连接，**ESTABLISHED** 表示已建立连接，**TIME_WAIT** 表示关闭后的等待状态。

**应用复用**：多个进程通过不同端口和 Socket 共享同一端系统的网络连接。端口和 Socket 使不同应用的数据能被正确分派。

**总结**：**UDP 和 TCP 都是传输层协议，提供端到端通信。UDP 无连接、不可靠但快；TCP 面向连接、可靠但慢。二者都通过 Socket 编程**，**端口允许应用和服务复用**。

**展望**：**UDP 性能高但不可靠，TCP 可靠但性能开销较大。可以用 UDP 在应用中重新实现部分 TCP 功能，若不需要 TCP 的所有特性，可降低开销**。

## 六、提高性能：窗口、流水线与可靠性


**改进目标**：提高通信性能，包括更高吞吐和更好利用率；提高可靠性，包括流量控制和拥塞/过载控制。目标是在更低错误率下获得更高性能。

序列号示例回顾：作为后续窗口技术的基础，消息编号和确认让发送方知道哪些数据已被接收。

回顾 Stop-and-Wait：发送方每发一个数据单元就等待确认，简单但吞吐受 RTD 限制。

**Stop-and-Wait 可达吞吐示例一**：由于一次只能有一个未确认分组在网络中，链路利用率可能很低，尤其在 RTD 很大时。

**Stop-and-Wait 可达吞吐示例二**：即使链路带宽很高，等待确认也会造成长时间空闲。

**Stop-and-Wait 可达吞吐示例三**：吞吐受分组大小和往返时延影响，性能通常低于链路理论带宽。

**改进方法**：不用 Stop-and-Wait，而是一次发送多个 PDU，让发送和确认重叠进行，从而提高链路利用率。

**滑动窗口协议/流水线**：Tanenbaum 称为 Schiebefensterprotokoll，Kurose/Ross 称为 Pipelining。**本质是在等待确认的同时继续发送一定窗口内的数据**。

**窗口技术动画**：发送窗口中包含允许发送但尚未全部确认的数据，确认到达后窗口向前滑动。

消息示例：“The quick brown fox jumps over the lazy frog.” 用于说明按字节或片段编号和窗口传输。

检查序列号：问题是接收方如何判断某个 SeqNr 是新消息还是重复消息。

新消息与重复消息范围图：序列号空间有限时，必须选择足够大的序列号空间和合适窗口大小，否则新旧数据会混淆。

窗口技术优点：它是通用机制，不只适用于本课程中的网络与分布式系统，也广泛用于其他通信和同步场景。

**Kurose/Ross 参考页**：介绍流水线、Go-Back-N 等可靠数据传输机制。

**Selective Repea**t，选择重传：接收方可缓存乱序到达的数据，发送方只重传丢失或出错的数据，而不是重传整个窗口。

**Kurose/Ross** 参考页：继续介绍选择重传和窗口协议细节。

**总结**：序列号是连续编号，确认用于反馈接收进度。**窗口技术允许多个未确认分组同时在网络中，提高吞吐；但窗口大小、序列号空间和最大生命周期必须匹配**。

## 七、流量控制与拥塞控制

**流量控制与拥塞控制区别**：流量控制是发送方防止接收方过载的措施；拥塞控制是防止多个独立通信关系共同使用中转网络时压垮网络的措施。

**队列图示**：发送方、网络节点、队列和中转网络共同决定是否出现等待和过载。

**流量控制方案**：Stop-and-Wait 会振荡，大 RTD 下很差；固定或动态窗口会遇到窗口缩小问题；发送前请求缓冲区预留会增加消息和延迟；接收方通过分配消息提供 Credit，如 TCP，但不够稳健；也可按时间栅格发送。

**流量控制一**：接收方缓冲区示例。发送方写入 2K，再写入 2K，4K 缓冲区变满；当应用读取 2K 后，接收方才能继续接收更多数据。

**流量控制二**：发送方变量包括 RcvWindow、LastByteSent、LastByteAcked；接收方变量包括 RcvBuffer、LastByteRcvd、LastByteRead。`RcvWindow = RcvBuffer - (LastByteRcvd - LastByteRead)`，且 `LastByteSent - LastByteAcked <= RcvWindow`。可对多个段使用累计确认。

**流量控制与拥塞控制回顾：流量控制解决接收方过载；拥塞控制解决中转网络过载。**

**网络过载**：**太多端系统上的太多传输实例同时向网络注入太多分组。**

**网络过载结果**：网络层路由器过载。传输层需要判断中转网络是否过载，并决定如何响应。

**发现过载**：一种方式是使用专门协议，让过载路由器显式通知端点，例如 XCP。另一种方式是通过间接指标判断，例如 Segment 丢失。

**拥塞控制方案**：虚拟信道安全但浪费缓冲并需要建立时间；限制进入网络的流量，如限制每个主机/进程连接数、限制发送速率；确保网络流出，如接收方保证接收率或在主机重组；保持网络恒定负载；按连接调节队列填充；Internet 中通常在过载时丢弃 Segment。

**对指标的反应**：丢失 Segment，如 Timeout 或 ACK 重复，意味着网络拥塞，应降低发送速率。收到 ACK 表示数据流动良好，可提高发送速率。TCP 通过持续试探提高速率，直到丢失发生，再降低速率。

**TCP 拥塞控制三步**：Slow Start 慢启动，指数增长；Congestion Avoidance 拥塞避免，线性增长；Fast Recovery 快速恢复。

**拥塞控制准备**：TCP 实例需要 CongWindow 控制注入网络的数据量，满足 `LastByteSent - LastByteAcked <= CongWindow`。有效窗口是 `min(RcvWindow, CongWindow)`。MSS 是最大 Segment 数据量。Threshold 是 CongWindow 快速增长的上界。RT 是重传计时器，超时后重发未确认段。

**TCP Tahoe 无负载时发送速率一**：Slow Start 从小拥塞窗口开始，`CongWin = 1 MSS`。到 Threshold 前指数增长，每收到一个 ACK 增加 1 MSS，每轮约为 1 RTD。


TCP Tahoe 慢启动图示：引用 Kurose/Ross，展示拥塞窗口按轮次加倍增长。

**TCP Tahoe 无负载时发送速率二**：Congestion Avoidance 中，达到 Threshold 后线性增长，每轮 `CongWin = CongWin + 1 MSS`，直到 Segment 丢失或 Timeout。

**TCP Tahoe 网络过载时发送速率三**：Timeout 时降低发送速率，把 Threshold 设为 `CongWin/2`，并用 Slow Start 重新开始，`CongWin = 1 MSS`。

丢失优化 RFC 5681：Segment 丢失后，后续收到的 Segment 会用同一确认号确认，形成重复 ACK。重复 ACK 可在 Timeout 前提前指示丢失。可能优化包括 Fast Retransmit 和 Fast Recovery。


Fast Retransmit：发送方在收到确认重复时，提前重发后续 Segment，不必等到 Timeout。

Fast Retransmit 示例图：展示通过重复 ACK 提前发现丢失并重传。

Fast Recovery：快速重传后，CongWin 仍被尊重但不直接回到 Slow Start。从第三个重复 ACK 起重新计算 Threshold 和 CongWin；每收到一个确认重复，就发送一个新数据 Segment，形成虚拟 CongWin。这对应 TCP Reno。

虚拟 CongWin：每收到一个重复 ACK，就发送一个新数据 Segment。因为丢失分组未被确认，CongWin 不会自然向前移动；为了继续发送，每个重复 ACK 让 CongWin 暂时加 1。Fast Recovery 后 CongWin 恢复为原值。

TCP Reno 网络过载时发送速率：Threshold 仍设为 `CongWin/2`。如果是三次重复 ACK，则 `CongWin = new Threshold`，继续 Congestion Avoidance；如果是 Timeout，则与 Tahoe 一样设置 `CongWin = 1 MSS` 并 Slow Start。

拥塞控制流程图：展示 TCP 拥塞控制不同状态和事件之间的转换。

拥塞控制状态图：展示 Slow Start、Congestion Avoidance、Fast Recovery 等状态关系。

TCP 拥塞避免的实用原则：AIMD，即 Additive Increase Multiplicative Decrease。没有 Slow Start 时，线性增加；发生丢包时，CongWin 按因子 2 减半。

TCP 吞吐计算：给定当前 `CongWin = w`、丢包时窗口 `W`、RTD。假设 W 在连接期间恒定，则吞吐在 `W/(2 RTD)` 和 `W/RTD` 之间摆动，近似 TCP 吞吐为 `0.75 W / RTD`。

**复习题一**：

- **哪些参数影响序列号空间大小？** Welche Parameter beeinflussen die Größe des Sequenznummernraums?

	1. 网络时延大小，以确保在消息的最大往返延迟（RTD）期间，序列号不会回绕并被重复使用
	2. TCP头部中序列号字段的位数

- **如何确保序列号唯一？** Wie wird die Eindeutigkeit der Sequenznummern sichergestellt?

	通过以下机制确保唯一性：选择足够大的序列号空间（32位提供约43亿个序列号）；使用初始序列号（ISN）随机化，每次连接使用不同的起始值；结合MSL确保旧报文段在序列号回绕前已过期；现代系统还使用时间戳选项（PAWS）来区分新旧报文段。

	Die Eindeutigkeit wird durch folgende Mechanismen gewährleistet: Verwendung eines ausreichend großen Sequenznummernraums (32 Bit bieten etwa 4,3 Milliarden Sequenznummern); Randomisierung der initialen Sequenznummer (ISN), wobei jede Verbindung einen anderen Startwert verwendet; Kombination mit MSL, um sicherzustellen, dass alte Segmente vor dem Wraparound ablaufen; moderne Systeme verwenden zusätzlich die Timestamp-Option (PAWS), um alte von neuen Segmenten zu unterscheiden.

- **窗口技术用于什么？** Wofür wird die Fenstertechnik verwendet?

	窗口技术主要用于**流量控制**和**提高传输效率**。它允许发送方在收到确认前**发送多个报文段**，从而**充分利用网络带宽**，避免"**停止等待**"协议的**低效问题**。同时，接收方通过通告窗口大小来防止缓冲区溢出，**实现端到端的流量控制**。

	Die Fenstertechnik dient hauptsächlich der Flusskontrolle und der Verbesserung der Übertragungseffizienz. Sie ermöglicht dem Sender, mehrere Segmente zu senden, bevor eine Bestätigung empfangen wird, wodurch die Netzwerkbandbreite optimal genutzt und die Ineffizienz des Stop-and-Wait-Protokolls vermieden wird. Gleichzeitig verhindert der Empfänger durch die Ankündigung der Fenstergröße einen Pufferüberlauf und realisiert so eine Ende-zu-Ende-Flusskontrolle.

- **窗口太小有什么影响？** Was sind die Auswirkungen eines zu kleinen Fensters?

	窗口太小会严重影响传输效率。发送方需要频繁等待确认，无法充分利用可用带宽，导致吞吐量下降。特别是在高带宽延迟积（BDP）的网络中，小窗口会成为瓶颈，因为管道无法被填满，网络资源被浪费。

	Ein zu kleines Fenster beeinträchtigt die Übertragungseffizienz erheblich. Der Sender muss häufig auf Bestätigungen warten und kann die verfügbare Bandbreite nicht vollständig nutzen, was zu einem verringerten Durchsatz führt. Besonders in Netzwerken mit hohem Bandbreiten-Verzögerungs-Produkt (BDP) wird ein kleines Fenster zum Engpass, da die Pipeline nicht gefüllt werden kann und Netzwerkressourcen verschwendet werden.

- **序列号空间太小有什么影响？** Was sind die Auswirkungen eines zu kleinen Sequenznummernraums?

	序列号空间太小会导致序列号回绕问题：在高速网络中，序列号可能在旧报文段仍在网络中传输时就已回绕并被重新使用，导致接收方无法区分新旧数据，造成数据混淆或损坏。这也是为什么TCP引入时间戳和PAWS机制来应对现代高速网络的原因。

	Ein zu kleiner Sequenznummernraum führt zu Wraparound-Problemen: In Hochgeschwindigkeitsnetzwerken können Sequenznummern bereits umlaufen und wiederverwendet werden, während alte Segmente noch im Netzwerk unterwegs sind. Dies führt dazu, dass der Empfänger alte und neue Daten nicht unterscheiden kann, was zu Datenverwechslung oder -beschädigung führt. Aus diesem Grund hat TCP Timestamps und den PAWS-Mechanismus eingeführt, um modernen Hochgeschwindigkeitsnetzwerken gerecht zu werden.


**复习题二**：

- **为什么传输层安全连接建立需要三次握手？** Warum erfordert der Aufbau einer sicheren Transportschichtverbindung einen Drei-Wege-Handshake?

	三次握手确保双方都能确认对方的接收和发送能力。第一次握手：客户端证明能发送；第二次握手：服务器证明能接收和发送；第三次握手：客户端证明能接收，同时确认服务器的初始序列号。这防止了旧的重复连接请求被误接受，确保双方初始序列号同步。

	Der Drei-Wege-Handshake stellt sicher, dass beide Seiten die Empfangs- und Sendefähigkeit der Gegenseite bestätigen können. Erster Handshake: Der Client beweist seine Sendefähigkeit. Zweiter Handshake: Der Server beweist seine Empfangs- und Sendefähigkeit. Dritter Handshake: Der Client beweist seine Empfangsfähigkeit und bestätigt die initiale Sequenznummer des Servers. Dies verhindert, dass alte duplizierte Verbindungsanfragen fälschlicherweise akzeptiert werden, und gewährleistet die Synchronisation der initialen Sequenznummern beider Seiten.

- **什么条件下二次握手足够？** **Unter welchen Bedingungen reicht ein Zwei-Wege-Handshake aus?**

	当网络完全可靠（无丢包、无延迟、无重复）且不存在旧连接请求干扰时，二次握手理论上足够。在实践中，这适用于某些简单的、无连接状态的协议，或者通信双方通过其他机制已建立信任关系的场景。但在不可靠网络中，二次握手无法防止旧请求导致的半开连接问题。

	Wenn das Netzwerk vollständig zuverlässig ist (kein Paketverlust, keine Verzögerung, keine Duplikate) und keine alten Verbindungsanfragen stören, reicht ein Zwei-Wege-Handshake theoretisch aus. In der Praxis gilt dies für einfache, zustandslose Protokolle oder Szenarien, in denen beide Kommunikationspartner bereits durch andere Mechanismen eine Vertrauensbeziehung aufgebaut haben. In unzuverlässigen Netzwerken kann ein Zwei-Wege-Handshake jedoch das Problem halboffener Verbindungen durch alte Anfragen nicht verhindern.

- **为什么连接释放必须双向发生？** Warum muss die Verbindungsfreigabe bidirektional erfolgen?

	**TCP是全双工协议，两个方向的数据流相互独立**。一方关闭只表示它不再发送数据，但仍可接收。四次挥手**允许半关闭状态**：A发送FIN表示A完成发送，但B可能还有数据要发给A。只有双方都发送并确认FIN后，连接才完全释放，确保所有数据都被传输完毕。

	TCP ist ein Vollduplex-Protokoll, bei dem die Datenströme in beiden Richtungen unabhängig voneinander sind. Wenn eine Seite schließt, bedeutet dies nur, dass sie keine Daten mehr sendet, aber noch empfangen kann. Der Vier-Wege-Handshake ermöglicht einen halb geschlossenen Zustand: A sendet FIN, um anzuzeigen, dass A das Senden beendet hat, aber B kann noch Daten an A senden. Erst wenn beide Seiten FIN gesendet und bestätigt haben, wird die Verbindung vollständig freigegeben, um sicherzustellen, dass alle Daten übertragen wurden.

- **流量控制和拥塞控制有什么区别？** Was ist der Unterschied zwischen Flusskontrolle und Staukontrolle?

	流量控制是端到端机制，**防止发送方压垮接收方的缓冲区**，通过接收窗口实现，**保护的是接收端**。拥塞控制是网络层面机制，**防止发送方压垮网络本身**，通过拥塞窗口实现，**保护的是整个网络**。

	Flusskontrolle ist ein Ende-zu-Ende-Mechanismus, der verhindert, dass der Sender den Puffer des Empfängers überlastet. Sie wird durch das Empfangsfenster (rwnd) realisiert und schützt den Empfänger. Staukontrolle ist ein Mechanismus auf Netzwerkebene, der verhindert, dass der Sender das Netzwerk selbst überlastet. Sie wird durch das Staufenster (cwnd) realisiert und schützt das gesamte Netzwerk.

**复习题三**：

- **哪些应用需求会支持使用 UDP？** Welche Anwendungsanforderungen sprechen für die Verwendung von UDP?

	实时性要求高、能容忍少量丢包的应用适合UDP：如视频会议、在线游戏、VoIP语音通话、直播流媒体。此外，简单的请求-响应模式（如DNS查询）、广播/多播通信、以及需要自定义可靠性机制的应用也选择UDP。UDP开销小、无连接建立延迟，适合对延迟敏感但对完整性要求不严格的场景。

	Anwendungen mit hohen Echtzeitanforderungen, die einen geringen Paketverlust tolerieren können, eignen sich für UDP: z.B. Videokonferenzen, Online-Spiele, VoIP-Telefonie und Live-Streaming. Darüber hinaus wählen auch einfache Anfrage-Antwort-Muster (wie DNS-Abfragen), Broadcast-/Multicast-Kommunikation und Anwendungen, die eigene Zuverlässigkeitsmechanismen implementieren, UDP. UDP hat geringen Overhead und keine Verbindungsaufbauverzögerung, was es für latenzempfindliche Szenarien mit weniger strengen Integritätsanforderungen geeignet macht.

- **使用 UDP 的应用如何保证可靠数据传输？** Wie gewährleisten Anwendungen, die UDP verwenden, eine zuverlässige Datenübertragung?

	应用层自行实现可靠性机制：添加序列号检测丢包和乱序；实现确认和重传机制；使用校验和验证数据完整性；设置超时定时器触发重传。典型例子是QUIC协议，在UDP之上构建了完整的可靠传输、拥塞控制和加密功能。游戏应用可能只对关键数据实现可靠传输，非关键数据允许丢失。

	Die Anwendungsschicht implementiert selbst Zuverlässigkeitsmechanismen: Hinzufügen von Sequenznummern zur Erkennung von Paketverlust und Neuordnung; Implementierung von Bestätigungs- und Neuübertragungsmechanismen; Verwendung von Prüfsummen zur Überprüfung der Datenintegrität; Setzen von Timeout-Timern zum Auslösen von Neuübertragungen. Ein typisches Beispiel ist das QUIC-Protokoll, das auf UDP vollständige zuverlässige Übertragung, Staukontrolle und Verschlüsselung aufbaut. Spieleanwendungen implementieren möglicherweise nur für kritische Daten zuverlässige Übertragung und erlauben den Verlust unkritischer Daten.

- **Alice 发送序列号 90 和 110 的两个 TCP Segment，如何判断第一个 Segment 的数据长度和 Bob 的确认号？** Alice sendet zwei TCP-Segmente mit Sequenznummern 90 und 110. Wie bestimmt man die Datenlänge des ersten Segments und Bobs Bestätigungsnummer?

	第一个Segment的数据长度 = 第二个序列号 - 第一个序列号 = 110 - 90 = **20字节**。
	Bob收到第一个Segment后的确认号 = 下一个期望接收的字节序列号 = 90 + 20 = **110**。
	确认号表示"我已收到110之前的所有数据，下一个期望收到序列号110"。
	
	Datenlänge des ersten Segments = zweite Sequenznummer - erste Sequenznummer = 110 - 90 = **20 Bytes**.
	
	Bobs Bestätigungsnummer nach Empfang des ersten Segments = nächste erwartete Byte-Sequenznummer = 90 + 20 = **110**.
	
	Die Bestätigungsnummer bedeutet: „Ich habe alle Daten vor 110 empfangen und erwarte als nächstes Sequenznummer 110."

**复习题四**：

- **TCP 使用哪些机制实现可靠端到端连接？** Welche Mechanismen verwendet TCP für eine zuverlässige Ende-zu-Ende-Verbindung?

	TCP通过多种机制实现可靠性：序列号确保数据按序重组并检测丢失；确认机制（ACK）让发送方知道数据已到达；超时重传处理丢包；校验和检测数据损坏；流量控制（滑动窗口）防止接收方溢出；三次握手建立可靠连接；四次挥手确保完整释放。

	TCP erreicht Zuverlässigkeit durch mehrere Mechanismen: Sequenznummern stellen sicher, dass Daten in der richtigen Reihenfolge zusammengesetzt werden und Verluste erkannt werden; Bestätigungsmechanismen (ACK) informieren den Sender über angekommene Daten; Timeout-Neuübertragung behandelt Paketverluste; Prüfsummen erkennen Datenbeschädigungen; Flusskontrolle (Schiebefenster) verhindert Überlauf beim Empfänger; Drei-Wege-Handshake baut zuverlässige Verbindungen auf; Vier-Wege-Handshake gewährleistet vollständige Freigabe.

- **传输层如何实现应用复用？** Wie realisiert die Transportschicht Anwendungsmultiplexing?

	通过端口号实现复用和解复用。每个应用绑定唯一端口号，传输层用源端口+目标端口+源IP+目标IP的四元组标识每个连接。发送时，传输层将多个应用的数据复用到同一网络连接；接收时，根据端口号将数据解复用分发给对应应用进程。

	Multiplexing und Demultiplexing werden durch Portnummern realisiert. Jede Anwendung bindet sich an eine eindeutige Portnummer, und die Transportschicht identifiziert jede Verbindung durch das Vier-Tupel aus Quellport + Zielport + Quell-IP + Ziel-IP. Beim Senden multiplext die Transportschicht Daten mehrerer Anwendungen auf dieselbe Netzwerkverbindung; beim Empfangen werden die Daten anhand der Portnummer demultiplext und an den entsprechenden Anwendungsprozess verteilt.

- **TCP Tahoe 如何进行拥塞控制？** Wie führt TCP Tahoe die Staukontrolle durch?

	Tahoe包含三个阶段：慢启动阶段，cwnd从1开始指数增长（每RTT翻倍），直到达到ssthresh；拥塞避免阶段，cwnd线性增长（每RTT加1）；检测到丢包时（超时或三次重复ACK），ssthresh设为当前cwnd的一半，cwnd重置为1，重新进入慢启动。Tahoe对所有丢包事件反应相同，都重置cwnd为1。

	Tahoe umfasst drei Phasen: In der Slow-Start-Phase wächst cwnd exponentiell ab 1 (Verdopplung pro RTT) bis ssthresh erreicht wird; in der Stauvermeidungsphase wächst cwnd linear (plus 1 pro RTT); bei erkanntem Paketverlust (Timeout oder drei doppelte ACKs) wird ssthresh auf die Hälfte des aktuellen cwnd gesetzt und cwnd auf 1 zurückgesetzt, dann beginnt erneut Slow-Start. Tahoe reagiert auf alle Paketverlustereignisse gleich und setzt cwnd immer auf 1 zurück.

- **窗口技术中发送窗口和接收窗口能否大小不同？** Können Sendefenster und Empfangsfenster bei der Fenstertechnik unterschiedlich groß sein?

	可以不同，而且通常确实不同。接收窗口（rwnd）由接收方根据其缓冲区空间通告；发送窗口受rwnd和拥塞窗口（cwnd）共同限制，取两者最小值。接收方缓冲区大小、网络拥塞状况、以及双方处理能力不同都会导致窗口大小差异。发送方必须尊重接收方通告的窗口大小。

	Ja, sie können unterschiedlich sein und sind es normalerweise auch. Das Empfangsfenster (rwnd) wird vom Empfänger basierend auf seinem Pufferspeicher angekündigt; das Sendefenster wird sowohl durch rwnd als auch durch das Staufenster (cwnd) begrenzt und nimmt das Minimum beider Werte. Unterschiedliche Puffergrößen beim Empfänger, Netzwerkstausituationen und unterschiedliche Verarbeitungskapazitäten beider Seiten führen zu unterschiedlichen Fenstergrößen. Der Sender muss die vom Empfänger angekündigte Fenstergröße respektieren.


## 专用词汇表

| 德语/英文 | 中文 | 说明 |
|---|---|---|
| Transport | 传输 | 第三章主题，关注端到端通信 |
| Transportschicht | 传输层 | ISO/OSI 第 4 层 |
| Ende-zu-Ende | 端到端 | 端系统之间的通信视角 |
| Best Effort | 尽力而为 | 网络层常见服务模型，不保证可靠性 |
| Zuverlässigkeit | 可靠性 | 不丢失、不重复、不乱序等保证 |
| Verbindungslos | 无连接 | 不建立连接状态的通信方式 |
| Verbindungsorientiert | 面向连接 | 建立逻辑连接/会话的通信方式 |
| Handshake | 握手 | 建立或释放连接的控制消息交换 |
| 2-Way Handshake | 二次握手 | 请求与确认两步建立 |
| 3-Way Handshake | 三次握手 | TCP 连接建立基本机制 |
| SYN | 同步/建立请求 | TCP 建连标志 |
| ACK | 确认 | 确认收到或确认序列号 |
| NAK | 否认 | 表示未正确收到 |
| FIN / CLS | 关闭请求 | 连接释放相关控制信息 |
| Sequenznummer / SeqNr | 序列号 | 标识数据顺序 |
| Quittung / AckNr | 确认号 | 表示下一期望数据 |
| RTD | 往返时延 | Round Trip Delay |
| Maximale Lebensdauer | 最大生命周期 | 消息在网络中可能存在的最长时间 |
| Port | 端口 | 传输层端点地址 |
| TSAP | 传输服务接入点 | Transport Service Access Point |
| Well Known Port | 知名端口 | 0-1023 |
| Registered Port | 注册端口 | 1024-49151 |
| Dynamic/Private Port | 动态/私有端口 | 49152-65535 |
| Socket | 套接字 | 程序访问传输服务的 API |
| Datagram Socket | 数据报 Socket | 常用于 UDP |
| Stream Socket | 流 Socket | 常用于 TCP |
| Client-Server | 客户端-服务器 | 常见网络编程模型 |
| UDP | 用户数据报协议 | 无连接、不可靠、低开销 |
| TCP | 传输控制协议 | 面向连接、可靠、字节流 |
| Datagramm | 数据报 | UDP 传输单位 |
| Segment | 段 | TCP PDU |
| Header | 头部 | 协议控制信息 |
| Checksum | 校验和 | 检测错误 |
| Pseudoheader | 伪头部 | TCP/UDP 校验中使用的 IP 字段集合 |
| Multiplexing | 复用 | 多应用共享传输层服务 |
| Stop-and-Wait | 停等 | 发一个等一个确认 |
| Sliding Window | 滑动窗口 | 允许多个未确认数据在途 |
| Pipelining | 流水线 | 发送与确认重叠进行 |
| Selective Repeat | 选择重传 | 只重传丢失/出错数据 |
| Flow Control | 流量控制 | 防止接收方过载 |
| Congestion Control | 拥塞控制 | 防止网络中转部分过载 |
| RcvWindow | 接收窗口 | 接收方可用缓存对应窗口 |
| CongWindow / CongWin | 拥塞窗口 | 控制发送方注入网络的数据量 |
| MSS | 最大段大小 | Maximum Segment Size |
| Threshold | 阈值 | 慢启动转拥塞避免的边界 |
| Timeout | 超时 | 计时器到期触发重传 |
| Slow Start | 慢启动 | 拥塞窗口指数增长 |
| Congestion Avoidance | 拥塞避免 | 拥塞窗口线性增长 |
| Fast Retransmit | 快速重传 | 重复 ACK 触发提前重传 |
| Fast Recovery | 快速恢复 | TCP Reno 的丢包恢复策略 |
| TCP Tahoe | TCP Tahoe | 慢启动、拥塞避免、超时后重启 |
| TCP Reno | TCP Reno | 加入快速重传和快速恢复 |
| AIMD | 加性增、乘性减 | 拥塞避免基本原则 |
| QUIC | QUIC 协议 | 基于 UDP 的多路复用安全传输 |

## 讲义常见句式

| 原句式 | 中文译法 | 讲义中的用法 |
|---|---|---|
| `Dienste der ... erwarten` | ……层的服务期望 | 传输层需求说明 |
| `Dienste der ... bieten` | ……层的服务提供 | 网络层能力说明 |
| `mit oder ohne ...` | 带或不带…… | 可靠性保证 |
| `ohne Zuverlässigkeitsgarantien` | 不提供可靠性保证 | Best Effort |
| `Mögliche Probleme bei ...` | ……中可能的问题 | 消息投递问题 |
| `Zielkonflikt` | 目标冲突 | 可靠性与性能 |
| `Trade-Off erforderlich` | 需要权衡 | 成本收益权衡 |
| `Keine Verbindung zwischen ...` | ……之间没有连接 | 无连接服务 |
| `Aufbau einer logischen Verbindung` | 建立逻辑连接 | 面向连接服务 |
| `Wie wird ... mitgeteilt?` | 如何通知……？ | 连接建立问题 |
| `Teil der Protokollbeschreibung` | 协议描述的一部分 | 控制消息类型 |
| `Zur Erinnerung` | 回顾 | TCP、序列号等复习页 |
| `im Speziellen` | 特别是 | TCP 握手细节 |
| `bezieht sich auf ...` | 针对/涉及…… | TCP 字节流序列号 |
| `wird verworfen` | 被丢弃 | 旧消息、过载分组 |
| `wird neu gesendet` | 被重新发送 | 超时或快速重传 |
| `je ACK` | 每收到一个 ACK | CongWin 增长 |
| `bei Timeout` | 超时时 | TCP Tahoe/Reno 行为 |
| `identisch zu ...` | 与……相同 | Reno 与 Tahoe 比较 |
| `gegeben durch ...` | 由……给出 | Socket 对、连接 |
| `obliegen ...` | 由……负责 | UDP 应用层责任 |
| `Zusammenfassung` | 总结 | 小节总结 |
| `Fragen` | 问题/复习题 | 章节末复习 |

