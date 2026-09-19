# 第一章 Motivation & Organisatorisches 中文翻译

来源文件：`lecture/Kapitel 1 Motivation Organisatorisches.pdf`

说明：本译稿按 PDF 页序逐页整理。第一章主要是课程动机与组织信息，图示页保留页号并翻译可抽取文字；日期、邮箱、链接等按讲义原文保留。

## 一、课程动机：为什么要理解 Internet

### 第 1 页
《计算机网络与分布式系统》课程封面。慕尼黑大学，2026 夏季学期，Prof. Dr. D. Kranzlmueller。课程主页：https://www.ifi.lmu.de/rnvs

### 第 2 页
第 1 章：动机与组织事项。

### 第 3 页
“慕尼黑的早晨”。本页以慕尼黑文化图片作为开场，包括白香肠、椒盐卷饼、芥末和啤酒等 Oktoberfest 元素。

### 第 4 页
乌克兰战争与 Internet：基于 RIPE Atlas 的可达性地图，时间为 2022 年 11 月。图例包括连接可用、无数据、连接不可用。

### 第 5 页
乌克兰战争与 Internet：同样是 RIPE Atlas 可达性地图，时间为 2023 年 3 月。用于对比战争期间网络连接状态的变化。

### 第 6 页
来自乌克兰的数据流量：展示自 2022 年 2 月 24 日以来，从乌克兰到 Cloudflare 的 Internet 流量变化，来源为 Cloudflare 关于乌克兰战争一周年的报告。

### 第 7 页
韧性与 Internet：Internet 的基本思想之一是容错，例如 TCP、BGP 等协议。Internet 由许多独立自治系统组成，具有灵活路径选择和冗余能力。

### 第 8 页
连接如何保持可达：答案是去中心化和 Peering。讲义用 2022 年 11 月与 2023 年 3 月的乌克兰 peer-to-peer 连接情况进行对比。

### 第 9 页
ASN 14593，Starlink：展示乌克兰到 Cloudflare 的 Internet 流量中，经由 ASN 14593 的比例。ASN 14593 对应 Starlink。

### 第 10 页
相关报道：本页列出关于乌克兰战争期间 Internet 韧性的阅读材料，包括 Cloudflare、RIPE Labs、The Economist、First Monday、Der Standard、Spiegel、Washington Post 等来源。

### 第 11 页
计算机科学学习：每一位计算机科学学习者都必须知道“Internet 是如何工作的”。

### 第 12 页
“Internet 对我们所有人来说都是 Neuland”。本页引用“Neuland”这一说法，并给出 YouTube 链接，用作引入或调侃：即使 Internet 已很普遍，人们对其底层机制仍常不熟悉。

### 第 13 页
Internet。本页为图示/过渡页，用于引出对 Internet 的实际理解。

### 第 14 页
The IT Crowd：课程引用英剧《The IT Crowd》第 3 季第 4 集 “Die Rede”，并列出 Channel 4 链接，作为关于 Internet 误解的幽默例子。

### 第 15 页
关于课程团队自身。本页引出后续关于 LRZ、MWN 和 MNM 团队的介绍。

### 第 16 页
空白/过渡页。

## 二、LRZ、MWN 与 MNM 团队

### 第 17 页
Leibniz-Rechenzentrum，巴伐利亚科学院的莱布尼茨计算中心。地点位于 Garching 研究中心园区。

### 第 18 页
SuperMUC-NG，Next Generation：德国国家级高性能计算机。本页展示超级计算机图片。

### 第 19 页
慕尼黑科研网 MWN。本页为 MWN 图示。

### 第 20 页
MWN 数据概览：MWN 是服务慕尼黑高校和科研机构的网络，覆盖约 136,000 名学生和 30,000 名员工。关键数据包括 14 个核心路由器，计划 21 个；64 个站点路由器；2,700 台交换机和 6,400 个接入点；83 条租用 Dark Fibre 和 40 多条私有 Dark Fibre；超过 200,000 台设备；90 个地点和 650 栋建筑。2024 年 3 月数据传输量为每月 6,300/2,700 TB 进/出，Backbone 每月约 70 PB。

### 第 21 页
负载情况：Garching 的 X-WiN 连接利用率。本页为流量/负载图。

### 第 22 页
MWN 用户：包括慕尼黑艺术学院、大学和高校建设部门、巴伐利亚科学院、巴伐利亚州立图书馆、巴伐利亚州立博物馆与收藏机构、剧院学院、植物园、德国心脏中心、多所高校、Fraunhofer、Max Planck、LMU、TUM、学生宿舍和多个研究/文化机构等。

### 第 23 页
MWN 结构。本页展示慕尼黑科研网的结构图。

### 第 24 页
MWN、WiN、GEANT：面向科学研究的网络。图中提到 DFN 的光纤平台，以及连接成本约 0.85 百万欧元/年。

### 第 25 页
Munich Network Management Team：慕尼黑网络管理团队，研究和管理网络与分布式系统。负责人 Prof. Dr. Dieter Kranzlmueller。

### 第 26 页
MNM 团队。本页为团队图片/过渡页。

### 第 27 页
MNM 团队成员：列出 Prof. Dr. Dieter Kranzlmueller、Prof. Dr. Heinz-Gerd Hegering，以及多位教授、博士、研究人员和团队成员。

### 第 28 页
选定研究主题：包括高性能和最高性能计算、并行化、算法与数据结构、软件工具、分析与调试、能源效率、可重构计算与加速器；分布式计算、Grid Computing、Cloud Computing、Edge/Fog Computing、数据保护与伦理；Future Computing、Exascale Computing、量子计算；网络监控、网络中立性、软件定义网络；物联网通信协议和标准化；虚拟化；IT 安全、后量子安全、安全管理和风险方法；IT 服务管理；VR 与可视化等。

### 第 29 页
相关课程：包括本课程“计算机网络与分布式系统”，重点是分层架构和 Internet Protocol Suite；Grid and Cloud Computing；IT-Management；IT-Sicherheit；Parallel Computing；计算机网络实践课；IT 安全实践课。

## 三、课程组织信息

### 第 30 页
上课时间：每周五 09:15-11:45，10:30-10:45 休息。答疑课预计 17.07.2026。考试在学期结束后不久，时长 90 分钟。可能会有客座讲者。

### 第 31 页
组织事项：关于课件的问题或提示发送至 `rnvs-skript@nm.ifi.lmu.de`。习题课负责人邮箱为 `rnvs@nm.ifi.lmu.de`，负责人包括 Daniel Diefenthaler 和 Fabian Dreer。更多课程信息在 Moodle 课程中：https://moodle.lmu.de/course/view.php?id=44443

### 第 32 页
Moodle 与录播：课程加入密钥为 `HotPotato`。Moodle 提供课件、公告和讨论论坛、习题负责人/Prof. Kranzlmueller/助教交流、习题表和提交入口。往年视频在 LMUCast：https://cast.itunes.uni-muenchen.de/vod/playlists/PBZ51y45zp.html

### 第 33 页
习题课安排：Tutorien 从 2026 年 4 月 20 日星期一开始。习题每周一发布。没有固定习题组分配。往年助教视频链接：https://cast.itunes.uni-muenchen.de/vod/playlists/dzIS6o2dJR.html

### 第 34 页
每周习题表：完成的习题可自愿通过 Moodle 提交；提交的习题由助教批改；标有 H 的题一定会讨论；其他题目会在明确、具体提问并展示自己解法思路时讨论。

### 第 35 页
习题课时间：所有时间为 c.t.。时间包括周一 10-12、12-14、14-16；周二 14-16、16-18、18-20；周三 10-12、14-16、16-18。地点为 Geschwister-Scholl-Pl. 1, D Z007，对应不同 Tutor。

### 第 36 页
成功参与课程：考试计划包括学期考试，预计 7 月底；补考预计 10 月初。考试内容包括全部讲课内容、全部习题内容以及习题课讨论过的所有主题。讲座和习题课构成一个整体。

### 第 37 页
答疑课：讲义原文写为预计 2025 年 7 月 17 日星期五，课程结束时举行。问题可在 2025 年 7 月 14 日前随时提交至 `rnvs-fragen@nm.ifi.lmu.de`。答疑课只讨论已提交的问题，不处理关于考试内容范围的问题。注意：本课程封面为 2026 夏季学期，此处年份可能是讲义中的遗留日期。

### 第 38 页
计算机网络文献：参见 https://www.nm.ifi.lmu.de/rn.html。本页为网络方向参考书目引入页。

### 第 39 页
计算机网络文献。本页为参考文献图示。

### 第 40 页
分布式系统文献。本页为参考文献图示。

## 专用词汇表

| 德语/英文 | 中文 | 说明 |
|---|---|---|
| Motivation | 动机 | 课程为什么重要 |
| Organisatorisches | 组织事项 | 课程安排、Moodle、习题、考试 |
| Rechnernetze | 计算机网络 | 本课程主题之一 |
| Verteilte Systeme | 分布式系统 | 本课程主题之一 |
| Internet | 互联网 | 本章动机核心对象 |
| Resilienz | 韧性 | 系统在故障、攻击、战争等情况下保持运行的能力 |
| Fehlertoleranz | 容错 | Internet 设计基本思想之一 |
| TCP | 传输控制协议 | 本页作为容错相关协议例子 |
| BGP | 边界网关协议 | 自治系统间路由协议，本页作为韧性例子 |
| Autonomes System | 自治系统 | 独立管理的网络系统 |
| Wegewahl | 路径选择/路由 | Internet 能灵活选择路径 |
| Redundanz | 冗余 | 通过多条路径或备用资源提升韧性 |
| Dezentralisierung | 去中心化 | 没有单一中心，提升可用性 |
| Peering | 对等互联 | 网络之间直接交换流量 |
| ASN | 自治系统编号 | Autonomous System Number |
| Starlink | 星链 | ASN 14593，对乌克兰网络可达性案例有影响 |
| Datenverkehr | 数据流量 | 网络中传输的数据量 |
| Cloudflare | Cloudflare | 讲义中用于观察乌克兰流量的网络服务商 |
| RIPE Atlas | RIPE Atlas | 用于测量网络可达性的全球探针平台 |
| Leibniz-Rechenzentrum / LRZ | 莱布尼茨计算中心 | 巴伐利亚科学院计算中心 |
| SuperMUC-NG | SuperMUC 下一代 | 国家级高性能计算机 |
| MWN | 慕尼黑科研网 | Muenchner Wissenschaftsnetz |
| X-WiN | 德国科研骨干网 | DFN 运营的科学网络连接 |
| GEANT | 欧洲科研教育网络 | 欧洲范围研究网络 |
| DFN | 德国科研网协会 | Deutsches Forschungsnetz |
| Dark Fibre | 暗光纤 | 已铺设但由用户自行点亮/运营的光纤 |
| Backbone | 骨干网 | 网络核心传输部分 |
| Access Point | 接入点 | 无线网络接入设备 |
| Core-Router | 核心路由器 | 骨干网络中的核心路由设备 |
| Standort-Router | 站点路由器 | 连接具体地点/校区的路由器 |
| MNM-Team | 慕尼黑网络管理团队 | Munich Network Management Team |
| Hochleistungsrechnen | 高性能计算 | HPC 相关研究方向 |
| Hoechstleistungsrechnen | 最高性能计算 | 超算/极高性能计算 |
| Exascale Computing | 百亿亿次级计算 | 未来高性能计算方向 |
| Quantencomputing | 量子计算 | Future Computing 主题 |
| Software-Defined Networks | 软件定义网络 | 网络研究主题 |
| Internet-of-Things | 物联网 | 网络与分布式系统应用领域 |
| Virtualisierung | 虚拟化 | 管理虚拟基础设施 |
| IT-Sicherheit | IT 安全 | 课程相关领域 |
| Post-Quanten-Sicherheit | 后量子安全 | 抗量子计算攻击的安全技术 |
| IT-Servicemanagement | IT 服务管理 | 运营联网系统的管理方法 |
| Vorlesung | 讲座/课程 | 正课 |
| Tutorium | 习题课/辅导课 | 助教带领的练习课 |
| Uebungsblatt | 习题表 | 每周练习 |
| Klausur | 考试 | 课程考核 |
| Nachholpruefung | 补考 | 第二次考试机会 |
| Fragestunde | 答疑课 | 课程末问题讨论 |
| Moodle | Moodle 平台 | 课程资料和提交系统 |
| Einschreibeschluessel | 加课密钥 | Moodle 课程注册密码 |
| c.t. | 延后 15 分钟开始 | 德语大学时间习惯，cum tempore |

## 讲义常见句式

| 原句式 | 中文译法 | 讲义中的用法 |
|---|---|---|
| `Kapitel 1: ...` | 第 1 章：…… | 章节标题 |
| `Krieg in der Ukraine und das Internet` | 乌克兰战争与 Internet | 动机案例标题 |
| `Verbindung moeglich` | 连接可用 | RIPE Atlas 图例 |
| `Keine Daten` | 无数据 | RIPE Atlas 图例 |
| `Keine Verbindung moeglich` | 连接不可用 | RIPE Atlas 图例 |
| `Fehlertoleranz als Grundgedanke` | 容错作为基本思想 | Internet 韧性说明 |
| `Wie bleiben ... erreichbar?` | ……如何保持可达？ | 可达性讨论 |
| `Anteil des Internetverkehrs ...` | ……的 Internet 流量占比 | Starlink 流量图 |
| `Berichte` | 报道/阅读材料 | 资料列表 |
| `Jede/r ... muss wissen ...` | 每一位……都必须知道…… | 课程动机 |
| `In eigener Sache` | 关于我们自己/课程团队 | 团队介绍引入 |
| `in Zahlen` | 用数字看…… | MWN 数据页 |
| `Netz fuer ...` | 面向……的网络 | MWN 服务对象 |
| `Ausgewaehlte Forschungsthemen` | 选定研究主题 | MNM 研究方向 |
| `Lehrveranstaltungen im Umfeld` | 相关课程 | 课程体系说明 |
| `Jeweils freitags ...` | 每周五…… | 上课时间 |
| `vsl.` | 预计/可能 | voraussichtlich 的缩写 |
| `kurz nach Semesterende` | 学期结束后不久 | 考试时间 |
| `Fragen/Hinweise zu ...` | 关于……的问题/提示 | 联系邮箱 |
| `Weitere Informationen ... gibt es ...` | 更多信息见…… | Moodle 说明 |
| `Beginn der Tutorien ...` | 习题课开始于…… | 习题安排 |
| `Freiwillige Abgabe` | 自愿提交 | 习题提交 |
| `in jedem Fall besprochen` | 一定会讨论 | H 标记习题 |
| `auf explizite/konkrete Nachfrage` | 在明确/具体提问时 | 非 H 题讨论条件 |
| `Vorlesung und Uebung bilden eine Einheit` | 讲座和习题课构成整体 | 考试范围提醒 |
| `Siehe auch` | 另见/参见 | 参考资料页 |

