# RaMAx 项目开发任务清单

## 项目目标

RaMAx 旨在实现基于聚类的局部 all-all 渐进式多基因组比对，对标 [progressive cactus](https://github.com/ComparativeGenomicsToolkit/cactus)。本阶段目标为完成基础代码框架，并确保每项功能有独立的子程序实现。主要算法流程见下：

---

### 算法流程

#### 预处理阶段

1. 读取 seqfile（与 cactus 格式一致），获取进化树和所有基因组名称及路径（支持超链接）。
2. 解析进化树，将所有基因组名称与路径存入字典，并把本地/云端文件同步到工作目录。
3. 清洗序列，仅保留 AGCTN，执行[重复序列掩蔽](https://github.com/BioinformaticsToolsmith/Red)，并生成对应的 fai 索引。
4. 基于进化树分组，每组约 5 个基因组，构建比对指导树，并对查询序列预切片（掩蔽序列不参与查询）。

#### 比对阶段

1. 对组内找到中心星序列构建索引，把其他序列比对到中心星序列行
2. 过滤比对结果，得到初步结果。对于剩下基因组的剩下部分，继续完成星比对，依次迭代
3. 最后输出为hal文件

#### 合并阶段

1. 按指导树完成全部比对，生成所有 hal 文件。
2. 按指导树顺序合并 hal 文件。

---

## 当前优先任务（🔴🟠🟢）

### 预处理阶段

### 索引阶段
- [ ] 🟢 优化FM_index类中的bisectAnchors和findAnchorsAccurate函数，因为现在找到的锚点不是严格意义的MUM，还可以向左扩展，并确保输出准确


### 比对阶段
- [ ] 🔴 目前哺乳类动物的比对效果很差，尝试进行优化
- [ ] 🟢 增加对重复序列的比对支持
- [ ] 🟢 尝试引入MUMmer的后缀数组＋LCP索引方案
- [ ] 🟢 使用slaping来加速后缀数组查找

### 过滤阶段
- [ ] 🔴 目前灵长类动物比对的相关指标很差，需要优化比对的过滤算法
- [ ] 🟠 完成Anchor两端的延申
- [ ] 🟢 将repeat Match补充到比对结果中

### 构图阶段
- [ ] 🔴 基于比对图构建公共祖先
- [ ] 🟢 segment存在0的长度
- [ ] 🟠 对最终图进行润色优化
- [ ] 🟢 串行构图优化为并行构图

### 结果输出
- [ ] 🔴 多基因组比对完成hal文件的输出
- [ ] 🟢 双基因组比对支持sam，maf，paf，delta，lav，axt格式的输出
- [ ] 🟢 实现自定义比对格式，开发自定义格式与其他格式的转换工具

### 整体框架
- [ ] 🔴 支持 group 模式（基于聚类的渐进式比对）和支持 star 模式（所有序列完成星比对）
- [ ] 🟢 实现预处理模块子程序 RaMA-preprocess
- [ ] 🟢 实现比对模块 RaMA-G（预留子命令 build、align、map）
- [ ] 🟢 实现合并模块 RaMA-merge
- [ ] 🟢 实现pipeline模块RaMA-prepare，输出所有命令完成比对
- [ ] 🟠 完善单元测试模块
- [ ] 🟢 增加conda和docker的安装方式
- [ ] 🟢 完善README文档，包含安装、使用、开发等说明
- [ ] 🟢 单独把fmidex作为一个轮子开源到一个新仓库
- [ ] 🟢 单独把SeqPro作为一个轮子开源到一个新仓库

### 优化相关
- [ ] 🔴 增加对contig, Scaffold基因组的支持，还没有测试
- [ ] 🟢 把能改为unique的shared指针都进行修改
- [ ] 🔴 优化内存泄露

### bug修复
- [ ] 🟠 query的相同位点存在两个相同位置的匹配

### 用户相关
- [ ] 🟢 开发可视化软件，支持自定义格式的可视化
- [ ] 🟢 为RaMAx做一个主页


## 已经完成的任务
- [x] 🔴 在data_process.h中对所有基因组完成[重复掩蔽](https://github.com/BioinformaticsToolsmith/Red)
- [x] 🔴 在rare_aligner.h中新增MultipeRareAligner类，用于多基因组比对，PairRareAligner暂时放弃开发。
- [x] 🟢 实现进化树读取与解析，构建比对指导树
- [X] 🟠 确保被掩蔽的序列不在参考序列中参与比对
- [x] 🟠 在index.h中完成fastamanger的优化，只有部分序列的fasta快速索引到原有位置
- [X] 🟠 支持反向链的寻找
- [x] 🔴 完成rare_aligner.h中starAlignGroup，确认后续的开发任务
- [x] 🔴 完成初步的星比对
- [X] 🟠 match的匹配不对应
- [x] 🟠 fasta_manager类中的charidxmap不支持重复掩蔽
- [x] 🔴 完成多基因组比对的main函数RaMAx.cpp
- [x] 🔴 完成benchmark文件夹的构建，包括双基因组比对和多基因组比对
- [X] 🟠 支持反向链的锚点过滤和聚类
- [X] 🟠 将fastamanager拆分为两个类，第一个类是fastaprocessor，用于处理fasta文件，第二个类是fastamanager，用于索引fasta文件
- [X] 🔴 在anchor.h中参考 MUMmer4 mgaps.cc 实现锚点聚类这个类的实现，并完成锚点过滤的功能。
- [X] 🟢 在index.h中完成fmindex的优化，支持拼接的序列用\01作为分割符号
- [x] 🔴 完成并行图的合并
   - [x] 🔴 并行查找所有重叠区域
      - [x] 🔴 高层并行，按染色体划分
      - [x] 🟠 中层并行，利用采样表构建局部区间树查询
   - [x] 🔴 并行分割与合并，为每一对重叠创建并行任务
      - [x] 🔴 分割，将seg分割成三部分，重叠前，重叠区，重叠后，需要根据CIGAR同步分割
      - [x] 🔴 合并，创建新的block，并将原先的block
- [X] 🟠 基于贪婪算法完成Match聚类的选择，构建RaMesh图
- [X] 🔴 完成RaMesh的图结构的设计和开发
- [X] 🔴 完成星比对结果的图结构合并
- [X] 🟢 完成迭代补充星比对的图结构结果
- [X] 🟠 Match在cluster的过程中，可能会把overlap的两个锚点其中一个删掉，但实际应该合并他们
- [X] 🟢 优化cmake(如sdsl的编译速度,windowmasker的安装（权限)，使用新的cmake来管理submodule）
- [X] 🟠 在align.h中完成双序列比对类的开发，支持锚点的延伸和细化
- [X] 🟠 findQueryFile函数线程池现在是初步使用片段数量切分chunk，可以优化为片段长度之和,把多个片段合并为一个序列
---

## 协作与开发规范

1. **命名规范**

   * 变量：下划线分隔
   * 函数：小驼峰（如 `getAnchorList`）
   * 类名：大驼峰（如 `FastaManager`）
2. **分支协作**

   * 每人维护一条开发分支，独立开发，功能完成后通过 PR 合并至主分支，合并后及时 pull 更新
3. **代码质量**

   * 合并前确认无明显 bug，重要模块补充英文注释
