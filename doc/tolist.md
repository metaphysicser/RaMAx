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

1. 为每组基因组构建 FM-index 索引。
2. 将切片序列比对到组内所有基因组索引。
3. 过滤比对所得锚点，完善比对结果。
4. 合并所有比对结果为 hal 文件。

#### 合并阶段

1. 按指导树完成全部比对，生成所有 hal 文件。
2. 按指导树顺序合并 hal 文件。

---

## 当前优先任务（🔴🟠🟢）

### 对接架构开发
- [x] 🔴在rare_aligner.h中新增MultipeRareAligner类，用于多基因组比对，PairRareAligner暂时放弃开发。
- [ ] 🔴在anchor.h中完成AnchorPtrListVec和MUMmer的delta格式的互相转换，完成后两个人可以各自独立开发，zpl使用MUMmer的结果进行all2all的开发，tqz完成整体项目架构的开发
### 索引和比对功能开发
- [ ] 🟠在rare_aligner.h中MultipeRareAligner中增加buildIndex功能，完成所有基因组fmidex索引的构建
- [ ] 🟠在rare_aligner.h中MultipeRareAligner中增加all-all比对函数，完成簇内高效并行all-all比对
- [ ] 🟢在data_process.h中对所有基因组完成[重复掩蔽](https://github.com/BioinformaticsToolsmith/Red),并确保被掩蔽的序列不在参考序列中参与比对
- [ ] 🟢在index.h中优化 MUM（最大唯一匹配）搜索性能

### 锚点过滤功能开发
- [ ] 🔴在anchor.h中参考 MUMmer4 mgaps.cc 实现锚点聚类这个类的实现，并完成锚点过滤的功能。
- [ ] 🔴支持 MUMmer delta 输出（结果转换为delta格式，完成两人的对接），支持 lastz lav/axt 输出
- [ ] 🟢 实现自定义比对格式，开发自定义格式与其他格式的转换工具

### 完善整体框架
- [ ] 🟢 实现进化树读取与解析，构建比对指导树
- [ ] 🟢 支持 group 模式（基于聚类的渐进式比对）和支持 all 模式（全 all-all 比对）
- [ ] 🟢 实现预处理模块子程序 RaMA-preprocess
- [ ] 🟢 实现比对模块 RaMA-G（子命令 build、align、map）
- [ ] 🟢 实现合并模块 RaMA-merge
- [ ] 🟠 完善单元测试模块

### all-all 合并阶段开发
- [ ] 🔴 完成rare_aligner.h中mergeAll2AllHSP，确认后续的开发任务

### 性能与结构优化

- [ ] 🟢优化 MUM 搜索相关算法
- [ ] 🟠优化 fasta\_manager 类，将清洗与索引功能解耦

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
