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

### 对接架构开发
- [x] 🔴 在rare_aligner.h中新增MultipeRareAligner类，用于多基因组比对，PairRareAligner暂时放弃开发。
### 索引和比对功能开发
- [ ] 🟠 在rare_aligner.h中MultipeRareAligner中增加group比对函数，完成簇内使用星比对完成比对
- [ ] 🔴 在data_process.h中对所有基因组完成[重复掩蔽](https://github.com/BioinformaticsToolsmith/Red)
- [ ] 🟠 在index.h中完成fmindex的优化，支持拼接的序列用\01作为分割符号
- [ ] 🟠 确保被掩蔽的序列不在参考序列中参与比对
- [ ] 🟠 在index.h中完成fastamanger的优化，只有部分序列的fasta快速索引到原有位置
### 锚点过滤功能开发
- [ ] 🔴 在anchor.h中参考 MUMmer4 mgaps.cc 实现锚点聚类这个类的实现，并完成锚点过滤的功能。
- [ ] 🟢 支持 MUMmer delta 输出（结果转换为delta格式，完成两人的对接），支持 lastz lav/axt 输出
- [ ] 🟢 实现自定义比对格式，开发自定义格式与其他格式的转换工具

### 完善整体框架
- [ ] 🟢 实现进化树读取与解析，构建比对指导树
- [ ] 🟢 支持 group 模式（基于聚类的渐进式比对）和支持 star 模式（所有序列完成星比对）
- [ ] 🟢 实现预处理模块子程序 RaMA-preprocess
- [ ] 🟢 实现比对模块 RaMA-G（预留子命令 build、align、map）
- [ ] 🟢 实现合并模块 RaMA-merge
- [ ] 🟢 实现pipeline模块RaMA-prepare，输出所有命令完成比对
- [ ] 🟠 完善单元测试模块

### 星比对合并阶段开发
- [ ] 🔴 完成rare_aligner.h中starAlignGroup，确认后续的开发任务
- [ ] 🔴 完成初步的星比对
- [ ] 🔴 完成星比对结果的合并
- [ ] 🟢 完成迭代补充星比对结果

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
