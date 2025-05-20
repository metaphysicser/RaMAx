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

## 当前优先任务

### 功能开发

1. **完善比对过程（src/rare\_aligner.cpp, alignQueryFile）**

   * 增加反向链比对，正反链结果分开处理
   * 掩蔽序列不参与查询（后续可切换策略，建议多版本实现）
   * 优化 MUM（最大唯一匹配）搜索性能

2. **完善锚点过滤（src/rare\_aligner.cpp, filterAnchors）**

   * 参考 [MUMmer4 mgaps.cc](https://github.com/mummer4/mummer) 实现锚点聚类
   * 参考 [lastz](https://lastz.github.io/lastz/#stage_chaining) 或 MUMmer4，实现比对结果筛选
   * 比对结果支持输出 MUMmer delta、lastz lav/axt 及自定义格式
   * 实现自定义格式与其他格式的相互转换（开发独立可执行程序）

3. **完善整体框架**

   * 进化树读取与解析，构建渐进式比对指导树
   * 支持两种模式：基于聚类的渐进式比对（group 子命令）和全 all-all 比对（all 子命令）
   * 每个阶段开发对应子程序，形成标准 pipeline

     * 预处理阶段：RaMA-preprocess
     * 双基因组比对：RaMA-G（子命令：build、align、map，预留扩展）
     * 合并阶段：RaMA-merge
   * 补充完整的单元测试模块

4. **完成 all-all 合并阶段开发**

### 性能与结构优化

* 优化 MUM 搜索相关算法
* 优化 fasta\_manager 类，将清洗与索引功能解耦

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
