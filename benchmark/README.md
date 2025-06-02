
### ✅ 环境配置说明（Python 2.7 + numpy + scipy）

#### 📦 方式一：使用 Conda（推荐）

```bash
# 创建 Python 2.7 环境并安装依赖
conda create -n benchmark-env python=2.7 numpy scipy -y

# 激活环境
conda activate benchmark-env
```

---

#### 🐍 方式二：使用 pip + virtualenv（适用于没有 conda 的情况）

```bash
# 安装 virtualenv（如尚未安装）
pip install virtualenv

# 创建并激活 Python 2.7 虚拟环境
virtualenv -p python2.7 benchmark-env
source benchmark-env/bin/activate

# 安装依赖
pip install numpy scipy
```

### ⚙️ 编译项目

确保你已经激活了上面创建的虚拟环境（`benchmark-env`），然后执行以下步骤来编译依赖工具和分析主程序：

```bash
bash setup_benchmark.sh
```
运行后，相关可执行文件和脚本将被编译并放置在 `bin/` 目录下。

### 📦 获取评估数据集

项目提供了用于下载基准测试数据的脚本。请运行以下命令：

```bash
# 下载哺乳动物测试数据
bash ./downloadMammals.sh

# 下载灵长类测试数据
bash ./downloadPrimates.sh
```

这些命令将在当前目录下分别生成两个数据文件夹：

* `packageMammals/`
* `packagePrimates/`

### 📁 数据集结构说明

每个数据包文件夹遵循如下结构，用于 Alignathon 分析工具评估比对质量：

```
package/
├── README.txt
├── annotations/     # 内部使用，不要修改
├── predictions/     
├── sequences/       # 测试集使用的原始参考序列
├── truths/          # 已知的真实对齐（用于评估）
```


### 🔍 获得双基因组比对结果文件并进行评估

本节以 `simHuman` 和 `simChimp` 为例，说明如何获取真实比对结果与预测比对结果，并使用评估工具计算比对准确性。

---

#### ✅ 1. 提取真实比对结果（Ground Truth）

下载的数据集提供的 `ancestor.maf` 文件是多个物种的祖先对齐结果。你可以使用 `filter_maf_by_species` 工具提取目标物种之间的真实对齐关系：

```bash
./bin/filter_maf_by_species \
    --ref simHuman \
    --query simChimp \
    -i ./packagePrimates/truth/ancestor.maf \
    -o ./simHuman-simChimp-truth.maf
```

该命令将从 `ancestor.maf` 中提取 `simHuman` 和 `simChimp` 的比对区域，生成一个仅包含它们的 MAF 文件，作为后续评估的“参考答案”。

---

#### 🔁 2. 获取预测比对结果

如果使用 **MUMmer** 进行比对，其输出为 `.delta` 文件。需要通过 `delta2maf` 工具转换为标准 MAF 格式，才能参与评估：

```bash
./bin/delta2maf \
    --ref ./packagePrimates/sequences/simHuman.fa \
    --query ./packagePrimates/sequences/simChimp.fa \
    --input ./simHuman-simChimp.delta \
    --output ./simHuman-simChimp-mummer.maf
```

对于其他比对工具（如 LASTZ），如果输出已是 MAF 格式，可直接使用。

---

#### 📊 3. 评估预测比对质量

使用 `mafComparator` 工具对比预测比对与真实比对，生成 XML 格式的比对评估文件：

```bash
./bin/mafComparator \
    --maf1 ./simHuman-simChimp-truth.maf \
    --maf2 ./simHuman-simChimp-mummer.maf \
    --out mummer.xml
```

接着，使用 Python 2.7 环境下的 `comparatorSummarizer.py` 脚本解析 XML 输出，提取关键性能指标（如准确率、召回率、F1 值）：

```bash
# 请确保已激活 Python 2.7 虚拟环境 benchmark-env
python ./bin/comparatorSummarizer.py --xml mummer.xml
```

输出将显示准确率，召回率和 F1 值等指标，这些指标可以作为评估不同比对工具性能的依据。



