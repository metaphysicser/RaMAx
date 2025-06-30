# simHuman-simChimp结果

基于 `xml` 文件，通过 `comparatorSummarizer.py` 脚本生成，展示了以 human 为参考基因组，chimp 为查询基因组时的比对评估结果。

## 运行命令

```bash
python ./mwgAlignAnalysis/evaluations/src/comparatorWrapper/comparatorSummarizer.py --xml result.xml
````

## 指标说明

* **Precision**（精确率）：在所有被识别为匹配的区域中，真正正确的比例。
* **Recall**（召回率）：在所有真实匹配区域中，被成功识别的比例。
* **F-score**：精确率与召回率的调和平均数，是一个综合评价指标。
* **TP (A)**：参考基因组（human）中正确匹配的区域数。
* **TP (B)**：查询基因组（chimp）中正确匹配的区域数。
* **FP (B)**：查询基因组中错误地匹配为正例的区域数。
* **FN (A)**：参考基因组中未能识别出的真实正例区域数。

下面补充了用来生成上述对齐结果的具体命令示例，你可以根据实际文件名和参数做相应调整：

---

## 对齐工具运行命令

### 1. MUMmer (`nucmer` + `delta2maf`)

```bash
# 1) 用 nucmer 比对 human（参考） vs chimp（查询）
nucmer \
  -t 32 \                             # 使用 16 个线程
  -p simHuman-simChimp \              # 输出前缀
  simHuman.fa \                       # 参考基因组 FASTA
  simChimp.fa                         # 查询基因组 FASTA

# 2) 过滤 Delta（取最佳链）
delta-filter -1 simHuman-simChimp.delta > simHuman-simChimp.filter.delta

# 3) 转换为 MAF
./delta2maf \
  --ref simHuman.fa \
  --query simChimp.fa \
  --input simHuman-simChimp.filter.delta \
  --output simHuman-simChimp.nucmer.maf
```

### 2. wfmash (`wfmash` + `paf2maf`)

```bash
# 1) 用 wfmash 直接输出 PAF
wfmash \
  -t 32 \                             # 使用 32 个线程
  simHuman.fa \                       # 参考基因组 FASTA
  simChimp.fa \                       # 查询基因组 FASTA
  > simHuman-simChimp.wfmash.paf

# 2) 转换为 MAF
./paf2maf \
  --ref simHuman.fa \
  --query simChimp.fa \
  --input simHuman-simChimp.wfmash.paf \
  --output simHuman-simChimp.wfmash.maf
```
### 3. lastz
```bash
lastz  "/simHuman.fa[multiple]" 
	"simChimp.fa[multiple]" 
	--format=maf --notransition 
	--step=20   --chain --gapped 
	--progress=1 > lastz.maf
```

### 4. minimap2 (`minimap2` + `sam2maf`)
```bash
minimap2 -t 20 -ax asm5 --eqx 
simHuman.fa simChimp.fa > simHuman-simChimp.minimap2.sam

# 转换为 MAF
./sam2maf \
  --ref simHuman.fa \
  --query simChimp.fa \
  --input simHuman-simChimp.minimap2.sam \
  --output simHuman-simChimp.minimap2.maf
```


## 评估指标对比

### MUMmer 结果

| 比较方式               | Precision | Recall  | F-score | TP (A) | TP (B) | FP (B) | FN (A) |
| ------------------ | --------- | ------- | ------- | ------ | ------ | ------ | ------ |
| Overall (w/o self) | 0.59915   | 0.99941 | 0.74916 | 999347 | 598999 | 400757 | 594    |
| Overall (w/ self)  | 0.59915   | 0.99941 | 0.74916 | 999347 | 598999 | 400757 | 594    |
| simChimp-simHuman  | 0.59915   | 0.99941 | 0.74916 | 999347 | 598999 | 400757 | 594    |

### wfmash 结果
| 比较方式               | Precision | Recall  | F-score | TP (A) | TP (B) | FP (B) | FN (A) |
| ------------------ | --------- | ------- | ------- | ------ | ------ | ------ | ------ |
| Overall (w/o self) | 0.59775   | 0.99844 | 0.74780 | 998156 | 598270 | 402608 | 1559   |
| Overall (w/ self)  | 0.59775   | 0.99844 | 0.74780 | 998156 | 598270 | 402608 | 1562   |
| simChimp-simHuman  | 0.59775   | 0.99844 | 0.74780 | 998156 | 598270 | 402608 | 1559   |

### minimap2 结果

| 比较方式               | Precision | Recall  | F-score | TP (A) | TP (B) | FP (B) | FN (A) |
| ------------------ | --------- | ------- | ------- | ------ | ------ | ------ | ------ |
| Overall (w/o self) | 0.01431   | 0.02371 | 0.01785 | 23692  | 14311  | 986005 | 975421 |
| Overall (w/ self)  | 0.01431   | 0.02371 | 0.01785 | 23692  | 14311  | 986005 | 975425 |
| simChimp-simHuman  | 0.01431   | 0.02371 | 0.01785 | 23692  | 14311  | 986005 | 975421 |

### RaMA-G 结果（当前版本）

| 比较方式               | Precision | Recall  | F-score | TP (A) | TP (B) | FP (B) | FN (A) |
| ------------------ | --------- | ------- | ------- | ------ | ------ | ------ | ------ |
| Overall (w/o self) | 0.60367   | 0.97430 | 0.74546 | 604552 | 972441 | 396916 | 25649  |
| Overall (w/ self)  | 0.60367   | 0.97430 | 0.74546 | 604552 | 972441 | 396916 | 25653  |
| simChimp-simHuman  | 0.60367   | 0.97430 | 0.74546 | 604552 | 972441 | 396916 | 25649  |

