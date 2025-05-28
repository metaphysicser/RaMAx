
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
# 进入 mafTools 目录并编译
cd mafTools
make -j

# 返回上一级目录
cd ..

# 编译 mwgAlignAnalysis
cd mwgAlignAnalysis
make -j
```


