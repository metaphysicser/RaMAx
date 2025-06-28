#!/bin/bash

set -e  # 出错即退出
set -o pipefail

echo "[1/6] 编译 sonLib"
cd sonLib
make -j
cd ..

echo "[2/6] 编译 mafTools"
cd mafTools
make -j
cd ..

echo "[3/6] 编译 mwgAlignAnalysis"
cd mwgAlignAnalysis
make -j
cd ..

echo "[4/6] 编译 其他格式转换maf的软件"
g++ -std=c++11 -O3 -o delta2maf delta2maf.cpp
g++ -std=c++11 -O3 -o sam2maf sam2maf.cpp
g++ -std=c++11 -O3 -o paf2maf paf2maf.cpp

echo "[5/6] 编译 filter_maf_by_species"
g++ -std=c++11 -O3 -fopenmp -o filter_maf_by_species filter_maf_by_species.cpp

echo "[6/6] 收集所有可执行文件到 bin/"
mkdir -p bin
cp delta2maf sam2maf paf2maf filter_maf_by_species bin/
cp mafTools/bin/* bin/
cp ./mwgAlignAnalysis/evaluations/src/comparatorWrapper/comparatorSummarizer.py bin/


echo "✅ 所有步骤完成，可执行文件已放入 ./bin"
