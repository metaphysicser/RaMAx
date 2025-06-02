#!/bin/bash
set -e
set -o pipefail

name=Mammals

package=package$name
seqs=sim$name.seqs.tar.gz
annots=sim$name.annots.tar.gz
ancestor=sim$name.ancestor.maf.gz
burnin=sim$name.burnin.maf.gz
noparalogy=sim$name.noparalogyMafs.tar.gz

# curl -O "https://cgl.gi.ucsc.edu/data/alignathon/{README.txt,data/sim$name.annots.tar.gz,data/sim$name.seqs.tar.gz,data/sim$name.ancestor.maf.gz,data/sim$name.burnin.maf.gz,data/sim$name.noparalogyMafs.tar.gz}"
echo "bddd7ab44c51b45f79380f190fd7dfa0  sim$name.annots.tar.gz
a554a2151b3bbe269c2dcf6e07030ab7  sim$name.seqs.tar.gz
4bab2832a972a26a9a43af150096295e  sim$name.ancestor.maf.gz
0a4c595644a806e7342ec3be62893f39  sim$name.burnin.maf.gz
bffc18321f937a1eac183c946049d190  sim$name.noparalogyMafs.tar.gz" > md5sum.txt
md5sum --check md5sum.txt
rm md5sum.txt
mkdir -p $package/annotations $package/predictions $package/sequences $package/truths $package/regions
mv README.txt $package/

mv $seqs $package/sequences
pushd $package/sequences
tar -xvzf $seqs
rm $seqs
popd

mv $annots $package/annotations
pushd $package/annotations
tar -xvzf $annots
rm $annots
popd

mv $ancestor $package/truths
mv $burnin $package/truths
mv $noparalogy $package/truths
pushd $package/truths
gunzip -c $ancestor > ancestor.maf
gunzip -c $burnin > burnin.maf
tar -xvzf $noparalogy
rm $noparalogy
touch version_3
popd
