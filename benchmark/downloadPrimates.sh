#!/bin/bash
set -e
set -o pipefail

name=Primates

package=package$name
seqs=sim$name.seqs.tar.gz
annots=sim$name.annots.tar.gz
ancestor=sim$name.ancestor.maf.gz
burnin=sim$name.burnin.maf.gz
noparalogs=sim$name.noparalogyMafs.tar.gz

# curl -O "https://cgl.gi.ucsc.edu/data/alignathon/{README.txt,data/sim$name.annots.tar.gz,data/sim$name.seqs.tar.gz,data/sim$name.ancestor.maf.gz,data/sim$name.burnin.maf.gz,data/sim$name.noparalogyMafs.tar.gz}"
# echo "7d337b5e4f7c6eeb8eeeda95c2c21271  sim$name.annots.tar.gz
# d817e8739c10a0ddfcbe37200545b7f9  sim$name.seqs.tar.gz
# 1e2417d2ae8b4cf2743d5e740b7c5ed3  sim$name.ancestor.maf.gz
# 4fb72a9f14cf016c0d7b906d25e4731f  sim$name.burnin.maf.gz
# 3fc4fcb8fa64958f2a9d655b992387f7  simPrimates.noparalogyMafs.tar.gz" > md5sum.txt

# md5sum --check md5sum.txt
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
pushd $package/truths
gunzip -c $ancestor > ancestor.maf
gunzip -c $burnin > burnin.maf
touch version_3
popd

mv $noparalogs $package/truths
pushd $package/truths
tar -xvzf $noparalogs
rm $noparalogs
popd
