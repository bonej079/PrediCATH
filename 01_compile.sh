#!/usr/bin/env bash

# Replace with local anaconda path

. /home/jbon4/apps/anaconda3/etc/profile.d/conda.sh
conda activate PrediCATH

for i in ./GoGraph/classes/*.py; do mv -v "$i"{,x} ; done
for i in ./utilities/*.py; do mv -v "$i"{,x} ; done

python3 setupClasses.py build_ext --inplace
./02_patchdist.sh classes

python3 setupUtilities.py build_ext --inplace
./02_patchdist.sh utilities

ln -s ./GoGraph/classes/go_20191007.obo
ln -s ./GoGraph/classes/gene_ontology_edit_20160601.obo
ln -s ./GoGraph/classes/icd.pkl

conda deactivate