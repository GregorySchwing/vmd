# QwikFold VMD plugin

(Contains a non-cached version of alphafold and a slightly modified run_alphafold.py script. We'll change it ASAP!!!)

```bash
#!/bin/bash
# QwikFold install instructions for AlphaFold
#
# Diego E. B. Gomes  | dgomes@auburn.edu
# Rafael C. Bernardi | rcbernardi@auburn.edu
#

# Install alphafold 2 - Ubuntu 18.04.
# https://github.com/deepmind/alphafold

# Step 0 - Modify these paths according
export ALPHAFOLD_DATASETS='/data/alphafold_dbs/'

# Step 1 - create a conda environment
conda create -n af2 python=3.8 -y

# Step 2 - activate the environment
conda activate af2

# Step 3 - Clone Alphafold from GitHub into THIS folder.
git clone https://github.com/deepmind/alphafold.git

# Step 4 - Download chemical properties to the common folder
wget -q -P alphafold/common/ https://git.scicore.unibas.ch/schwede/openstructure/-/raw/7102c63615b64735c4941278d92b554ec94415f8/modules/mol/alg/src/stereo_chemical_props.txt

# Step 5 - Go to AlphaFold folder, and set path to current DIR
cd alphafold
alphafold_path=${PWD}

# Step 6 - Install required packages
conda install -y -c nvidia cudnn==8.0.4
conda install -y -c bioconda hmmer hhsuite==3.3.0 kalign2
conda install -y -c conda-forge openmm=7.5.1 pdbfixer pip

# Step 7 - Upgrade PIP 
pip3 install --upgrade pip

# Step 8 - Install additional AlphaFold requirements using PIP
pip3 install -r ./requirements.txt
pip3 install --upgrade "jax[cuda111]" -f https://storage.googleapis.com/jax-releases/jax_releases.html

# Step 9 - Apply patch to OpenMM
python_path=$(which python)
cd $(dirname $(dirname ${python_path}))/lib/python3.8/site-packages
patch -p0 < ${alphafold_path}/docker/openmm.patch

# Step 10 (very slow) - Download AlphaFold datasets (reduced datasets ~/400Gb)
if  [ ! -d ${ALPHAFOLD_DATASETS} ] ; then 
  mkdir -p ${ALPHAFOLD_DATASETS}
fi

# Step 11 - Download Alphafold reduced datasets ( ~400Gb)
cd ${alphafold_path}
bash scripts/download_all_data.sh ${ALPHAFOLD_DATASETS} reduced_dbs

# Step 12 - Download Alphafold complete datasets ( ~2.2Tb )
cd ${alphafold_path}
bash scripts/download_all_data.sh ${ALPHAFOLD_DATASETS}
```
