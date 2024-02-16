### Set up conda environment with alphafold/colabfold and necessary dependencies.
conda create -y -n alphafold python=3.9 pip
conda activate alphafold
condapip=${CONDA_PREFIX}/bin/pip
$condapip install --no-warn-conflicts 'colabfold[alphafold-minus-jax] @ git+https://github.com/sokrypton/ColabFold'
sed -i 's/weights = jax.nn.softmax(logits)/logits=jnp.clip(logits,-1e8,1e8);weights=jax.nn.softmax(logits)/g' ${CONDA_PREFIX}/lib/python3.9/site-packages/alphafold/model/modules.py
conda install -y -c conda-forge -c bioconda kalign2=2.04 hhsuite=3.3.0 openmm=7.7.0 python=3.9 pdbfixer
