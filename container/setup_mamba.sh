# set up miniforge
apt-get --allow-releaseinfo-change update && apt-get install -y curl
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p "/opt/conda"
rm Miniforge3-$(uname)-$(uname -m).sh
. /opt/conda/etc/profile.d/conda.sh
. /opt/conda/etc/profile.d/mamba.sh
mamba env update --name sc-dandelion-container -f environment.yml
echo ". /opt/conda/etc/profile.d/conda.sh" | tee -a $SINGULARITY_ENVIRONMENT
echo ". /opt/conda/etc/profile.d/mamba.sh" | tee -a $SINGULARITY_ENVIRONMENT
echo "conda activate sc-dandelion-container" | tee -a $SINGULARITY_ENVIRONMENT
