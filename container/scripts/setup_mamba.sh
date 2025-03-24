# set up miniforge
apt-get --allow-releaseinfo-change update && apt-get install -y curl git language-pack-en gcc
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh -b -p "/opt/conda"
rm Miniforge3-$(uname)-$(uname -m).sh
. /opt/conda/etc/profile.d/conda.sh
. /opt/conda/etc/profile.d/mamba.sh
# fix pathing
echo ". /opt/conda/etc/profile.d/conda.sh" | tee -a $SINGULARITY_ENVIRONMENT
echo ". /opt/conda/etc/profile.d/mamba.sh" | tee -a $SINGULARITY_ENVIRONMENT
echo "conda activate sc-dandelion-container" | tee -a $SINGULARITY_ENVIRONMENT
echo "export GERMLINE=/share/database/germlines/" | tee -a $SINGULARITY_ENVIRONMENT
echo "export IGDATA=/share/database/igblast/" | tee -a $SINGULARITY_ENVIRONMENT
echo "export BLASTDB=/share/database/blast/" | tee -a $SINGULARITY_ENVIRONMENT
echo "LANG=en_US.UTF-8" | tee -a $SINGULARITY_ENVIRONMENT
chmod +x /share/dandelion_preprocess.py
chmod +x /share/changeo_clonotypes.py
# install dependencies
mamba create -n sc-dandelion-container
mamba env update --name sc-dandelion-container -f environment.yml
mamba activate sc-dandelion-container
# Get igblast
URL="https://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/"
latest_version=$(curl -s $URL | grep -oE 'igblast-[0-9]+\.[0-9]+\.[0-9]+' | head -1 | awk -F'-' '{print $2}')
curl -L -O "ftp://ftp.ncbi.nih.gov/blast/executables/igblast/release/LATEST/ncbi-igblast-$latest_version-x64-linux.tar.gz"
tar -xzvf ncbi-igblast-$latest_version-x64-linux.tar.gz
mv ncbi-igblast-$latest_version /share/ncbi-igblast-$latest_version
rm ncbi-igblast-$latest_version-x64-linux.tar.gz
echo "export PATH=/share/ncbi-igblast-$latest_version/bin:$PATH" >/etc/profile.d/igblast.sh
echo "export IGBLAST_VERSION=$latest_version" >>/etc/profile.d/igblast.sh
chmod +x /etc/profile.d/igblast.sh
