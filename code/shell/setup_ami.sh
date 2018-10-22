# This script is adapted from that used by the BEBi103 Course taught by Justin Bois at Caltech
# http://www.github.com/justinbois/bebi103_ami_setup.git
#
sudo yum -y groupinstall "Development Tools";
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh;
bash ./miniconda.sh -b -p $HOME/miniconda;
export PATH="$HOME/miniconda/bin:$PATH";
echo 'export PATH="$HOME/miniconda/bin:$PATH"' >> ~/.bashrc;
source ~/.bashrc;
conda update -q -y conda;
conda install -y numpy cython pystan;
conda config --add channels conda-forge;
conda config --add channels defaults;
conda install -y scipy pandas numba scikit-image scikit-learn statsmodels bokeh altair holoviews watermark tqdm matplotlib seaborn ipython jupyter jupyterlab nodejs xarray netcdf4;
jupyter labextension install --no-build jupyterlab_bokeh;
jupyter labextension install --no-build @ijmbarr/jupyterlab_spellchecker;
jupyter labextension install --no-build @jupyter-widgets/jupyterlab-manager;
jupyter labextension install --no-build @pyviz/jupyterlab_pyviz;
jupyter labextension install --no-build  jupyterlab_vim;
jupyter lab build;
echo 'export LSCOLORS="gxfxcxdxCxegedabagacad"' >> ~/.bashrc;