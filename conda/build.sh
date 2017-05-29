# Commands to build and upload the packages for osx and linux, py2 and py3

conda build ipython-cluster-helper --python 2.7
conda build ipython-cluster-helper --python 3.6
cd /Users/vlad/miniconda3/conda-bld
conda convert -p linux-32 osx-64/ipython-cluster-helper-0.5.4-py36_0.tar.bz2
conda convert -p linux-64 osx-64/ipython-cluster-helper-0.5.4-py36_0.tar.bz2
conda convert -p linux-32 osx-64/ipython-cluster-helper-0.5.4-py27_0.tar.bz2
conda convert -p linux-64 osx-64/ipython-cluster-helper-0.5.4-py27_0.tar.bz2
anaconda upload osx-64/ipython-cluster-helper-0.5.4-py36_0.tar.bz2
anaconda upload osx-64/ipython-cluster-helper-0.5.4-py27_0.tar.bz2
anaconda upload linux-32/ipython-cluster-helper-0.5.4-py36_0.tar.bz2
anaconda upload linux-32/ipython-cluster-helper-0.5.4-py27_0.tar.bz2
anaconda upload linux-64/ipython-cluster-helper-0.5.4-py36_0.tar.bz2
anaconda upload linux-64/ipython-cluster-helper-0.5.4-py27_0.tar.bz2

conda build ngs_utils -c bcbio -c bioconda -c conda-forge --python 2.7
conda build ngs_utils -c bcbio -c bioconda -c conda-forge --python 3.6
cd /Users/vlad/miniconda3/conda-bld
conda convert -p linux-32 osx-64/ngs_utils-1.1.3-py36_0.tar.bz2
conda convert -p linux-64 osx-64/ngs_utils-1.1.3-py36_0.tar.bz2
conda convert -p linux-32 osx-64/ngs_utils-1.1.3-py27_0.tar.bz2
conda convert -p linux-64 osx-64/ngs_utils-1.1.3-py27_0.tar.bz2
anaconda upload osx-64/ngs_utils-1.1.3-py36_0.tar.bz2
anaconda upload osx-64/ngs_utils-1.1.3-py27_0.tar.bz2
anaconda upload linux-32/ngs_utils-1.1.3-py36_0.tar.bz2
anaconda upload linux-32/ngs_utils-1.1.3-py27_0.tar.bz2
anaconda upload linux-64/ngs_utils-1.1.3-py36_0.tar.bz2
anaconda upload linux-64/ngs_utils-1.1.3-py27_0.tar.bz2
