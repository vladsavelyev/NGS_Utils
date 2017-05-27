$PYTHON setup.py install --single-version-externally-managed --root=/
# `--single-version-externally-managed --root=/` is added to make setuptools avoid downloading dependencies. This resolves the conda-build
#       error "Setuptools downloading is disabled in conda build. Be sure to add all dependencies in the meta.yaml"
chmod -R o+r $PREFIX/lib/python*/site-packages/ngs_utils*
#$PYTHON -c "from ngs_utils import version ; print(version.__version__)" > __conda_version__.txt
