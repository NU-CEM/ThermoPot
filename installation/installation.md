ThermoPot can be installed using the Python package manader `pip`:

`pip install thermopot`

If you use conda/anaconda, the safest thing to do is to create a new environment and then install ThermoPot:

```
conda create -n thermopot python
conda activate thermopot
pip install thermopot
```

If you wish, you can install the latest developer version of ThermoPot from Github. For this you will may want to specify a particular branch. Note that although the latest Github versions may include more features, it may not be stable. To install the latest, possibly unstable version:

```
git clone https://github.com/NU-CEM/ThermoPot.git
cd ThermoPot
git checkout <branch>
pip install .
```
