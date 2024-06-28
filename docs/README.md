Required software installed as follows:

```shell
$ mamba create -n docs python=3.12
$ mamba activate docs
$ pip install sphinx==7.3.7 sphinx_rtd_theme==2.0.0
```

To test document build:

```shell
$ sphinx-build -M html docs/ .test-docs/
```