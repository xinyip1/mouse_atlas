Jupyter notebook files are binary and not very suited for git tracking (see this detailed [guide](https://nextjournal.com/schmudde/how-to-version-control-jupyter)). Here, we propose a method where you track .html and .py files as well as the .ipynb files, in order to see the difference in versions more easily. To save these files automatically, move the jupyter_notebook_config.py file into the .jupyter/ directory under your home directory:

```
> mv jupyter_notebook_config.py ~jarny/.jupyter/
```

After that, whenever you save your notebook, it will automatically create .html and .py files, which you can track with git.

Note that this method is creating a global config for all your notebooks - ie. this script will run for all your notebook saves regardless of which repo the notebook is in.