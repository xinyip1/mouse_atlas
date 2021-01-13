## Setting up git lfs

[git lfs](https://git-lfs.github.com/) is used to support large files - usually data files which are too large to track . Here we recommend that you let git lfs manage all files inside received/. To do this, you only need run a couple of git lfs commands, as shown in [here](https://git-lfs.github.com/).

First install git lfs
```
> git lfs install
```

Then select a folder (or file pattern) where you'll store the large files. If you're using this template repository as your starting point, remember to delete the existing history of this by removing .git directory first, otherwise the repo will have memory of its existing large files.
```
> git lfs track 'received/**'
```

> Note that it appears that git knows about the original files in received/ even after .git/ has been removed and a fresh repo has been initialised (Paul discovered this 23/11/2020). So remove all files in received first if you're working with your own input files.

To see the file which are tracked by git lfs, use

```
> git lfs ls-files
```