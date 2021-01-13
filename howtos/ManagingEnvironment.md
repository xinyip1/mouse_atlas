# Managing the environment

An environment is a collection of programs. [This page](https://angus.readthedocs.io/en/2019/conda_tutorial.html) has a nice introductory description about environments and conda. Here, we list some useful tips that the team has gathered over time, some of which aren't that easy to find.

## Recording your environment with sinfo (IV, 2020-11-23)

Results can change when you change software versions. It's good to keep a record of what environment you generated your results in, so you can recreate them in the future. [sinfo](https://gitlab.com/joelostblom/sinfo) can help with this, by reporting a manifest of your current environment inside your notebook.

<details>
<summary> Example usage (click to expand) </summary>

```python
import math

import natsort
import numpy
import pandas
from sinfo import sinfo

sinfo()
```

```
-----
natsort     5.3.3
numpy       1.17.3
pandas      0.25.1
sinfo       0.3.0
-----
Python 3.7.3 | packaged by conda-forge | (default, Dec  6 2019, 08:54:18) [GCC 7.3.0]
Linux-5.4.2-arch1-1-x86_64-with-arch
4 logical CPU cores
-----
Session information updated at 2019-12-14 16:14
```
</details>


## How to invoke conda activate automatically on Windows powershell (JC, 2020-11-20)

In modern windows, the default shell is powershell. When you open the terminal from VSCode, it starts the windows powershell but you can't run conda activate from there. To enable this, you have to run the following command as Administrator. Right click the Start Button, choose Windows Powershell (Admin) and run this command:

> set-executionpolicy remotesigned

After this, every time you open powershell, conda will automatically be activated on base (see https://answers.microsoft.com/en-us/windows/forum/all/whats-wrong-with-my-windows-powershell/f05e72f2-a429-4ee0-81fb-910c8c8a1306?auth=1)
