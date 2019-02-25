
histogram_mrc.py
===========
**histogram_mrc.py** is a graphical python program which displays the
histogram of voxel intensities contained in an MRC file.
It can be useful when deciding what thresholds to use
with in the "**filter_mrc**" and "**combine_mrc**" programs.

Typical usage:
```
   histogram_mrc  file.mrc
```
or
```
   histogram_mrc  -n 100  file.mrc  
```
## Arguments:

### -n  num_bins

The "**-n num_bins**" argument allows you to specify the number of
bins in the histogram.
(The default number of bins sometimes looks poor,
especially when the majority of the voxel intensities are narrowly distributed.)

### -rescale

The "**-rescale**" argument will shift and multiply all of the voxel intensities
so that the minimum and maximum voxel intensities are 0 and 1, respectively.

## Installation:


**histogram_mrc.py** requires python 2.7 or higher and depends
on the **matplotlib** and **mrcfile** python modules.

These modules can be installed using
```
pip install matplotlib
pip install mrcfile
```
To be safe, you may want to run
the following commands
beforehand:
```
virtualenv venv
source venv/bin/activate
```
*(If you do this, then you will have to run
"source venv/bin/activate" beforehand every time you
want to use this software. The virtualenv tool is
[explained in detail here](http://docs.python-guide.org/en/latest/dev/virtualenvs/)
)*

The *histogram_mrc.py* program is not currently installable using pip.
To use *histogram_mrc.py*, copy the "histogram_mrc.py"
file to somewhere in your path (such as /usr/local/bin/)
```
sudo cp histogram_mrc.py /usr/local/bin
```

Alternatively, you can edit your PATH variable manually to include
the subdirectory where the "histogram_mrc.py" script is located.
Suppose the directory with this README file is named ``visfd/doc''
and it is located in your home directory:

If you use the bash shell, typically you would edit your 
`~/.profile`, `~/.bash_profile` or `~/.bashrc` files 
to contain the following lines:

    export PATH="$PATH:$HOME/visfd/bin/histogram_mrc"

