Installation instructions for VISFD
==========================

This file explains how to install the VISFD tools
and their dependencies.

## Outline

1) Install VISFD
2) Install the python dependencies
3) Install PoissonRecon/SSDRecon
4) Windows and MacOS specific instructions

## Prerequisites:

Before you begin, you will need:

- A computer with a terminal running BASH (not TCSH).
  This is included in linux.
  Instructions for [#Windows-instructions] and [#Apple-MacOS-instructions]
  users is included below.
- At least 16GB of RAM memory.  (Preferrably 32GB or more.)
- A text editor.  (Such as vi, emacs, atom, VisualStudio, Notepad++, ...)
- A recent version of python (version 3.0 or later).
- The "pip3" (or "pip") python installer tool (usually included with python)
- A modern C++ compiler (C++17 or later).  This software has been tested
  using GCC (v9.3.0) and CLANG (v10.0.0).  Earlier versions may not work.
- OpenMP
  In linux (and Windows WSL) OpenMP support is included with the GCC and
  CLANG compilers by default.  *(I don't know if MacOS supports OpenMP.
  See [below.](#Apple-MacOS-instructions))*
- The "make" tool.
  (Type "make" into the terminal.  If it doesn't complain "command not found"
  then it is installed.)
- [IMOD/3dmod](https://bio3d.colorado.edu/imod/)
- [Meshlab](https://www.meshlab.net)
- git (Usually "git" can be installed using your OS' package management
  installer tool, such as "apt", "yum", or "brew", etc...)


### Optional, but recommended:

- [ChimeraX](https://www.rbvi.ucsf.edu/chimerax/download.html)
- A computer with *at least* 4 CPU cores (8 threads).
- A *hard drive*.  (An HDD, not an SSD.  Frequent rewrites to an SSD
  may reduce its lifetime.  An external/USB hard drive should be fast enough.)

*(This file does not cover installation of these prerequisites.)*

When you can satisfy these requirements, you are ready to proceed.


## STEP 1: Install the VISFD software tools

```
git clone https://github.com/jewettaij/visfd ~/visfd
```

If this fails because "git" is not installed, then:
- Open this web site with a web browser:
  https://github.com/jewettaij/visfd
- Click on the word "Releases" on the right-hand side of the page.
  (It's not at the top.  It's about half-way down the page on the right.)
- Click on the "Tags" tab button on the upper-left hand side of the page.
- Click on the top most icons with the word "zip" underneath it.
- Save the zip file.
- Unpack the zip file and save it in your home directory with the name "visfd".


## Now compile the VISFD source code

```
cd ~/visfd
```

The next step depends on which compiler you are using.

#### If you are using the GCC compiler, enter:
```
gcc --version
```
This will print the version of your compiler to the terminal.
If the version is earlier than 9.0.0, then stop.
Either update your compiler, or install CLANG.
If your compiler is 9.0.0 or later, then enter this:
```
source setup_gcc.sh
```
...and skip th the next section.

#### If you are using the CLANG compiler, enter this instead:
```
source setup_clang.sh
```

### Now compile the VISFD software using:

```
make clean
make
```

Copy the VISFD executable programs to somewhere in your PATH.
In these instructions, I assume you have a ~/bin directory in your path.
If you do not have a writeable directory in your PATH,
(or if you have no idea what I'm talking about), then enter this beforehand:

```
mkdir ~/bin
echo "export PATH=\"$PATH:$HOME/bin\" >> ~/.bashrc
```

Now copy these files to that ~/bin directory:
```
cp -f bin/filter_mrc/filter_mrc ~/bin/
cp -f bin/sum_voxels/sum_voxels ~/bin/
cp -r bin/combine_mrc/combine_mrc ~/bin/
cp -r bin/voxelize_mesh/voxelize_mesh.py ~/bin/
```

Then log out and log in again for the changes to take effect.



## STEP 2: Install the python modules that VISFD needs:

#### Optional, but recommended:

First setup a virtual environment where you can run the VISFD tools.
```
python3 -m venv ~/bin/venv_visfd    #(use "python" if "python3" is not found)
deactivate # Don't worry if this step prints an error message.
source ~/bin/venv_visfd/bin/activate
```

Then run these commands:
*(Note: If it complains "pip3: command not found",
then try using "pip" instead.)*

```
pip3 install numpy       # (Use "pip" if "pip3" fails.)
pip3 install mrcfile
pip3 install pyvista
pip3 install matplotlib  # optional but recommended
pip3 install skimage     # optional but recommended
pip3 install scipy       # optional but recommended
```


## STEP 3:  Install "PoissonRecon" (and "SSDRecon")

```
git clone https://github.com/mkazhdan/PoissonRecon ~/PoissonRecon
cd ~/PoissonRecon
```

*(Note: The code in that repository changes frequently.  If at some point,
the SSDRecon and PoissonRecon programs are no longer able to adequately
create closed surfaces, then try downloading the old version of that code in
[my backup fork](https://github.com/jewettaij/PoissonRecon)
instead.)*



If you are using the CLANG compiler, edit the "Makefile" file
and uncomment the line containing this text:
```
   #COMPILER ?= clang
```
Then enter:
```
make
```

Note: Compiling this software requires at least 12GB of RAM to compile.
If you run out of ram during compilation (ie. if your computer freezes),
then visit [this page](https://github.com/mkazhdan/PoissonRecon/issues/152).

The newly compiled programs are in the Bin/ subdirectory.
Copy them to somewhere in your PATH (ie. to your ~/bin/ directory).

```
cd Bin/*/
cp -f PoissonRecon ~/bin/
cp -f SSDRecon ~/bin/
```

Then exit this terminal.
(We need to restart the terminal for the changes in your PATH to take effect.)
```
exit
```


Installation is complete.



### *If you created a virtual environment*

If you followed the suggestion above and created a virtual environment
*(using "python -m venv ~/bin/venv_visfd")*, then you must do this
before using the VISFD tools:  Start a new terminal and enter:
```
source ~/bin/venv_visfd/bin/activate
```
You must do this every time before you run the VISFD software.
(Otherwise, the "*voxelize_mesh.py*" program will not work.)
Alternatively, you can add that command to your ~/.bashrc file this way:
```
echo "source ~/bin/venv_visfd/bin/activate" >> ~/.bashrc
```




## Windows instructions

It is recommended that you install the BASH shell environment on your computer,
along with *make* and either *gcc* or *clang*.  Once you have done that,
you can follow the instructions above for linux users.
There are several ways to to create a BASH environment,
but perhaps the easiest method is to install
[Windows Subsystem for Linux (WSL2)](https://docs.microsoft.com/en-us/windows/wsl/install-win10),
***or***
[virtualbox](https://www.virtualbox.org)
(In the later case, you will also need to install a linux distribution,
preferably with a lightweight
desktop such as [xubuntu](https://xubuntu.org).)
Alternatively, you can try 
[Hyper-V](https://www.nakivo.com/blog/run-linux-hyper-v/)
or (if you have an older version of windows)
[CYGWIN](https://www.cygwin.com/).

WSL and virtualbox are virtual machines that allow you to run an
alternate operating system from within windows.
In this case that operating system is linux.  The BASH shell and the
compiler tools that you need can be easily installed from within in linux.
Both WSL and virtualbox also create an alternate filesystem inside windows
where the linux operating system is stored.  Software (like *visfd/filter_mrc*)
that you download and install there can access the files in that filesystem.
So you may need to copy your tomograms and otherf files to this fileystem
beforehand.  If you **are using WSL or WSL2**, then you should
[use caution when using windows programs to edit files stored there](https://devblogs.microsoft.com/commandline/do-not-change-linux-files-using-windows-apps-and-tools/).
This includes text editors, image editors, and tomography processing software
(such as IMOD).
One possible way to avoid problems is to try to restrict yourself to using
programs which you downloaded and installed directly from within the
WSL or WSL2 environment *(if possible)*.  In particular, a couple of the
features of the *filter_mrc* program require you to learn and use a
(unix-style) text editor.  (Word, Wordpad, and Notepad will not work.)
Again, it is a good idea to install and run such programs from within WSL,
not windows.



## Apple MacOS instructions

You must install OpenMP, and a compiler that supports it.
Otherwise, this software will run at an intolerably slow speed.
Last I checked OpenMP is disabled in MacOS by default.
You must figure out how to install it.
*(If it helps, there are some old instructions for doing that
[here](https://iscinumpy.gitlab.io/post/omp-on-high-sierra/).)*

I do not have access to a computer running MacOS.
If you figure out how to install OpenMP, please contact me
(or create a pull-request) so that we can update this section.
