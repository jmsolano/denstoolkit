# DensToolKit

This is a first development distribution to the suite DensToolKit. For the time being, only private access will be granted.

## Git instructions:

Once you have recived an invitation you are able to see this page. Now you can install DensToolKit as follows. In your bash terminal type:

```
$cd /local/path/for/dtk
$git clone https://username@bitbucket.org/jmsolano/denstoolkitdevbb.git
```

Your bitbucket password will be requested to clone the repository. This is because this repository is not public.

After this, git will transfer the source files to ```/local/path/for/dtk/denstoolkit```

Optionally, you may want to first run
```
/local/path/for/dtk/src/checkdependencies
```
This will check that auxiliary programs (to plot the data computed by DensToolKit) are installed
in your system. You can skip this step if you are only interested in obtaining the numerical 
data. You may process the data later with your preferred plotting program.


For building the binaries, type:

```
$cd denstoolkit/src
$make
$sudo make install
```

This should compile and install the binaries into ```/usr/local/bin```

If you want to test the correct compilation of the program, type:

```
$make runtest
```
The above command will use the binaries present in /local/path/for/dtk/bin, therefore, the
tests can be made without running "sudo make install".

## Updating DensToolKit

If you only need to get the latest version of DensToolKi, and you are not using personal modifications to the source, in the local main directory (using the above example should be ```\local\path\for\dtk\denstoolkit```):

```
$git pull
$make distclean
$make
$sudo make install
```

This should update your binaries in the local installation directory.


#Compiling with OpenMP

DensToolKit contains paralellized implementation of almost all basic functions. For using the parallel version, you must edit the Makefile by changing the line

```
SETDTKNPROC=1
```

To use N processors, the above line should look like this:

```
SETDTKNPROC=N
```

In the current version, it is not advised to use more than 4 processors. There will be no speed improvement using more processors. Also, less than 3 makes no difference in the processing times.



