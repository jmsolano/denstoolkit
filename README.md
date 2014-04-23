# Welcome

This is a first development distribution to the suite DensToolKit. For the time being, only private access will be granted.

## Git instructions:

After your forking this distribution, in your bash terminal type (note that jmsolano should be replaced by your own account login):

```
$cd /local/path/for/dtk
$ git clone https://jmsolano@bitbucket.org/jmsolano/denstoolkitdevbb.git/wiki denstoolkit
```

This will transfer the source files to ```/local/path/for/dtk/denstoolkit``` Now type:

```
$cd denstoolkit/src
$make
$make install
```

This should compile and install the binaries into ```/usr/local/bin```

If you want to test the correct compilation of the program, type:

```
$make runtest
```

## Updating DensToolKit



