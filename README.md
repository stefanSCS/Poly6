#  Poly6

This code generates the polynomial Poly6-model of an orthotropic sheet metal based uniaxial and biaxial mechanical test data. Theory and algorithm are presented in:

- Preprint (first draft): *A parameter identification scheme for the orthotropic Poly6 yield function satisfying the convexity condition* (posted at [ResearchGate](http://dx.doi.org/10.13140/RG.2.2.18552.78082))

-     

There are two versions (with details described next):

- Version one (**aPoly6**) uses pure Python (only Python needs be installed).
- Version two (**aPoly6_v2**) uses external C functions for the computationally intensive parts (requires Python and a C-compiler; Linux comes packed with both; Windows10 users can also use the Linux-shell or the provided binary).



## Version specific notes

- Both versions require the external Python packages [numpy](https://numpy.org/install/) and [matplotlib](https://matplotlib.org/stable/users/installing.html).

- Both versions use the file 'mat000File.txt' (discussed below) as input for material data; it must be located in the same directory as 'aPoly.py' or 'aPoly6_v2.py' (or else, edit the latter files appropriately).  

### Version: **aPoly6_v2** (updates **aPoly6**)

- Focus on speed: the computationally intensive parts (convexity check and sections plots) are now implemented in C (file 'cFuncs.c');

- The functions 'cvxCheck' and 'plotPoly6' were reimplemented to call via the 'ctypes' module the corresponding C-functions.

- USAGE: The C-functions are implemented in 'cFuncs.c'. In order for 'aPoly6_v2.py' to be able to call them, 'cFuncs.c' must be compiled as a shared library. In shell (command prompt), navigate to your working directory (the location of 'aPoly6_v2.py') and compile  'cFuncs.c' as a shared object via:
```c
gcc -fPIC -shared -o cFuncs.so cFuncs.c
```
Then simply execute:
```python
python aPoly6_v2.py
```

### Version: **aPoly6** (initial)

- USAGE: In shell/cmd navigate to the location where both 'aPoly6.py' and 'mat000File.txt' are stored. Execute:
```python
python aPoly.py
```


## Input

Both scripts, 'aPoly6.py' and 'aPoly6_v2.py' require the input file 'mat000File.txt'.
Here the material data is recorded. The provided sample contains instructions and can be used as template. The rest of the provided 'mat*.txt' files, e.g., 'matDP980.txt', 'matAA2090T3.txt', etc, are containers for the materials used as tests (in the cited article). Simply copy the content of any of these files and replace with it (paste over) the content of 'mat000File.txt' to obtain the corresponding Poly6-model.

A special note regarding the shape parameter **bxShape** in 'mat000File.txt' is in order. This parameter controls the convexity of the biaxial (zero shear) plane stress section. Its theoretical values are within the interval [0,1] but, in general, the null value is not compatible with the overall convexity of the yield surface. If 'mat000File.txt' does not define **bxShape** (its value is '*') one of the following default values will be used:

- AL: 0.01
- Fe: 0.35

Aluminum alloys are on the boundary of the Poly6 convexity domain (on the 'tight' side) and for them the smallest admissible value for **bxShape** is 0.01. Should the solution found (for 0.01) be unacceptable, this value may be increased incrementally to 0.014, 0.018, etc, until a satisfactory solution is found. (This can be automated also, if desired, by enclosing the current algorithm into a loop on **bxShape**). Locations where convexity is not achieved are reported in the console (via spherical coordinates angles) and this info can be used to decide whether to modify any default values.   

## Output

If a solution is found (the provided data may not always be compatible with convexity), the script will generate (and save on disk, at the same location):

- a data file containing diagnosis and Poly6-coefficients (as in Tables 3 and 2 of the article);
- plot of directional properties;
- plot of several sections through the yield surface.


## License

This software is distributed under the MIT license.
