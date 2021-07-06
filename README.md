# Poly6
----LICENCE---
This software is distributed under the MIT license

----USAGE----
In shell/cmd navigate to the location where both 'aPoly6.py' and 'mat000File.txt' are stored.
Execute:

python  aPoly.py

Notes:
--tested with Python 3.7.9 (should work with any version 3.*);
--requires the Python packages 'numpy' and 'matplotlib'.

----INPUT---
The script 'aPoly6.py' requires the input file 'mat000File.txt'. 
Th latter contains the material data. 
The provided sample contains instructions and can be used as template.
The rest of the provided 'mat*.txt' files, e.g., 'matDP980.txt', 'matAA2090T3.txt', etc, 
are containers for the materials tested in the article supported by this code 
('A parameter identification scheme for the orthotropic Poly6 yield function satisfying the convexity condition', 
preprint at http://dx.doi.org/10.13140/RG.2.2.18552.78082). 
Simply copy the content of any of these files and replace with it (paste over) the content of 'mat000File.txt' 
to obtain the corresponding Poly6-model.


----OUTPUT---
If a solution is found, the script will generate (and save on disk, at the same location):
1) a data file containing diagnosis and Poly6-coefficients (as in Tables 3 and 2 of the article);
2) plot of directional properties;
3) plot of several sections through the yield surface.     
