# Hybrid-Sim
Beautiful python code to simulate hybrid engine performance

Unfortunately a virtual environment doesn't work because you need to install Fortran compilers on your system to compile RocketCEA.

**You need to install the following dependencies:**
-numpy
-matplotlib
-sympy
-rocketcea
-this list is probably incomplete


**Setup Instructions:** 
Install MinGW-w64 GCC compiler: http://mingw-w64.org/doku.php/download/mingw-builds
Follow setup instructions here for Windows: https://rocketcea.readthedocs.io/en/latest/installgfortran.html
Next add MinGW to the system's PATH environment variable --> C:\MinGW\mingw64\lib OR C:\MinGW\mingw32\lib must be added to system's PATH variable

AFTER you've done this, you can try to install RocketCEA --> pip install rocketcea
If you get an error, screw around with MinGW, try installing a different one, make sure you have gfortran installed.

install all dependencies --> pip install xxxxx

If you are using an IDE like PyCharm, make sure your .gitignore file is configured to ignore any IDE specific directories.
