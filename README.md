# NIFTY_SEG PACKAGE #

--------------------------------
1 WHAT DOES THE PACKAGE CONTAIN?
--------------------------------
This project, lead by Jorge Cardoso at King's College London, contains programs to perform EM based segmentation of images in nifti or analyse format. NiftySeg is an open-source toolkit licensed under the BSD license. It also contains a package of label fusion algorithms (MV, STAPLE, SBA) with different types of ranking strategies. Features: LoAd: Locally Adaptive Brain Segmentation General Purpose EM segmentation Single and Multi-label Fusion package To download the latest version, please check out the code by copying the following line to the termini;:

git clone git@github.com:KCL-BMEIS/NiftySeg.git niftyseg

A packaged stable release is also available in the files menu above. This release in only updated once in a while, thus, it does not have the latest developments. See below for installation instructions. 

-----------------------
2 HOW TO BUILD THE CODE
-----------------------

**Download**
The code can be easily build using cmake (http://www.cmake.org/). The latest 
version can be downloaded from http://www.cmake.org/cmake/resources/software.html

To download the latest version, please check out the code by copying the following line to the terminal:

~~~
git clone git@github.com:KCL-BMEIS/NiftySeg.git niftyseg
~~~

A packaged stable release is also available at this website. 

**Linux & OSX**

* **Build**

Assuming that the code source are in the source path folder, you will have to ﬁrst create a new folder, i.e. build path (step 1) and then to change directory to move into that folder (step 2).

~~~
mkdir build_path

cd build_path 
~~~

There you will need to call ccmake (step 3a) in order to ﬁll in the 
build options. If you don’t want to specify options, we could just use cmake 
(step 3b) and the default build values will be used.

~~~
ccmake source_path
cmake source_path
~~~

The main option in the ccmake gui are deﬁned bellow:
~~~
> CMAKE BUILD INSTALL options are Release, RelWithDebInfo or Debug 
> INSTALL_PRIORS Will install the population atlas for the segmentation pipeline
> INSTALL_PRIORS_DIRECTORY Directory where the population atlas is going to be installed
> INSTALL_NIFTYREG Will fetch and automatically configure and install the niftireg package. 
~~~

Once all the ﬂags are properly ﬁlled in, just press the ”c” to conﬁgure the Make- 
ﬁle and then the ”g” key to generate them.

* **Install**

In the prompt, you just have to make (step 4) ﬁrst and then make install (step 5).

~~~
make 
make install
~~~

**Windows**

* **Build**

The building process is the following:
1. Get the source
2. Create a new directory for the build: "niftyseg-build"
3. Launch CMake-Gui, set the source path to "niftyseg" and the build path to "niftyseg-build" then hit configure
4. Cmake will prompt you to select the generator, which means you'll need to select the Visual Studio version you have installed earlier
5. Make sure the Use OpenMP option is enabled.
6. Set the CMAKE_INSTALL_PREFIX to the folder where you want to install NiftySeg.
7. Note, that if you want to install NiftySeg under Program Files, you'll need to create the folder yourself and explicitly apply full write permissions. 
8. Once the flags are set, hit configure and generate. This will generate the Visual Studio project files.

* **Install**

1. Go to "niftyseg-build", and launch NiftySeg.sln. This will start Visual Studio.
2. In Visual Studio select build type, for generic use select Release and build the project (hit F7). 
3. Once the build finished Select and run the Install task (Right Click on Install > Project Only > Build only Install). This will install NiftySeg to the folder you selected earlier.
4. Probably you'll want to add the install folder to your system path.

---------
3 LICENSE
---------
Copyright (c) 2018, NiftySeg Development Team, United-Kingdom
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.
Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

Neither the name of the University College London nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
THE POSSIBILITY OF SUCH DAMAGE.

---------
5 CONTACT
---------
For any comments, please, feel free to contact M. Jorge Cardoso (manuel.cardoso@kcl.ac.uk).

