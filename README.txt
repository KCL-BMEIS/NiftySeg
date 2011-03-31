#####################
# NIFTY_SEG PACKAGE #
#####################

##############################################################################

--------------------------------
1 WHAT DOES THE PACKAGE CONTAIN?
--------------------------------
The code contains programs to perform EM based image segmentation, both on it's standard form [2] and as described in [1].

This code also contains a local copy of the nifty-reg 1.3 package available at http://sourceforge.net/projects/niftyreg/. Both projects will (hopefully) merge at some point in the future. 


##############################################################################

-----------------------
2 HOW TO BUILD THE CODE
-----------------------
The code can be easily build using cmake (http://www.cmake.org/). The latest 
version can be downloaded from http://www.cmake.org/cmake/resources/software.html
Assuming that the code source are in the source path folder, you will have 
to first create a new folder, i.e. build path (#1) and then to change 
directory to move into that folder (#2).
#1 >> mkdir build_path 
#2 >> cd build_path 

There you will need to call ccmake (#3a) in order to edit in the 
build options. If you don't want to specify options, we could just use cmake 
(#3b) and the default build values will be used.
#3a >> ccmake source_path
#3b >> cmake source_path

Once all the options are properly chosen, just press "c" to configure the Make- 
File and then the "G" key to generate it. In the prompt, you just have to 
make (#4) and then make install (#5).
#4 >> make 
#5 >> make install 

##############################################################################

---------
3 LICENSE
---------
Copyright (c) 2009, University College London, United-Kingdom
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

##############################################################################

---------
5 CONTACT
---------
For any comment, please, feel free to contact M. Jorge Cardoso (manuel.cardoso@ucl.ac.uk).

##############################################################################

------------
6 REFERENCES
------------
[1] M. Jorge Cardoso, Matthew J. Clarkson, Gerard R. Ridgway, Marc Modat, Nick C. Fox, Sebastien Ourselin, Alzheimer's Disease
Neuroimaging Initiative, LoAd: A locally adaptive cortical segmentation algorithm, Neuroimage
[2] K. Van Leemput, F. Maes, D. Vandermeulen, P. Suetens, Automated Model-Based Tissue Classification of MR Images of the Brain, 
IEEE Transactions on Medical Imaging, vol. 18, no. 10, pp. 897-908, October 1999


##############################################################################
##############################################################################