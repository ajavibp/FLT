%
% You must have installed the Template Numerical Toolkit (http://math.nist.gov/tnt),
% the GNU GSL (http://www.gnu.org/software/gsl) and Fuzzy Logic Tools on your build
% path, or you can use '-I' option
%
% GSL can be installed with 'libgsl0-dev' and 'libgsl0ldbl' packages in debian based distributions.
%
% The Mex environment of MATLAB must be configured with gcc (mex -setup).
%

% Install 'libgsl-dev'

% Add FLT and TNT headers path
setenv('CPATH',[':/home/javi/Datos/Investigacion/Código/include:/home/javi/Datos/Investigacion/Código/fuzzylogictools/flt/']);

directory = mfilename('fullpath');
directory = strcat(directory(1:end-length(mfilename)), '../src');
directory = cd(directory);

% Compile the base code

mex -largeArrayDims membership.cpp -c
mex -largeArrayDims rule.cpp -c
mex -largeArrayDims system.cpp -c
mex -largeArrayDims fuzzyIO.cpp -c
mex -largeArrayDims utilities.cpp -c
mex -largeArrayDims derivatives.cpp -c
mex -largeArrayDims Kalman.cpp -c -DHAVE_INLINE

!mv -f *.o ../matlab_utilities

cd ../matlab_utilities

mex -largeArrayDims aux_matlab.cpp -c

% Compile and link mex -largeArrayDims files

mex -largeArrayDims activation.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims antec2mat.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims aproxjac.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o derivatives.o
mex -largeArrayDims aproxlinear.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o derivatives.o
mex -largeArrayDims conseq2mat.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims extractPoints.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims initializeController.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims fis2txt.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims fuz2mat.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims fuzcomb.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims fuzderconseq.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o derivatives.o
mex -largeArrayDims fuzderaffine.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o derivatives.o
mex -largeArrayDims fuzderantec.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o derivatives.o
mex -largeArrayDims fuzderparam.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o derivatives.o
mex -largeArrayDims fuzeval.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims fuzjac.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o derivatives.o
mex -largeArrayDims fuzlinear.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o derivatives.o
mex -largeArrayDims fuzprint.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims Kalmanantec.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o derivatives.o Kalman.o -lgsl -lm -lgslcblas % Other options: -latlas -lcblas lgslcblas
mex -largeArrayDims Kalmanconseq.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o derivatives.o Kalman.o -lgsl -lm -lgslcblas % Other options: -latlas -lcblas lgslcblas
mex -largeArrayDims Kalmanfuz.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o derivatives.o Kalman.o -lgsl -lm -lgslcblas % Other options: -latlas -lcblas lgslcblas
mex -largeArrayDims mat2antec.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims mat2conseq.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims mat2fuz.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims fuzsubsystem.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o
mex -largeArrayDims txt2fis.cpp membership.o rule.o system.o utilities.o fuzzyIO.o aux_matlab.o

% Delete the object code
!rm *.o
cd(directory);clear directory

% #############################################################################
%
%    Copyright (C) 2004-2022
%    ANTONIO JAVIER BARRAGAN, antonio.barragan@diesia.uhu.es
%    http://uhu.es/antonio.barragan
%
%    Collaborators:
%    JOSE MANUEL ANDUJAR, andujar@diesia.uhu.es
%
%    DPTO. DE ING. ELECTRONICA, DE SISTEMAS INFORMATICOS Y AUTOMATICA
%    ETSI, UNIVERSITY OF HUELVA (SPAIN)
%
%    For more information, please contact with authors.
%
%    This software is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This software is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
% #############################################################################
