%
% You must have installed the Template Numerical Toolkit (http://math.nist.gov/tnt)
% and GNU GSL for Windows (http://gnuwin32.sourceforge.net/packages/gsl.htm) on
% your build path, or you can use '-I' option
%
% The Mex environment of MATLAB must be configured with gcc.
% See GNUMex (http://gnumex.sourceforge.net)
%
% GSL for Windows must be installed. 'libgsl.a' and 'libgslcblas.a' are in
% the 'lib' directory from GSL for Windows (usually C:\Program Files\GnuWin32\lib).
% You must add the path for '.a' files where they are used, for example:
% 'C:\Program Files\GnuWin32\lib\libgsl.a'
%
clc
directory = mfilename('fullpath');
directory = strcat(directory(1:end-length(mfilename)), '../src');
directory = cd(directory);

% Compile the base code
mex membership.cpp -c
mex rule.cpp -c
mex system.cpp -c
mex utilities.cpp -c
mex derivatives.cpp -c
mex Kalman.cpp -c -DHAVE_INLINE

!move /Y *.obj ..\matlab_utilities

cd ..\matlab_utilities

mex aux_matlab.cpp -c

% Compile and link MEX files
mex activation.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex antec2mat.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex aproxjac.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj derivatives.obj
mex aproxlinear.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj derivatives.obj
mex conseq2mat.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex extractPoints.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex initializeController.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex fis2txt.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex fuz2mat.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex fuzcomb.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex fuzderconseq.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj derivatives.obj
mex fuzderaffine.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj derivatives.obj
mex fuzderantec.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj derivatives.obj
mex fuzderparam.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj derivatives.obj
mex fuzeval.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex fuzjac.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj derivatives.obj
mex Kalmanconseq.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj derivatives.obj Kalman.obj 'C:\Archivos de programa\GnuWin32\lib\libgsl.a' 'C:\Archivos de programa\GnuWin32\lib\libgslcblas.a'
mex Kalmanantec.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj derivatives.obj Kalman.obj 'C:\Archivos de programa\GnuWin32\lib\libgsl.a' 'C:\Archivos de programa\GnuWin32\lib\libgslcblas.a'
mex Kalmanfuz.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj derivatives.obj Kalman.obj 'C:\Archivos de programa\GnuWin32\lib\libgsl.a' 'C:\Archivos de programa\GnuWin32\lib\libgslcblas.a'
mex fuzlinear.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj derivatives.obj
mex fuzprint.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex mat2antec.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex mat2conseq.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex mat2fuz.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex fuzsubsystem.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj
mex txt2fis.cpp membership.obj rule.obj system.obj utilities.obj aux_matlab.obj

% Delete the object code
!del *.obj -Q
cd(directory);clear directory

% #############################################################################
%
%    Copyright (C) 2004-2015
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
