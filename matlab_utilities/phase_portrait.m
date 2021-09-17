% This script prints Phase Portrait for close loop system in range X1,X2=[-10:10]

Plant = readfis('Plant.fis');
Controller = readfis('Controller.fis');

% Other option:
% Plant = txt2fis('Plant.txt');
% Controller = txt2fis('Controller.txt');

accuracy = 0.5;
tend = 0.6;
Limit = 10;
X1 = -Limit:accuracy:Limit;
X2 = X1;

figure, hold on
X0 = [0;0];

for j = 1:length(X2)
    if (j==1 || j==length(X2))
        index = 1:length(X1);
    else
        index = [1,length(X1)];
    end
    for i = index
        X0 = [X1(i),X2(j)];
        [t,dX] = ode45(@fuzeval,[0,tend],X0,[],Plant,Controller);
        plot(dX(:,1),dX(:,2))
    end
end

plot([-Limit-1,Limit+1],[0,0],'k-.',[0,0],[-Limit-1,Limit+1],'k-.')
axis([min(X1)-1 max(X1)+1 min(X2)-1 max(X2)+1]);
xlabel('X1'),ylabel('X2')
title('Phase Portrait')

% #############################################################################
%
%    Copyright (C) 2004-2011
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

