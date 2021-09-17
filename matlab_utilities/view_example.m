clc, clear all, close all

% This script simulates a closed-loop system (Plant + Controller) usign fuzzy models
% Plant = readfis('Plant.fis');
% Controller = readfis('Controller.fis');

Plant = txt2fis('Plant.txt'); % The previous models in TXT format
Controller = txt2fis('Controller.txt');

tf = 3.5;
X0 = [10;10];

[t,X] = ode45(@fuzeval,[0,tf],X0,[],Plant,Controller); % Closed loop

subplot(221)
plot(t,X(:,1)),title('x_1(t)'),box off
xlabel('Time (s)'),ylabel('x_1(t)')
subplot(223),title('x_2(t)')
plot(t,X(:,2)),box off
xlabel('Time (s)'),ylabel('x_2(t)')

% Control Signals
U = zeros(2,size(X,1));
for i=1:size(X,1)
	U(:,i) = fuzeval(X(i,:)',Controller);
end
subplot(122)
plot(t,U(1,:),t,U(2,:)),box off
title('Control Signals')
legend('u_1(t)','u_2(t)')
xlabel('Time (s)'),ylabel('u(t)')
set(gcf, 'Position', get(0,'Screensize'));
% #############################################################################
%
%	Copyright (C) 2004-2011
%	ANTONIO JAVIER BARRAGAN, antonio.barragan@diesia.uhu.es
%	http://uhu.es/antonio.barragan
%
%	Collaborators:
%	JOSE MANUEL ANDUJAR, andujar@diesia.uhu.es
%
%	DPTO. DE ING. ELECTRONICA, DE SISTEMAS INFORMATICOS Y AUTOMATICA
%	ETSI, UNIVERSITY OF HUELVA (SPAIN)
%
%	For more information, please contact with authors.
%
%	This software is free software: you can redistribute it and/or modify
%	it under the terms of the GNU General Public License as published by
%	the Free Software Foundation, either version 3 of the License, or
%	(at your option) any later version.
%
%	This software is distributed in the hope that it will be useful,
%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%	GNU General Public License for more details.
%
%	You should have received a copy of the GNU General Public License
%	along with this program.  If not, see <http://www.gnu.org/licenses/>.
% #############################################################################
