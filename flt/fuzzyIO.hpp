/*  Copyright (C) 2004-2015
    ANTONIO JAVIER BARRAGAN, antonio.barragan@diesia.uhu.es
    http://www.uhu.es/antonio.barragan

    Collaborators:
    JOSE MANUEL ANDUJAR, andujar@diesia.uhu.es

    DPTO. DE ING. ELECTRONICA, DE SISTEMAS INFORMATICOS Y AUTOMATICA
    ETSI, UNIVERSITY OF HUELVA (SPAIN)

    For more information, please contact with authors.

    This software is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This software is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef _FUZZYIO_HPP_
#define _FUZZYIO_HPP_

/**
* \file fuzzyIO.hpp
* \brief Defines Fuzzy Logic Tools streams operators.
*/

#include <string.h>
#include <iostream>
#include <flt/system.hpp>

namespace std{
	/// Ostream operator for TYPE_MF enumeration
	ostream &operator<<(ostream &F, FLT::TYPE_MF &T);

	/// Istream operator for TYPE_MF enumeration
	istream &operator>>(istream &F, FLT::TYPE_MF &T);

	/// Ostream operator for Membership class
	ostream &operator<<(ostream &F, FLT::Membership &P);

	/// Istream operator for Membership class
	istream &operator>>(istream &F, FLT::Membership &P);

	/// Ostream operator for Rule class
	ostream &operator<<(ostream &F, FLT::Rule &R);

	/// Istream operator for Rule class
	istream &operator>>(istream &F, FLT::Rule &R);

	/// Ostream operator for System class
	ostream &operator<<(ostream &F, FLT::System &S);

	/// Istream operator for System class
	istream &operator>>(istream &F, FLT::System &S);
} // std

#endif
