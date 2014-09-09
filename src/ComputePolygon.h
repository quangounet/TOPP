// -*- coding: utf-8 -*-
// Copyright (C) 2014 Quang-Cuong Pham <cuong.pham@normalesup.org>
//
// This file is part of the Time-Optimal Path Parameterization (TOPP) library.
// TOPP is free software: you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef ComputePolygon_H
#define ComputePolygon_H


#include "TOPP.h"
#include "soplex.h"

namespace TOPP {

class Vertex {
public:
    Vertex();
    Vertex(dReal x0, dReal y0);
    dReal x,y;
    bool expanded;
    Vertex* next;
    dReal length();
    bool expand(soplex::SoPlex& lp, Vertex* vresg, std::chrono::duration<double>& duration);
    dReal area(std::pair<dReal,dReal> v1, std::pair<dReal,dReal> v2);
    std::string toString();

};

class Polygon {
public:
    Polygon(std::list<Vertex*>& vertices0);
    std::list<Vertex*> vertices;
    Vertex* v1, v2, v3;
    bool all_expanded();
    void iter_expand(soplex::SoPlex& lp, std::chrono::duration<double>& duration);
    std::string toString();
};


bool OptimizeDirection(std::pair<dReal,dReal> v, soplex::SoPlex& lp, std::pair<dReal,dReal>& z, std::chrono::duration<double>& duration);

 bool ComputePolygon(const std::string& lp_q,const std::string& lp_G,const std::string& lp_A,const std::string& lp_b,const std::string& lp_h, std::vector<std::pair<dReal,dReal> >& resvect);

}


#endif
