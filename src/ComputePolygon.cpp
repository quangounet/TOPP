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


#include "ComputePolygon.h"


namespace TOPP {



/*********************************************************************/
/********************************** Vertex ***************************/
/*********************************************************************/

Vertex::Vertex(){
    std::cout << "Error : should not get here\n";
}

Vertex::Vertex(dReal x0, dReal y0){
    x = x0;
    y = y0;
    expanded = false;
    next = NULL;
}

dReal Vertex::length(){
    return sqrt(pow(next->x-x,2)+pow(next->y-y,2));
}

bool Vertex::expand(glp_prob& lp, Vertex& vres){
    // Calculate the unit vector normal to the edge [v,vnext]
    std::pair<dReal,dReal> vortho(next->y-y,x-next->x);
    dReal vnorm = sqrt(pow(vortho.first,2)+pow(vortho.second,2));
    vortho.first = vortho.first/vnorm;
    vortho.second = vortho.second/vnorm;
    // Optimize in the direction of vortho
    std::pair<dReal, dReal> z;
    bool res = OptimizeDirection(vortho,lp,z);
    if(!res) {
        expanded = true;
        return false;
    }
    dReal xopt = z.first;
    dReal yopt = z.second;
    if(abs) {
        expanded = true;
        return false;
    }
    vres = *(new Vertex(xopt,yopt));
    vres.next = next;
    next = &vres;
    expanded = false;
    return true;
}



/*********************************************************************/
/********************************* Polygon ***************************/
/*********************************************************************/

Polygon::Polygon(std::list<Vertex*>& vertices0){
    vertices = vertices0;
}

bool Polygon::all_expanded(){
    std::list<Vertex*>::iterator ivertex = vertices.begin();
    while(ivertex!=vertices.end()) {
        if(!(**ivertex).expanded) {
            return false;
        }
    }
    return true;
}

void Polygon::iter_expand(glp_prob& lp){
    int maxiter = 15;
    int niter = 0;
    Vertex v = *(vertices.front());
    Vertex vnew;
    while((!all_expanded()) && niter < maxiter) {
        if(!v.expanded) {
            bool res = v.expand(lp, vnew);
            if(res) {
                vertices.push_back(&vnew);
                niter++;
            }
        }
        else{
            v = *(v.next);
        }
    }
}


/*********************************************************************/
/********************************* Main ** ***************************/
/*********************************************************************/

bool OptimizeDirection(std::pair<dReal,dReal> v, glp_prob& lp, std::pair<dReal,dReal>& z){
    // Modify the q matrix

    // Solve the LP

    // Return results
    return true;
}




std::vector<std::pair<dReal,dReal> > ComputePolygon(const std::string& lpstring){
    std::string buff;
    std::istringstream iss(lpstring);
    std::vector<dReal> tmpvect;

    glp_prob *lp;
    lp = glp_create_prob();

    while(iss.good()) {
        getline(iss, buff, '\n');
        // VectorFromString(buff, tmpvect);
    }
    // Build the matrices



}



}
