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
#include <ctime>

namespace TOPP {



/*********************************************************************/
/********************************** Vertex ***************************/
/*********************************************************************/

Vertex::Vertex(){
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
dReal Vertex::area(std::pair<dReal,dReal> v1, std::pair<dReal,dReal> v2){
    dReal out=sqrt(pow(v1.first*v2.second - v1.second*v2.first,2));
    return out;
}

bool Vertex::expand(soplex::SoPlex& lp, Vertex* vres){
    // Calculate the unit vector normal to the edge [v,vnext]
    std::pair<dReal,dReal> vortho(next->y-y,x-next->x);
    dReal vnorm = sqrt(pow(vortho.first,2)+pow(vortho.second,2));
    // if(vnorm<1e-5) {
    //     expanded = true;
    //     return false;
    // }
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

    //check if expansion adds area to the polygon
    dReal abs;
    abs = Vertex::area(std::pair<dReal,dReal> (xopt - x, yopt - y),std::pair<dReal,dReal> (x - next->x, y - next->y));
    if(abs < 0.01) { //TODO: shouldn't be this constant parametrizible somewhere?
        expanded = true;
        return false;
    }

    // Create new vertex
    vres->x=xopt;
    vres->y=yopt;
    vres->next = next;
    vres->expanded = false;
    next = vres;
    expanded = false;
    return true;
}

std::string Vertex::toString()
{
    std::ostringstream s("");
    s << std::setprecision (11);
    s << x << " " << y;
    return s.str();
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
        ivertex++;
    }
    return true;
}

void Polygon::iter_expand(soplex::SoPlex& lp){
    int maxiter = 20;
    int niter = 0;
    Vertex* v = (vertices.front());
    Vertex* vnew;
    while((!all_expanded()) && niter < maxiter) {
        if(!v->expanded) {
            vnew = new Vertex();
            bool res = v->expand(lp, vnew);
            if(res) {
                vertices.push_back(vnew);
                niter++;
            }
        }
        else{
            v = (v->next);
        }
    }
}

std::string Polygon::toString(){
    std::ostringstream s("");
    Vertex* v = vertices.front();
    int n = 0;
    while(n < int(vertices.size())) {
        s << v->toString() << "\n";
        v = v->next;
        n++;
    }
    return s.str();
}


/*********************************************************************/
/********************************* Main ** ***************************/
/*********************************************************************/
bool read(soplex::SoPlex& mysoplex, std::istringstream& q, std::istringstream& G, std::istringstream& A, std::istringstream& b, std::istringstream& h){
    int numcol=0;
    double a;
    soplex::DSVector dummycol(0);
    while(q >> a) {
//        std::cout << "Reading column" << numcol <<"\n";
//        std::cout << " Value" << a << "\n";
        numcol++;
        mysoplex.addColReal(soplex::LPCol(a,dummycol,soplex::infinity,-soplex::infinity));
    }
//   std::cout << "Together " << numcol << " variables";


    double tmp;

    while(h >> a) {
        soplex::DSVector row1(numcol);
        for (int i=0; i<numcol; i++) {
            if(G >> tmp) {

                row1.add(i,tmp);
            }
            else{
                std::cerr << "unable to read constraint";
                return false;
            }
        }
        //       std::cout << "Added G-type constraint \n";
        mysoplex.addRowReal(soplex::LPRow(-soplex::infinity,row1,a));
    }

    while(b >> a) {
        soplex::DSVector row1(numcol);
        for (int i=0; i<numcol; i++) {
            if(A >> tmp) {

                row1.add(i,tmp);
            }
            else{
                std::cerr << "unable to read constraint";
                return false;
            }
        }
//        std::cout << "Added A-type constraint \n";
        mysoplex.addRowReal(soplex::LPRow(a,row1,a));
    }
    return true;
}

bool OptimizeDirection(std::pair<dReal,dReal> v, soplex::SoPlex& mysoplex, std::pair<double,double>&z){
    int numcol=mysoplex.numColsReal();
    soplex::SPxSolver::Status stat;
    soplex::DVector prim(numcol);
    // Modify the q matrix
    mysoplex.changeObjReal(numcol-2, -v.first);
    mysoplex.changeObjReal(numcol-1, -v.second);
    // Solve the LP
    stat = mysoplex.solve();
    // Return results
    if( stat == soplex::SPxSolver::OPTIMAL )
    {
        mysoplex.getPrimalReal(prim);
        z = std::pair<dReal,dReal>(prim[numcol-2],prim[numcol-1]);
        return true;
    }
    else{
        std::cerr << "Unable to optimize for direction: ("<< v.first<<", "<<v.second<<")\n";
        return false;
    }
}




bool ComputePolygon(const std::string& lp_q, const std::string& lp_G, const std::string& lp_A, const std::string& lp_b, const std::string& lp_h, std::string& resstring){

    std::chrono::time_point<std::chrono::system_clock> t0, t1, t2;
    std::chrono::duration<double> chronoduration;

    // Load the LP matrices
    t0 = std::chrono::system_clock::now(); // Start chrono

    std::istringstream q(lp_q);
    std::istringstream G(lp_G);
    std::istringstream A(lp_A);
    std::istringstream b(lp_b);
    std::istringstream h(lp_h);
    std::vector<dReal> tmpvect;

    soplex::SoPlex lp;
    lp.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MINIMIZE);
    lp.setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_WARNING);
    read(lp, q, G, A, b, h);

    t1 = std::chrono::system_clock::now(); // Finished reading

    // Compute the first three directions
    bool res;
    std::pair<dReal,dReal> z1;
    std::pair<dReal,dReal> z2;
    std::pair<dReal,dReal> z3;
    res = OptimizeDirection(std::pair<dReal,dReal>(1,0), lp, z1);
    if(!res)
        return false;
    res = OptimizeDirection(std::pair<dReal,dReal>(cos(2*M_PI/3),sin(2*M_PI/3)), lp, z2);
    if(!res)
        return false;
    res = OptimizeDirection(std::pair<dReal,dReal>(cos(4*M_PI/3),cos(4*M_PI/3)), lp, z3);
    if(!res)
        return false;

    // Make the first three vertices
    Vertex v1(z1.first,z1.second);
    Vertex v2(z2.first,z2.second);
    Vertex v3(z3.first,z3.second);
    v1.next = &v2;
    v2.next = &v3;
    v3.next = &v1;

    // Populate the list of vertices
    std::list<Vertex*> l;
    l.resize(0);
    l.push_back(&v1);
    l.push_back(&v2);
    l.push_back(&v3);

    // Create the polygon and expand
    Polygon p(l);
    p.iter_expand(lp);
    resstring = p.toString();


    t2 = std::chrono::system_clock::now(); // Finish polygon computation

    // Display
    chronoduration = t1-t0;
    //std::cout << "Matrix reading: " << chronoduration.count() << std::endl;
    chronoduration = t2-t1;
    //std::cout << "Polygon computation: " << chronoduration.count() << std::endl;

    return true;
}



}
