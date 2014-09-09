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
dReal Vertex::area(std::pair<dReal,dReal> v1, std::pair<dReal,dReal> v2){
    dReal out=sqrt(pow(v1.first*v2.second - v1.second*v2.first,2));
    return out;
}

  bool Vertex::expand(soplex::SoPlex& lp, Vertex* vres, std::chrono::duration<double>& duration){
  //    std::cout<<"expanding vertex:\n"<<toString()<<"\n";
  //std::cout<<"in direction: ";
    // Calculate the unit vector normal to the edge [v,vnext]
    std::pair<dReal,dReal> vortho(next->y-y,x-next->x);
    //std::cout << "("<<vortho.first<<", "<<vortho.second<<") normalized to";
    dReal vnorm = sqrt(pow(vortho.first,2)+pow(vortho.second,2));
    if(vnorm<1.e-16){
      std::cerr << "size is zero\n";
      exit(0);
    }

    vortho.first = vortho.first/vnorm;
    vortho.second = vortho.second/vnorm;
    //std::cout << "("<<vortho.first<<", "<<vortho.second<<")\n";
    // Optimize in the direction of vortho
    std::pair<dReal, dReal> z;
    bool res = OptimizeDirection(vortho,lp,z,duration);
    if(!res) {
      //std::cout << "expanding finished \n";
        expanded = true;
	//	exit(0);
        return false;
    }
    dReal xopt = z.first;
    dReal yopt = z.second;
    //std::cout << "abs is " <<abs<<"\n";
    //check if expansion adds volume to the polygon
    dReal abs;
    abs = Vertex::area(std::pair<dReal,dReal> (xopt - x, yopt - y),std::pair<dReal,dReal> (x - next->x, y - next->y));
    //std::cout << "area is:\n"<<abs<<"\n";
    if(abs < 0.01) {//TODO: shouldn't be this constant parametrizible somewhere?
        expanded = true;
	//std::cout << "vertex already expanded to maximum\n";
        return false;
    }
    //std::cout << "creating new vertex with (" << xopt << ", " << yopt << "):\n";
    // Vertex *ver = new Vertex(xopt,yopt);
    vres->x=xopt;
    vres->y=yopt;
    //    std::cout << "ver addr " << vres <<"\n";
    //    std::cout << "vertex created " << vres<<"\n";
    //std::cout << "current next is "<< next <<"\n";
    vres->next = next;
    //std::cout << "created next is "<< vres->next <<"\n";
    vres->expanded=false;
    //std::cout << vres->toString()<<"\n";
    next = vres;
    expanded = false;
    //std::cout << "updating current vertex:\n" << toString()<<"\n";
    
    return true;
}

  std::string Vertex::toString()
{
  std::ostringstream o("");
  o<<"("<<x<<", "<<y<<") to ("<<next->x<<", "<<next->y<<") and expanded:"<<expanded<<" next "<<next;
  return o.str();
}




/*********************************************************************/
/********************************* Polygon ***************************/
/*********************************************************************/

Polygon::Polygon(std::list<Vertex*>& vertices0){
    std::cout << "initializing\n";
    vertices = vertices0;

    std::cout << "initializing ended \n";
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

  void Polygon::iter_expand(soplex::SoPlex& lp,std::chrono::duration<double>& duration){
    int maxiter = 15;
    int niter = 0;
    Vertex* v = (vertices.front());
    Vertex* vnew;
    std::chrono::duration<double> dur_temp;
    while((!all_expanded()) && niter < maxiter) {
      //std::cout<<"Starting iteration "<<niter<<" with vertex "<< v->toString()<<"\n";
        if(!v->expanded) {
	  vnew = new Vertex();
	  bool res = v->expand(lp, vnew,dur_temp);
	  duration+=dur_temp;
            if(res) {
                vertices.push_back(vnew);
                niter++;
            }
            else {
	      //std::cerr << "polygon expansion failture";
            }
        }
        else{
            v = (v->next);
        }
	//niter++;
	//xstd::cout << "at iteration: "<<niter<<"polygon looks like:\n"<<toString()<<"\n";
    }
}

  std::string Polygon::toString(){
    std::ostringstream o("");
    o << "Polygon contains vertices:\n";
    for(Vertex* v : vertices){
      o <<"  "<< v->toString() << "\n";
    }
    return o.str();
  }


/*********************************************************************/
/********************************* Main ** ***************************/
/*********************************************************************/
bool read(soplex::SoPlex& mysoplex, std::istringstream& q, std::istringstream& G, std::istringstream& A, std::istringstream& b, std::istringstream& h){
  int numcol=0;
  double a;
    soplex::DSVector dummycol(0);
  while(q >> a){
//        std::cout << "Reading column" << numcol <<"\n";
//        std::cout << " Value" << a << "\n";
        numcol++;
        mysoplex.addColReal(soplex::LPCol(a,dummycol,soplex::infinity,-soplex::infinity));
   }
//   std::cout << "Together " << numcol << " variables"; 

   
   double tmp;

   while(h >> a){
        soplex::DSVector row1(numcol);
        for (int i=0;i<numcol;i++){
            if(G >> tmp){
                
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
   
   while(b >> a){
        soplex::DSVector row1(numcol);
        for (int i=0;i<numcol;i++){
            if(A >> tmp){
                
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

  bool OptimizeDirection(std::pair<dReal,dReal> v, soplex::SoPlex& mysoplex, std::pair<double,double>& z, std::chrono::duration<double>& duration){
    std::chrono::time_point<std::chrono::system_clock> t0,t1;
    int numcol=mysoplex.numColsReal();
    soplex::SPxSolver::Status stat;
    soplex::DVector prim(numcol);
    // Modify the q matrix
    mysoplex.changeObjReal(numcol-2, -v.first);
    mysoplex.changeObjReal(numcol-1, -v.second);
    // Solve the LP
    t0=std::chrono::system_clock::now();
    stat = mysoplex.solve();
    t1=std::chrono::system_clock::now();
    duration=t1-t0;
    // Return results
    if( stat == soplex::SPxSolver::OPTIMAL )
    {
        mysoplex.getPrimalReal(prim);
        /*std::cout << "LP solved to optimality.\n";
        std::cout << "Objective value is " << mysoplex.objValueReal() << ".\n";
        std::cout << "Primal solution is [" << prim[numcol-2]<<", "<< prim[numcol-1]<<"].\n";*/
        z = std::pair<dReal,dReal>(prim[numcol-2],prim[numcol-1]);
    }
    else{
        std::cerr << "unable to optimize for direction: ("<< v.first<<", "<<v.second<<")\n";
        return false;
    }
    return true;
}




bool  ComputePolygon(const std::string& lp_q,const std::string& lp_G,const std::string& lp_A,const std::string& lp_b,const std::string& lp_h, std::vector<std::pair<dReal,dReal> >& resvect){
    std::chrono::duration<double> duration, dur_temp;
    std::chrono::time_point<std::chrono::system_clock> t0,t1, t2;
    
    double tic_0 = clock();
    std::istringstream q(lp_q);
    std::istringstream G(lp_G);
    std::istringstream A(lp_A);
    std::istringstream b(lp_b);
    std::istringstream h(lp_h);
    std::vector<dReal> tmpvect;
    //read(lp, q,G,A,b,h);
    //itic_t = clock();nitialize the lp problem
 t0=std::chrono::system_clock::now();
    soplex::SoPlex lp;
    lp.setIntParam(soplex::SoPlex::OBJSENSE, soplex::SoPlex::OBJSENSE_MINIMIZE);
    lp.setIntParam(soplex::SoPlex::VERBOSITY, soplex::SoPlex::VERBOSITY_WARNING);
    read(lp, q,G,A,b,h);
    t2=std::chrono::system_clock::now();
    double tic_t = clock();
    //read(lp, q,G,A,b,h);

    

    //TODO: read the matrixes
    
    //compute the first three directions
    bool res;
    std::pair<dReal,dReal> z1;
    std::pair<dReal,dReal> z2;
    std::pair<dReal,dReal> z3;
    res = OptimizeDirection(std::pair<dReal,dReal>(1,0), lp, z1,dur_temp);
    duration=dur_temp;
    if(!res)
        return false;
    res = OptimizeDirection(std::pair<dReal,dReal>(cos(2*M_PI/3),sin(2*M_PI/3)), lp, z2,dur_temp);
    duration+=dur_temp;
    if(!res)
        return false;
    res = OptimizeDirection(std::pair<dReal,dReal>(cos(4*M_PI/3),cos(4*M_PI/3)), lp, z3,dur_temp);
    duration+=dur_temp;
    if(!res)
        return false;
    
//    std::cout << "pair1: (" << z1.first<<", "<< z1.second <<"), pair2: (" << z2.first<<", "<< z2.second <<"), pair3: (" << z3.first<<", "<< z3.second <<")\n";    
//    return true;

    //TODO: create initial vertexes
    Vertex v1(z1.first,z1.second);
    Vertex v2(z2.first,z2.second);
    Vertex v3(z3.first,z3.second);
    v1.next=&v2;
    // std::cout<<"first vector followed by" << v1.next << "\n";
    v2.next=&v3;
    v3.next=&v1;
    //    std::cout<<"created vertexes\n";
    //TODO: create polygon
    std::list<Vertex*> l;
    l.resize(0);
    //std::cout<<"list initialized\n";
    l.push_back(&v1);
    l.push_back(&v2);
    l.push_back(&v3);
    //std::cout<<"list filled\n";
    Polygon p(l);
    //std::cout<<"created polygon\n";
    
    //TODO: expand polygon
    p.iter_expand(lp,dur_temp);
    duration+=dur_temp;
 t1=std::chrono::system_clock::now();
    std::cout << "chrono lp time " << duration.count() << std::endl;
    std::cout << (clock() - tic_t)/CLOCKS_PER_SEC << " seconds" << std::endl;
    std::cout << (clock() - tic_0)/CLOCKS_PER_SEC << " seconds" << std::endl;
    duration=t1-t0;
    std::cout << "whole time: " << duration.count() << std::endl;
    duration=t2-t0;
    std::cout << "reading matrixes took "<<duration.count()<<std::endl;
    std::cout<<"polygon is"<<p.toString();
    std::cout<<"polygon expanded\n";
    return true;

    //TODO: return... confilcting types for Polygon and compute_polygon?
}



}
