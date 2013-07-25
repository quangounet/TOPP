// -*- coding: utf-8 -*-
// Copyright (C) 2013 Quang-Cuong Pham <cuong.pham@normalesup.org>
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

#include <iostream>
#include <vector>
#include <list>
#include <math.h>

typedef double dReal;

#define TINY 1e-10;

namespace TOPP {


enum IntegrationReturnType {
    RT_MAXSTEPS = 0,
    RT_END = 1,
    RT_BOTTOM = 2,
    RT_MVC = 3
};

enum SwitchPointType {
    SP_TANGENT = 0,
    SP_SINGULAR = 1,
    SP_DISCONTINUOUS = 2
};


class Trajectory {
public:
    size_t dimension;
    dReal duration;
    void Eval(dReal t, std::vector<dReal>& q);
    void Evald(dReal t, std::vector<dReal>& qd);
    void Evaldd(dReal t, std::vector<dReal>& qdd);
};


class Tunings {
public:
    dReal mvctimestep;
    dReal integrationtimestep;
};


class SwitchPoint {
    dReal s;
    int switchpointtype;
};


class Profile {
private:
    std::vector<dReal> svect, sdvect, sddvect;
    dReal integrationtimestep;
public:
    Profile(std::list<dReal>& slist, std::list<dReal>& sdlist, std::list<dReal>&  sddlist, dReal integrationtimestep);
    dReal duration;
    int nsteps;
    bool Evalall(dReal t, dReal& s, dReal& sd, dReal& sdd);
};


class Constraints {
public:
    Constraints(Trajectory trajectory, Tunings tunings);

    Trajectory trajectory;
    Tunings tunings;
    size_t dimension;

    virtual void SampleDynamics(){
        throw "Virtual method not implemented";
    }
    virtual void ComputeMVC(){
        throw "Virtual method not implemented";
    }
    virtual void ComputeSwitchPoints(){
        throw "Virtual method not implemented";
    }

    dReal SdLimitMVC(dReal s);  // upper limit on sd given by the MVC
    dReal SdLimitDirect(dReal s);  // upper limit on sd given by the direct velocity constraints
    std::pair<dReal,dReal> SddLimits(dReal s, dReal sd);  // lower and upper limits on sdd
    void GetSwitchPoints(std::list<SwitchPoint>& switchpointslist);  // get the list of switch points
    void HandleSwitchPoint(SwitchPoint sp);
};

}
