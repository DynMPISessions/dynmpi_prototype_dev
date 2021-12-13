/**
 * @file SSWE_simple_scenarios_3d.hh
 *
 * This file contains the 3d versions of the 
 * scenarios in SWE_simple_scenarios.hh
 */

#ifndef __SWE_SIMPLE_SCENARIOS_3D_H
#define __SWE_SIMPLE_SCENARIOS_3D_H

#define PI 3.14159265

#include <cmath>
#include <cfloat>
#include <vector>
#include <array>
//#include <netcdf.h>
#include "SWE_Scenario_3d.hh"


using namespace std;


/**
 * a simple radial dambreak scenario
 */
class SWE_RadialDamBreakScenario : public SWE_Scenario_3d {

  public:

    float getWaterHeight(float x, float y, float z) {
       return (( sqrt( (x-500.f)*(x-500.f) + (y-500.f)*(y-500.f) + (z-500.f)*(z-500.f) ) < 100.f ) ? 15.f: 10.0f);
    };

	virtual float endSimulation() { return (float) 200; };

    /** Get the boundary positions
     *
     * @param i_edge which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(BoundaryEdge i_edge) {
       if (i_edge == BND_LEFT)
         return (float)0;
       else if (i_edge == BND_RIGHT)
         return (float)1000;
       else if (i_edge == BND_BOTTOM)
         return (float)0;
       else if (i_edge == BND_TOP)
         return (float)1000;
       else if (i_edge == BND_BACK)
           return (float)0;
       else
           return (float)1000;
    };
};

#endif
