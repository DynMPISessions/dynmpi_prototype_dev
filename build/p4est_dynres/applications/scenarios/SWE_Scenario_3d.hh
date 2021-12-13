/**
 * @file three-dimensional simulation scenario.
 *
 * Adapted from the SWE_Scenario.hh, see this file 
 * for original author information.
 *
 */

#ifndef __SWE_SCENARIO3D_H
#define __SWE_SCENARIO3D_H


/**
 * enum type: numbering of the boundary edges
 */
typedef enum BoundaryEdge {
   BND_LEFT, BND_RIGHT, BND_BOTTOM, BND_TOP, BND_FRONT, BND_BACK
} BoundaryEdge;

/**
 * SWE_Scenario defines an interface to initialise the unknowns of a 
 * shallow water simulation - i.e. to initialise water height, velocities,
 * and bathymatry according to certain scenarios.
 * SWE_Scenario can act as stand-alone scenario class, providing a very
 * basic scenario (all functions are constant); however, the idea is 
 * to provide derived classes that implement the SWE_Scenario interface
 * for more interesting scenarios.
 */
class SWE_Scenario_3d {

 public:

    virtual float getWaterHeight(float x, float y, float z) { return 10.0f; };
    virtual float getVeloc_u(float x, float y, float z) { return 0.0f; };
    virtual float getVeloc_v(float x, float y, float z) { return 0.0f; };
    virtual float getVeloc_w(float x, float y, float z) { return 0.0f; };
    virtual float getBathymetry(float x, float y, float z) { return 0.0f; };
    
    virtual float endSimulation() { return end_time; };
    virtual void setSimulationEnd(float end){end_time = end;};

    virtual float getBoundaryPos(BoundaryEdge edge) {
       if (edge==BND_LEFT || edge==BND_BOTTOM)
          return 0.0f;
       else
          return 1.0f; 
    };
    virtual ~SWE_Scenario_3d() {};

   private:

   float end_time = 0.1;

};


#endif
