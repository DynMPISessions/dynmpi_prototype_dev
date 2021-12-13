/**
 * @file
 * This file is part of SWE.
 *
 * @author Michael Bader, Kaveh Rahnema, Tobias Schnabel
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 * @author Praktikum Tsunami-Simulation WS 20 Gruppe 4
 *
 * @section LICENSE
 *
 * SWE is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SWE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SWE.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 *
 * This file includes various scenarios for simulating.
 */

#ifndef __SWE_SIMPLE_SCENARIOS_H
#define __SWE_SIMPLE_SCENARIOS_H

#define PI 3.14159265

#include <cmath>
#include <cfloat>
#include <vector>
#include <array>
#ifdef NETCDF
#include <netcdf.h>
#endif
#include "SWE_Scenario.hh"


using namespace std;

#ifdef NETCDF
/**
 * @brief Read bathymetry and displacement data from two netCDF files
 *
 * This scenario reads two netCDF files containing batymetry data and the
 * initial displacement caused by an earthquake. Additionally, the boundary
 * conditions and simulation time can be set by the user.
 *
 * @author Maximilian Streubel, Atamert Rahma, Daniel Strauss
 *
 * @see SWE_Scenario, swe_dimensionalsplitting.hh
 */
class SWE_TsunamiScenario : public SWE_Scenario {

    public:

        /**
         * create a tsunami scenario
         *
         * @param bathymetry path to a file containing bathymetry data
         * @param displacement path to a file containing displacement data,
         *        may have smaller dimensions than the bathymetry file
         * @param i_simulation_time time to simulate
         * @param b_top top boundary condition
         * @param b_right right boundary condition
         * @param b_bot bottom boundary condition
         * @param b_left left boundary condition
         */
        SWE_TsunamiScenario(const char* bathymetry, const char* displacement, float i_simulation_time,
                            BoundaryType b_top, BoundaryType b_right, BoundaryType b_bot, BoundaryType b_left) {
            read_netcdf(bathymetry, displacement);
            simulation_time = i_simulation_time;
            bnd_top = b_top;
            bnd_right = b_right;
            bnd_bottom = b_bot;
            bnd_left = b_left;
        }

        ~SWE_TsunamiScenario(){
            delete[] x;
            delete[] y;
            for(unsigned int i = 0; i < n_y; i++){
                delete[] h[i];
            }
            delete[] buf;
            delete[] h;
        }

        float getBathymetry(float i_x, float i_y) {
            int index_x, index_y;
            getIndices(i_x, i_y, index_x, index_y);
            if(index_x > -1 && index_y > -1){
                return b[index_y][index_x];
            }
            return 0;
        }

        float getWaterHeight(float i_x, float i_y) {
            int index_x, index_y;
            getIndices(i_x, i_y, index_x, index_y);
            if(index_x > -1 && index_y > -1){
                return h[index_y][index_x];
            }
            return 0;
        }

    virtual float endSimulation() { return simulation_time; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) {
        switch(edge){
        case BND_TOP:
            return bnd_top;
        case BND_RIGHT:
            return bnd_right;
        case BND_BOTTOM:
            return bnd_bottom;
        default:
            return bnd_left;
        }
    }

    /** Get the boundary positions
     *
     * @param i_edge which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(BoundaryEdge i_edge) {
       if ( i_edge == BND_LEFT ){
         return x[0];
        }else if ( i_edge == BND_RIGHT){
         return x[n_x-1];
       }else if ( i_edge == BND_BOTTOM ){
         return y[0];
       }else{
         return y[n_y-1];
       }
    };

    private:

        size_t n_x;
        size_t n_y;
        size_t nd_x;
        size_t nd_y;

        float *x;
        float *y;
        float **b;
        float **h;
        float *buf;

        float simulation_time;

        BoundaryType bnd_top;
        BoundaryType bnd_right;
        BoundaryType bnd_bottom;
        BoundaryType bnd_left;

        /**
         * read the input files and create scenario with the bathymetry and
         * initial displacement read, calculate initial water height.
         *
         * @param bathymetry the name of the file containing the bathymetry data
         * @param displacement the name of the file containing the initial displacement data
         */
        void read_netcdf(const char* bathymetry, const char* displacement) {

            //read bathymetry data

            //open nc file
            int ncid_b;
            if(nc_open(bathymetry, NC_NOWRITE, &ncid_b))
                cout << "\t\tERROR opening " << bathymetry << "\n";

            //read dimensions
            int dimid;

            //read x dimension
            nc_inq_dimid(ncid_b, "x", &dimid);
            nc_inq_dimlen(ncid_b, dimid, &n_x);

            //read y dimension
            nc_inq_dimid(ncid_b, "y", &dimid);
            nc_inq_dimlen(ncid_b, dimid, &n_y);

            x = new float[n_x];
            y = new float[n_y];
            h = new float*[n_y];
            b = new float*[n_y];
            for(unsigned int i = 0; i < n_y; i++){
                h[i] = new float[n_x];
                b[i] = new float[n_x];
            }
            //create matrix of correct size.
            buf = new float[n_x*n_y];  // a contiguous block of memory

            //get varid of x
            int varid;
            if(nc_inq_varid(ncid_b, "x", &varid))
                cout << "\t\tERROR get varid x\n";

            //read x
            if(nc_get_var_float(ncid_b, varid, &x[0]))
                cout << "\t\tERROR reading x\n";

            //get varid of y
            if(nc_inq_varid(ncid_b, "y", &varid))
                cout << "\t\tERROR get varid y\n";

            //read y
            if(nc_get_var_float(ncid_b, varid, &y[0]))
                cout << "\t\tERROR reading y\n";

            //get varid of b
            if(nc_inq_varid(ncid_b, "z", &varid))
                cout << "\t\tERROR get varid b\n";

            //read b
            if(nc_get_var_float(ncid_b, varid, &buf[0]))
                cout << "\t\tERROR reading b\n";

            //close input file
            nc_close(ncid_b);

            //copy bathymetry and calc height
            for(unsigned int i = 0; i < n_y; i++){
                for(unsigned int j = 0; j < n_x; j++){
                    b[i][j] = buf[i*n_x + j];
                    if(-20 < b[i][j] && b[i][j] < 0){
                        b[i][j] = -20;
                    } else if(0 <= b[i][j] && b[i][j] < 20){
                        b[i][j] = 20;
                    }
                }
            }


            //calc height
            for(unsigned int i = 0; i < n_y; i++){
                for(unsigned int j = 0; j < n_x; j++) {
                    if(b[i][j] > 0){
                        h[i][j] = 0;
                    } else {
                        h[i][j] = -b[i][j];
                    }
                 }
            }


            //initial displacement


            //open nc file
            int ncid_d;
            if(nc_open(displacement, NC_NOWRITE, &ncid_d))
                cout << "\t\tERROR opening" << displacement << "\n";


            //read dimensions

            //read x dimension
            nc_inq_dimid(ncid_d, "x", &dimid);
            nc_inq_dimlen(ncid_d, dimid, &nd_x);

            //read y dimension
            nc_inq_dimid(ncid_d, "y", &dimid);
            nc_inq_dimlen(ncid_d, dimid, &nd_y);

            float x_d[nd_x];
            float y_d[nd_y];
            float d[nd_y][nd_x];

            //get varid of x
            if(nc_inq_varid(ncid_d, "x", &varid))
                cout << "\t\tERROR get varid x\n";

            //read x
            if(nc_get_var_float(ncid_d, varid, &x_d[0]))
                cout << "\t\tERROR reading x\n";

            //get varid of y
            if(nc_inq_varid(ncid_d, "y", &varid))
                cout << "\t\tERROR get varid y\n";

            //read y
            if(nc_get_var_float(ncid_d, varid, &y_d[0]))
                cout << "\t\tERROR reading y\n";

            //get varid of d
            if(nc_inq_varid(ncid_d, "z", &varid))
                cout << "\t\tERROR get varid d\n";

            //read d
            if(nc_get_var_float(ncid_d, varid, &d[0][0]))
                cout << "\t\tERROR reading d\n";

            //close input file
            nc_close(ncid_d);

            //map displacement to bathymetry and update b
#pragma omp parallel for schedule(static)
            for(unsigned int i = 0; i < nd_y; i++) {
                for(unsigned int j = 0; j < nd_x; j++) {
                    int index_x, index_y;
                    getIndices(x_d[j], y_d[i], index_x, index_y);
                    if(index_x > -1 && index_y > -1){
                        b[index_y][index_x] += d[i][j];
                    }
                }
            }
        }

        /**
         * map x and y coordinates to the corresponding cell indices.
         * If no suitable cell was found, the indices will be set to -1.
         *
         * @param i_x the x coordinate
         * @param i_y the y coordinate
         * @param index_x will be set to the cells x index
         * @param index_y will be set to the cells y index
         */
        void getIndices(float i_x, float i_y, int &index_x, int &index_y) {
            index_x = -1;
            //find neighbor indices
            for(unsigned int i = 0; i < n_x - 1; i++) {
                if(i_x >= x[i] && i_x <= x[i+1]) {
                    //find closest neighbor
                    if(i_x -x[i] < x[i+1] - i_x) {
                        index_x = i;
                    } else {
                        index_x = i+1;
                    }
                }
            }

            index_y = -1;
            for(unsigned int i = 0; i < n_y - 1; i++) {
                if(i_y >= y[i] && i_y <= y[i+1]) {
                    if(i_y -y[i] < y[i+1] - i_y) {
                        index_y = i;
                    } else {
                        index_y = i+1;
                    }
                }
            }
        }
};

#endif


/**
 * @brief Artificial Tsunami
 *
 * This scenario constructs an artificial tsunami consisting of
 * constant bathymetry with a small displacement in the center of the domain
 * and constant water height.
 *
 * @author Maximilian Streubel, Atamert Rahma, Daniel Strauss
 *
 * @see SWE_Scenario
 */
class SWE_ArtificialTsunamiScenario : public SWE_Scenario {

  public:

    float getBathymetry(float x, float y) {
       x -= 5000;
       y -= 5000;
       if((x < -500 || x > 500) || (y < -500 || y > 500))
           return -100.0;
       return -100.0 + (5* (sin(((x/500.0) + 1) * 3.14159265)) * (-(y/500.0)*(y/500.0) + 1));
    };

    float getWaterHeight(float x, float y) {
       return 100.0;
    };

    virtual float endSimulation() { return (float) 200; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; };

    /** Get the boundary positions
     *
     * @param i_edge which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(BoundaryEdge i_edge) {
       if ( i_edge == BND_LEFT )
         return (float)0;
       else if ( i_edge == BND_RIGHT)
         return (float)10000;
       else if ( i_edge == BND_BOTTOM )
         return (float)0;
       else
         return (float)10000;
    };
};

/**
 * Scenario "Radial Dam Break":
 * elevated water in the center of the domain
 */
class SWE_RadialDamBreakScenario : public SWE_Scenario {

  public:

    float getBathymetry(float x, float y) {
       return 0;
    };

    float getWaterHeight(float x, float y) {
       return (( sqrt( (x-500.f)*(x-500.f) + (y-500.f)*(y-500.f) ) < 100.f ) ? 150.f: 100.0f)-getBathymetry(x, y);
    };

	virtual float endSimulation() { return (float) 200; };

    virtual BoundaryType getBoundaryType(BoundaryEdge edge) { return OUTFLOW; };

    /** Get the boundary positions
     *
     * @param i_edge which edge
     * @return value in the corresponding dimension
     */
    float getBoundaryPos(BoundaryEdge i_edge) {
       if ( i_edge == BND_LEFT )
         return (float)0;
       else if ( i_edge == BND_RIGHT)
         return (float)1000;
       else if ( i_edge == BND_BOTTOM )
         return (float)0;
       else
         return (float)1000;
    };
};

#endif
