/**
 * @file examples/p8est_mpidynres.cpp
 * 
 * @brief example demonstrating the simulation of three-dimensional
 *        shallow water equations on a p8est with dynamic resource
 *        management using libmpidynres
 * 
 * This file contains the main simulation functionality, including
 * functions for initializing and manipulating the p8est, managing
 * dynamic resources, and the main simulation loop.
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../swe_solvers/src/solver/AugRie.hpp"

#include "scenarios/SWE_Scenario_3d.hh"
#include "scenarios/SWE_simple_scenarios_3d.hh"

#include "writer/p8est_vtk_writer.hh"

#include "tools/args.hh"

#include <../p4est/src/p8est_geometry.h>
#include <../p4est/src/p8est_extended.h>
#include <../p4est/src/p8est_communication.h>

extern "C" {
#include <mpidynres.h>
#include <mpidynres_sim.h>
}

/*
 * @brief state of the simulation
 */
struct simulation_state {
    float t;                        // the current simulation time
    float t_end;                    // the time to be simulated
    float *phases;                  // an array containing all phases in time
    int num_phases;                 // the number of phases in time
    int current_phase;              // the index of the current phase
    char scenario[100];             // the name of the scenario
    char output_name[100];          // the name of the output file written to
    p8est_t *p8est;                 // the p8est simulation object
    MPIDYNRES_Session session;            // the MPI Session object
    char pset[MPI_MAX_PSET_NAME_LEN];   // the curent process set
    sc_MPI_Comm mpicomm;            // the MPI communicator of all working processes
    bool is_new;                    // flag showing whether this process has been 
                                    // running for more than one iteration
};


/**
 * state of a single cell, stored inside every p8est quadrant
 */
struct cell_state {
    //cell state
    double h;
    double hu;
    double hv;
    double hw;

    //derivatives / updates
    double update_h;
    double update_hu;
    double update_hv;
    double update_hw;
};


/**
 * helper struct for passing the ghost layer and the timestep 
 * to the p8est_iterate function simulatanously
 */
struct ghost_and_data {
   cell_state *ghost;
   double timestep;
};


// The riemann solver and the simulated scenario are static variables
// in order to make them accessible from the p8est callback functions below.
// An alternative way whould be storing then in the p8est->user_pointer.
// This is done for the scenario in the benchmarksin src/examples/benchmarks/

static solver::AugRie<double> riemann_solver;
static SWE_Scenario_3d *scenario;


//////////////////////////////////////////////////////////////////////
////////////////////// simulation on the p8est ///////////////////////
//////////////////////////////////////////////////////////////////////


/************************************************
 *      Functions for manipulating the p8est    *
 ***********************************************/

/**
 * @brief Callback function for refining the mesh.
 */
static int refine_fn(p8est_t *p8est, p4est_topidx_t tree, 
        p8est_quadrant_t *quadrant){

    cell_state *state = (cell_state *) quadrant->p.user_data;

    // weighting factor depending on the quadrant's level, 
    // for iterative refining
    double factor = 1;
    factor = pow(2, ((int) quadrant->level - 9));

    return (quadrant->level < 6 || 
            (quadrant->level < 11 && state->update_h / (double) factor > 5));
}

/**
 * @brief Callback function for coarsening the mesh.
 */
static int coarsen_fn(p8est_t *p8est, p4est_topidx_t tree, 
        p8est_quadrant_t *quadrant[]){

    //coarsen, if the maximum h update is too low
    double update_h_max = 0;

    for(int i = 0; i < 8; i++){
        p8est_quadrant_t *quad = quadrant[i];
        cell_state *state = (cell_state *) quad->p.user_data;
        double u_h = state->update_h;
        update_h_max = std::max(update_h_max, std::abs(u_h));
    }

    // weighting factor depending on the quadrant's level, 
    // for iterative coarsening
    double factor = 1;
    factor = pow(2, ((int) quadrant[0]->level - 9));

    return (quadrant[0]->level >= 11 ||
            ((quadrant[0]->level >= 6) && update_h_max / (double) factor < 4));
}

/*
 * @brief Function for initializing cells after adapting the grid
 */
static void replace_fn(p8est_t * p8est, p4est_topidx_t which_tree, 
        int num_outgoing, p8est_quadrant_t * outgoing[],
        int num_incoming, p8est_quadrant_t * incoming[]) {

    if(num_outgoing > 1) {
        //coarsening
        cell_state *new_state = (cell_state *) (incoming[0]->p.user_data);

        new_state->h = 0;
        new_state->hu = 0;
        new_state->hv = 0;
        new_state->hw = 0;

        new_state->update_h = 0;
        new_state->update_hu = 0;
        new_state->update_hv = 0;
        new_state->update_hw = 0;

        //average over all outgoing cells
        for(int i = 0; i < num_outgoing; i++){
            cell_state *old_state = (cell_state *) (outgoing[i]->p.user_data);

            new_state->h += old_state->h;
            new_state->hu += old_state->hu;
            new_state->hv += old_state->hv;
            new_state->hw += old_state->hw;

            new_state->update_h += old_state->update_h;
            new_state->update_hu += old_state->update_hu;
            new_state->update_hv += old_state->update_hv;
            new_state->update_hw += old_state->update_hw;
        }

        new_state->h /= num_outgoing; 
        new_state->hu /= num_outgoing;
        new_state->hv /= num_outgoing;
        new_state->hw /= num_outgoing;

        new_state->update_h /= num_outgoing; 
        new_state->update_hu /= num_outgoing;
        new_state->update_hv /= num_outgoing;
        new_state->update_hw /= num_outgoing;

    } else {
        //refine
        
        //simply copy the state of the old cell to all new cells.
        cell_state *old_state = (cell_state *) (outgoing[0]->p.user_data);

        for(int i = 0; i < num_incoming; i++){
            cell_state *new_state = (cell_state *) (incoming[i]->p.user_data);

            new_state->h = old_state->h;
            new_state->hu = old_state->hu;
            new_state->hv = old_state->hv;
            new_state->hw = old_state->hw;

            new_state->update_h = old_state->update_h;
            new_state->update_hu = old_state->update_hu;
            new_state->update_hv = old_state->update_hv;
            new_state->update_hw = old_state->update_hw;
        }
    }
}

/**
 * @brief initialize a cell using the scenario data
 */
static void init_cell(p8est_t *p8est, p4est_topidx_t tree, 
        p8est_quadrant_t *quadrant){

    if(scenario == NULL)
        return;

    cell_state *state = (cell_state*) quadrant->p.user_data;
    float x = quadrant->x;
    float y = quadrant->y;
    float z = quadrant->z;

    double xyz[3];
    p8est_qcoord_to_vertex(p8est->connectivity, tree, x, y, z, xyz); 

    xyz[0] *= scenario->getBoundaryPos(BND_RIGHT) - 
                  scenario->getBoundaryPos(BND_LEFT);
    xyz[1] *= scenario->getBoundaryPos(BND_FRONT) - 
                  scenario->getBoundaryPos(BND_BACK);
    xyz[2] *= scenario->getBoundaryPos(BND_TOP) - 
                  scenario->getBoundaryPos(BND_BOTTOM);
    xyz[0] += scenario->getBoundaryPos(BND_LEFT);
    xyz[1] += scenario->getBoundaryPos(BND_BACK);
    xyz[2] += scenario->getBoundaryPos(BND_BOTTOM);

    state->h = scenario->getWaterHeight(xyz[0], xyz[1], xyz[2]);
    state->hu = scenario->getVeloc_u(xyz[0], xyz[1], xyz[2]);
    state->hv = scenario->getVeloc_v(xyz[0], xyz[1], xyz[2]);
    state->hw = scenario->getVeloc_w(xyz[0], xyz[1], xyz[2]);

    state->update_h = 0;
    state->update_hu = 0;
    state->update_hv = 0;
    state->update_hw = 0;
}



/*****************************************
 *       Functions for simulation        *
 ****************************************/


/**
 * @brief reset all updates from previous timesteps
 */
static void resetUpdates(p8est_iter_volume_info_t * info, void *user_data){

    p8est_quadrant_t *quad = (p8est_quadrant_t *) info->quad;
    cell_state *state = (cell_state *) quad->p.user_data;
    
    state->update_h = 0;
    state->update_hu = 0;
    state->update_hv = 0;
    state->update_hw = 0;
}


/**
 * @brief compute updates for one face
 *
 * This function uses the AugRie solver to calculate updates 
 * at all faces of the tree.
 * For simplicity, the resulting updates are stored inside the quadrants.
 */
static void computeNumericalFluxes(p8est_iter_face_info_t * info, 
        void *user_data){

    cell_state *ghost_data = ((ghost_and_data *) user_data)->ghost;

    sc_array_t *sides = &(info->sides);
    p8est_iter_face_side_t *first_side = 
            p8est_iter_fside_array_index_int (sides, 0);
    p8est_iter_face_side_t *second_side = 
            p8est_iter_fside_array_index_int (sides, 1);

    cell_state *first;
    cell_state *second;

    //handle hanging faces
    if(first_side->is_hanging) {
        first = new cell_state {0, 0, 0, 0, 0, 0, 0};
        
        //average over small cells
        for(int i = 0; i < P8EST_HALF; i++){
            cell_state *small_cell;

            //handle ghost layer
            if(first_side->is.hanging.is_ghost[i]) {
                small_cell = (cell_state *) &ghost_data[first_side->is.hanging.quadid[i]];
            } else {
                small_cell = (cell_state *) first_side->is.hanging.quad[i]->p.user_data;
            }

            first->h += small_cell->h;
            first->hu += small_cell->hu;
            first->hv += small_cell->hv;
            first->hw += small_cell->hw;
        }

        first->h /= P8EST_HALF;
        first->hu /= P8EST_HALF;
        first->hv /= P8EST_HALF;
        first->hw /= P8EST_HALF;

    } else {
        //handle ghost layer
        if(first_side->is.full.is_ghost) {
            first = (cell_state *) &ghost_data[first_side->is.full.quadid];
        } else {
            first = (cell_state *) first_side->is.full.quad->p.user_data;
        }
    }

    //handle hanging faces
    if(second_side->is_hanging) {
        second = new cell_state {0, 0, 0, 0, 0, 0, 0};
        
        //average over small cells
        for(int i = 0; i < P8EST_HALF; i++){
            cell_state *small_cell;

            //handle ghost layer
            if(second_side->is.hanging.is_ghost[i]) {
                small_cell = (cell_state *) &ghost_data[second_side->is.hanging.quadid[i]];
            } else {
                small_cell = (cell_state *) second_side->is.hanging.quad[i]->p.user_data;
            }

            second->h += small_cell->h;
            second->hu += small_cell->hu;
            second->hv += small_cell->hv;
            second->hw += small_cell->hw;
        }

        second->h /= P8EST_HALF;
        second->hu /= P8EST_HALF;
        second->hv /= P8EST_HALF;
        second->hw /= P8EST_HALF;

    } else {
        //handle ghost layer
        if(second_side->is.full.is_ghost) {
            second = (cell_state *) &ghost_data[second_side->is.full.quadid];
        } else {
            second = (cell_state *) second_side->is.full.quad->p.user_data;
        }
    }

    //ordering of the quadrants
    if(first_side->face % 2 != 0){
        std::swap(first, second);
        std::swap(first_side, second_side);
    }


    //calculate maximum timestep from the minimal cell width

    double *maxTimestep = &((ghost_and_data *) user_data)->timestep;
    int lowestLevel;
    if(first_side->is_hanging){
        if(second_side->is_hanging){
            lowestLevel = std::max(first_side->is.hanging.quad[0]->level, 
                    second_side->is.hanging.quad[0]->level);
        } else {
            lowestLevel = std::max(first_side->is.hanging.quad[0]->level, 
                    second_side->is.full.quad->level);
        }
    } else {
        if(second_side->is_hanging){
            lowestLevel = std::max(first_side->is.full.quad->level, 
                    second_side->is.hanging.quad[0]->level);
        } else {
            lowestLevel = std::max(first_side->is.full.quad->level, 
                    second_side->is.full.quad->level);
        }
    }


    // compute updates
    double update_first_h = 0;
    double update_second_h = 0;
    double current_maxWaveSpeed;

    if(first_side->face < 2){
        //x-direction

        double update_first_hu = 0;
        double update_second_hu = 0;

        riemann_solver.computeNetUpdates(first->h, second->h,
                                 first->hu, second->hu,
                                 0, 0,
                                 update_first_h, update_second_h,
                                 update_first_hu, update_second_hu,
                                 current_maxWaveSpeed);

        first->update_h += update_first_h;
        second->update_h += update_second_h;
        first->update_hu += update_first_hu;
        second->update_hu += update_second_hu;

        double min_size = P8EST_QUADRANT_LEN (lowestLevel);
        min_size *= (scenario->getBoundaryPos(BND_RIGHT) - 
                        scenario->getBoundaryPos(BND_LEFT)) / 
                    (double) P8EST_ROOT_LEN;
        double current_maxTimestep = 0.4 * min_size / current_maxWaveSpeed;

        if(current_maxTimestep < *maxTimestep)
            *maxTimestep = current_maxTimestep;

        // if one face is hanging, we have to write the updates back 
        // to all of the small cells
        if(first_side->is_hanging){
            for(int i = 0; i < P8EST_HALF; i++){
                cell_state *small_cell;

                //handle ghost cells
                if(first_side->is.hanging.is_ghost[i]) {
                    small_cell = (cell_state *) &ghost_data[first_side->is.hanging.quadid[i]];
                } else {
                    small_cell = (cell_state *) first_side->is.hanging.quad[i]->p.user_data;
                }

                //write data
                small_cell->update_h += update_first_h;
                small_cell->update_hu += update_first_hu;
            }
        }

        if(second_side->is_hanging){
            for(int i = 0; i < P8EST_HALF; i++){
                cell_state *small_cell;

                //handle ghost cells
                if(second_side->is.hanging.is_ghost[i]) {
                    small_cell = (cell_state *) &ghost_data[second_side->is.hanging.quadid[i]];
                } else {
                    small_cell = (cell_state *) second_side->is.hanging.quad[i]->p.user_data;
                }

                //write data
                small_cell->update_h += update_second_h;
                small_cell->update_hu += update_second_hu;
            }
        }

    } else if(first_side->face < 4) {
        //y-direction

        double update_first_hv = 0;
        double update_second_hv = 0;

        riemann_solver.computeNetUpdates(first->h, second->h,
                                 first->hv, second->hv,
                                 0, 0,
                                 update_first_h, update_second_h,
                                 update_first_hv, update_second_hv,
                                 current_maxWaveSpeed);

        first->update_h += update_first_h;
        second->update_h += update_second_h;
        first->update_hv += update_first_hv;
        second->update_hv += update_second_hv;


        double min_size = P8EST_QUADRANT_LEN (lowestLevel);
        min_size *= (scenario->getBoundaryPos(BND_TOP) - 
                            scenario->getBoundaryPos(BND_BOTTOM)) /
                        (double) P8EST_ROOT_LEN;
        double current_maxTimestep = 0.4 * min_size / current_maxWaveSpeed;

        if(current_maxTimestep < *maxTimestep)
            *maxTimestep = current_maxTimestep;

        // if one face is hanging, we have to write the updates back to 
        // the small cells
        if(first_side->is_hanging){
            for(int i = 0; i < P8EST_HALF; i++){
                cell_state *small_cell;

                //handle ghost cells
                if(first_side->is.hanging.is_ghost[i]) {
                    small_cell = (cell_state *) &ghost_data[first_side->is.hanging.quadid[i]];
                } else {
                    small_cell = (cell_state *) first_side->is.hanging.quad[i]->p.user_data;
                }

                //write data
                small_cell->update_h += update_first_h;
                small_cell->update_hv += update_first_hv;
            }
        }

        if(second_side->is_hanging){
            for(int i = 0; i < P8EST_HALF; i++){
                cell_state *small_cell;

                //handle ghost cells
                if(second_side->is.hanging.is_ghost[i]) {
                    small_cell = (cell_state *) &ghost_data[second_side->is.hanging.quadid[i]];
                } else {
                    small_cell = (cell_state *) second_side->is.hanging.quad[i]->p.user_data;
                }

                //write data
                small_cell->update_h += update_second_h;
                small_cell->update_hv += update_second_hv;
            }
        }
    } else {
        //z-direction

        double update_first_hw = 0;
        double update_second_hw = 0;

        riemann_solver.computeNetUpdates(first->h, second->h,
                                 first->hw, second->hw,
                                 0, 0,
                                 update_first_h, update_second_h,
                                 update_first_hw, update_second_hw,
                                 current_maxWaveSpeed);

        first->update_h += update_first_h;
        second->update_h += update_second_h;
        first->update_hw += update_first_hw;
        second->update_hw += update_second_hw;

        double min_size = P8EST_QUADRANT_LEN (lowestLevel);
        min_size *= scenario->getBoundaryPos(BND_RIGHT) - 
                        scenario->getBoundaryPos(BND_LEFT) / 
                        (double) P8EST_ROOT_LEN;
        double current_maxTimestep = 0.4 * min_size / current_maxWaveSpeed;

        if(current_maxTimestep < *maxTimestep)
            *maxTimestep = current_maxTimestep;

        //if one face is hanging, we have to write the updates back to the small cells
        if(first_side->is_hanging){
            for(int i = 0; i < P8EST_HALF; i++){
                cell_state *small_cell;

                //handle ghost cells
                if(first_side->is.hanging.is_ghost[i]) {
                    small_cell = (cell_state *) &ghost_data[first_side->is.hanging.quadid[i]];
                } else {
                    small_cell = (cell_state *) first_side->is.hanging.quad[i]->p.user_data;
                }

                //write data
                small_cell->update_h += update_first_h;
                small_cell->update_hw += update_first_hw;
            }
        }

        if(second_side->is_hanging){
            for(int i = 0; i < P8EST_HALF; i++){
                cell_state *small_cell;

                //handle ghost cells
                if(second_side->is.hanging.is_ghost[i]) {
                    small_cell = (cell_state *) &ghost_data[second_side->is.hanging.quadid[i]];
                } else {
                    small_cell = (cell_state *) second_side->is.hanging.quad[i]->p.user_data;
                }

                //write data
                small_cell->update_h += update_second_h;
                small_cell->update_hw += update_second_hw;
            }
        }
    }

}

/**
 * @brief apply the updates previously calculated, 
 *        and reset the updates for the next timestep
 */
static void updateUnknowns(p8est_iter_volume_info_t * info, void *user_data){

    double *maxTimestep = (double *) user_data;
    p8est_quadrant_t *quad = (p8est_quadrant_t *) info->quad;
    cell_state *state = (cell_state *) quad->p.user_data;
    double dx = P8EST_QUADRANT_LEN (quad->level);
    dx *= (scenario->getBoundaryPos(BND_TOP) - 
                   scenario->getBoundaryPos(BND_BOTTOM)) / 
               P8EST_ROOT_LEN;

    //apply stored updates
    state->h -= *maxTimestep / dx  * state->update_h;
    state->hu -= *maxTimestep / dx * state->update_hu;
    state->hv -= *maxTimestep / dx * state->update_hv;
    state->hw -= *maxTimestep / dx * state->update_hw;
}



/**
 * @brief initializes the p8est and the initial state of all 
 *        grid cells according to the given SWE_Scenario.
 * 
 * @param [in] mpicomm     MPI communicator the p8est should be running on
 * @param [in] scenario  scenario, which is used during the setup.
 *
 * @return the p8est structure
 */
p8est_t *initP8est(sc_MPI_Comm mpicomm, SWE_Scenario_3d *scenario) {

    //create p8est
    p8est_connectivity_t *conn = p8est_connectivity_new_periodic();

    p8est_t *p8est = p8est_new(mpicomm, conn, sizeof(cell_state), 
            init_cell, NULL);

    //refine tree to minimum level
    for(int level = 0; level < 5; level++) {
        p8est_refine(p8est, 0, refine_fn, init_cell);
        p8est_balance(p8est, P8EST_CONNECT_FACE, init_cell);
        p8est_partition(p8est, 0, NULL);
    }

    //calculate one iteration to decide on further refinement,
    //refine by maximum three levels
    for(int i = 0; i < 3; i++) {
        p8est_ghost_t *ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);
        cell_state *ghost_data = P4EST_ALLOC(cell_state, 
                ghost->ghosts.elem_count);
        p8est_ghost_exchange_data(p8est, ghost, ghost_data);
        double maxTimestep = 100;
        ghost_and_data *param = new ghost_and_data {ghost_data, maxTimestep};
        
        p8est_iterate (p8est, ghost,
                (void *) param,
                NULL,
                computeNumericalFluxes,
                NULL, NULL);

        p8est_refine(p8est, 0, refine_fn, init_cell);
        p8est_balance(p8est, P8EST_CONNECT_FACE, init_cell);
        p8est_partition(p8est, 0, NULL);

        p8est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);
        ghost = NULL;
        ghost_data = NULL;
    }

    return p8est;
}

/**
 * @brief main simulation loop between two phases 
 *
 * @param state    the current simulation state
 *
 * @return the actual end time reached
 */
float simulate_interval(simulation_state *state) {

    p8est_t *p8est = state->p8est;
    float t = state->t;

    // repartition the p8est
    p8est_partition(p8est, 0, NULL);

    // create the ghost quadrants
    p8est_ghost_t *ghost = p8est_ghost_new(p8est, P8EST_CONNECT_FULL);

    // create space for storing the ghost data
    cell_state *ghost_data = P4EST_ALLOC(cell_state, 
            ghost->ghosts.elem_count);

    do {
        // reset updates
        p8est_iterate (p8est, NULL,
                NULL,
                resetUpdates,
                NULL,
                NULL, NULL);

        // synchronize the ghost data
        p8est_ghost_exchange_data(p8est, ghost, ghost_data);

        // reset maxTimestep
        double maxTimestep = 100;

        ghost_and_data *param = new ghost_and_data {ghost_data, maxTimestep};

        // compute numerical fluxes for every edge
        // updates maxTimestep
        p8est_iterate (p8est, ghost,
                (void *) param,
                NULL,
                computeNumericalFluxes,
                NULL, NULL);
        
        //collect the timestep size from all processes
        int mpiret = sc_MPI_Allreduce (&(param->timestep), &maxTimestep, 1, 
                sc_MPI_DOUBLE, sc_MPI_MIN, state->mpicomm);
        SC_CHECK_MPI (mpiret);

        // update unknowns accordingly
        p8est_iterate (p8est, NULL,
                (void *) &maxTimestep,
                updateUnknowns,
                NULL,
                NULL, NULL);

        delete param;
        t += maxTimestep;

    } while (t <= state->phases[state->current_phase]);

    p8est_ghost_destroy(ghost);
    P4EST_FREE(ghost_data);
    ghost = NULL;
    ghost_data = NULL;

    // adapt the grid
    for(int i = 0; i < 2; i++) {
        p8est_refine_ext(p8est, 0, 10, refine_fn, NULL, replace_fn);
    }
    for(int i = 0; i < 5; i++) {
        p8est_coarsen_ext(p8est, 0, 0, coarsen_fn, NULL, replace_fn);
    }

    p8est_balance_ext(p8est, P8EST_CONNECT_FACE, NULL, replace_fn);

    return t;
}



////////////////////////////////////////////////////////////////
//////////// Functions for ressource management ////////////////
////////////////////////////////////////////////////////////////


/*
 * @brief create an MPI_Info object containing the current 
 *        simulation state
 *
 * @param state the current simulation state
 *
 * @return an MPI_Info object containing the simulation state information
 */
MPI_Info create_info_from_state(simulation_state *state) {
   
    MPI_Info info;
    MPI_Info_create(&info);

    MPI_Info_set(info, "output_name", state->output_name);

    char t_end[100];
    snprintf(t_end, 99, "%f", state->t_end);
    MPI_Info_set(info, "end_time", t_end);

    char t[100];
    snprintf(t, 99, "%f", state->t);
    MPI_Info_set(info, "current_time", t);

    char num_phases[100];
    snprintf(num_phases, 99, "%d", state->num_phases);
    MPI_Info_set(info, "num_phases", num_phases);

    char current_phase[100];
    snprintf(current_phase, 99, "%d", state->current_phase);
    MPI_Info_set(info, "current_phase", current_phase);

    MPI_Info_set(info, "pset", state->pset);

    MPI_Info_set(info, "scenario", state->scenario);
    
    return info;
}


bool check_resource_change(simulation_state *state) {

    bool must_terminate = false;

    int mpirank = -1;
    if(state->p8est != NULL)
        mpirank = state->p8est->mpirank;

    int working_size = state->p8est->mpisize;
    int idle_size = -1;
    MPI_Comm_size(MPI_COMM_WORLD, &idle_size);
    idle_size -= working_size + 1;

    MPI_Info rc_info = MPI_INFO_NULL;
    MPIDYNRES_RC_type rc_type;
    char delta_pset[MPI_MAX_PSET_NAME_LEN];
    int rc_tag;

    if(mpirank == 0) {
        // decide on size and kind of resource change
        // using an optimal number of 10000 quadrants per core
        int optimal_num_cores = 
            state->p8est->global_num_quadrants / 10000 + 1;
        int delta = optimal_num_cores - state->p8est->mpisize;
        
        if(delta == 0) {
            rc_type = MPIDYNRES_RC_NONE;
        } else if(delta > 0) {
            rc_type = MPIDYNRES_RC_ADD;
            delta = std::min(delta, idle_size);
        } else {
            rc_type = MPIDYNRES_RC_SUB;
            delta = std::max(delta, -working_size + 1);
        }


        // request a resource change
        if(rc_type != MPIDYNRES_RC_NONE)
            DEBUG_MPIDYNRES_RC_user_request(state->session, delta, 
                    delta_pset, &rc_type, &rc_tag, &rc_info);
        // ignore rc_info
        if(rc_info != MPI_INFO_NULL)
            MPI_Info_free(&rc_info);
    }

    // bcast changes to all processes
    MPI_Bcast(&rc_type, 1, MPI_INT, 0, state->mpicomm);
    MPI_Bcast(&rc_tag, 1, MPI_INT, 0, state->mpicomm);
    MPI_Bcast(&delta_pset, MPI_MAX_PSET_NAME_LEN, MPI_CHAR, 0, state->mpicomm);


    // handle changes
    if(rc_type != MPIDYNRES_RC_NONE) {
        switch(rc_type) {
            case MPIDYNRES_RC_ADD:
                //create new pset
                if(mpirank == 0) 
                    MPIDYNRES_pset_create_op(state->session, MPI_INFO_NULL, 
                            state->pset, delta_pset, MPIDYNRES_PSET_UNION, 
                            state->pset);
                //send pset to all processes
                MPI_Bcast(state->pset, MPI_MAX_PSET_NAME_LEN, 
                        MPI_CHAR, 0, state->mpicomm);
                break;
            case MPIDYNRES_RC_SUB:
                //create new pset
                if(mpirank == 0) 
                    MPIDYNRES_pset_create_op(state->session, MPI_INFO_NULL, 
                            state->pset, delta_pset, MPIDYNRES_PSET_DIFFERENCE,
                            state->pset);
                //send pset to all processes
                MPI_Bcast(state->pset, MPI_MAX_PSET_NAME_LEN, 
                        MPI_CHAR, 0, state->mpicomm);
                //check whether this process has to terminate
                MPI_Info own_psets;
                MPIDYNRES_Session_get_psets(state->session, MPI_INFO_NULL, 
                        &own_psets);
                int contains_key, valuelen;
                MPI_Info_get_valuelen(own_psets, state->pset, 
                        &valuelen, &contains_key);
                if(!contains_key)
                    must_terminate = true;
                break;
            default:
                break;
        }


        // accept changes
        if(mpirank == 0) {
            MPI_Info state_info = create_info_from_state(state);
            MPIDYNRES_RC_accept(state->session, rc_tag, state_info);
        }

        // update communicator
        MPI_Group group;
        MPIDYNRES_Group_from_session_pset(state->session, state->pset, &group);
        if(!must_terminate) {
            MPIDYNRES_Comm_create_from_group(group, NULL, MPI_INFO_NULL,
                MPI_ERRORS_ARE_FATAL, &state->mpicomm);
        } else {
            state->mpicomm = MPI_COMM_NULL;
        }
        MPI_Group_free(&group);
    }
    return rc_type != MPIDYNRES_RC_NONE;
}




/////////////////////////////////////////////////////////////////
////////////////////// Main functionality ///////////////////////
/////////////////////////////////////////////////////////////////


/**
 * @brief initialize the simulation
 *
 * Initializes the simulation using cli parameters for the first
 * processes and an MPI_Info object for dynamically added ones.
 * Initializes the MPIDYNRES_Session environment, the p8est, the scenario,
 * and the simulation state.
 *
 * @param argc, argv    the cli parameters for this application, 
 *                      passed on by the main function and ignored
 *                      by processes allocated during runtime.
 *
 * @return the simulation state.
 */
simulation_state initialize(int argc, char* argv[]) {

    simulation_state state;

    int contains_key;
    int value_len;
    MPI_Info psets;

    //init MPI Session
    MPIDYNRES_Session_init(MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &state.session);
    MPIDYNRES_Session_get_psets(state.session, MPI_INFO_NULL, &psets);


    //check whether this process was spawned at startup time
    MPI_Info_get_valuelen(psets, "mpi://WORLD", &value_len, &contains_key);
    if(contains_key) {

        //initialize using cli parameters
        tools::Args args;

        args.addOption("scenario", 's', 
                "Load scenario, (artificial | dambreak), default: dambreak",
                args.Required, false); 
        args.addOption("output-basepath", 'o', 
                "Output base file name, mandatory argument, no default.");
        args.addOption("simulated-time", 't', 
                "Time in seconds to be simulated, default: 200", 
                args.Required, false);
        args.addOption("num-phases", 'c', 
                "Number of phases for visualization (at each phase\
                    in time, an output file is written), default: 200", 
                args.Required, false);

        tools::Args::Result ret = args.parse(argc, argv);

        switch (ret) {
            case tools::Args::Error:
                exit(1);
            case tools::Args::Help:
                exit(0);
            default:
                break; 
        }

        // read command line parameters
        strcpy(state.output_name, 
            args.getArgument<std::string>("output-basepath").c_str());
        state.num_phases = args.getArgument<int>("num-phases", 200);
        state.t_end = args.getArgument<float>("simulated-time", 200);

        std::string scenario_option = "";

        scenario_option = 
            args.getArgument<std::string>("scenario", "dambreak");


        strcpy(state.scenario, scenario_option.c_str());
        state.t = 0;
        state.current_phase = 0;
        strcpy(state.pset, "mpi://WORLD");

    } else {

        // initialize using the session info object
        MPI_Info session_info;
        MPIDYNRES_Session_get_info(state.session, &session_info);

        int exists;
        MPI_Info_get(session_info, "output_name", 99, state.output_name, &exists);
        if(!exists)
            std::cout << "ERROR output\n";

        char *num_phases = new char[100];
        MPI_Info_get(session_info, "num_phases", 99, 
                num_phases, &exists);
        if(!exists)
            std::cout << "ERROR num\n";
        state.num_phases = atoi(num_phases);

        char *t_end = new char[100];
        MPI_Info_get(session_info, "end_time", 99, t_end, &exists);
        if(!exists)
            std::cout << "ERROR end\n";
        state.t_end = atof(t_end);

        char *t = new char[100];
        MPI_Info_get(session_info, "current_time", 99, t, &exists);
        if(!exists)
            std::cout << "ERROR t\n";
        state.t = atof(t);

        char *current_phase = new char[100];
        MPI_Info_get(session_info, "current_phase", 99, 
                current_phase, &exists);
        if(!exists)
            std::cout << "ERROR phase\n";
        state.current_phase = atoi(current_phase);

        MPI_Info_get(session_info, "pset", MPI_MAX_PSET_NAME_LEN - 1, 
                state.pset, &exists);
        if(!exists)
            std::cout << "ERROR pset\n";

        MPI_Info_get(session_info, "scenario", 99, state.scenario, &exists);
    }

    state.phases = new float[state.num_phases + 1];
    for (int cp = 0; cp <= state.num_phases; cp++) {
        state.phases[cp] = cp * (state.t_end / state.num_phases);
    }

    state.is_new = true;

    //initialize scenario

    scenario = new SWE_RadialDamBreakScenario();


    //create mpicomm
    int mpiret;
    MPI_Group group;
    mpiret = MPIDYNRES_Group_from_session_pset(state.session, state.pset, &group);
    if(mpiret)
        MPIDYNRES_exit();

    mpiret = MPIDYNRES_Comm_create_from_group(group, NULL, MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &state.mpicomm);
    if(mpiret)
        MPIDYNRES_exit();

    // create and init the p8est simulation object
    state.p8est = NULL;
    if(contains_key) {
        state.p8est = initP8est(state.mpicomm, scenario);
        state.is_new = false;
    }

    MPI_Group_free(&group);

    return state;
}


/**
 * @brief free all allocated structures and destroy the p8est
 */
void finalize(simulation_state *state) {

    if(state->p8est != NULL) {
        p8est_connectivity_t *connectivity = state->p8est->connectivity;
        p8est_destroy(state->p8est);
        p8est_connectivity_destroy(connectivity);
    }
    delete[] state->phases;
}


/**
 * @brief main work loop
 *
 * calls simulate_interval and adapts the computing resources
 * writes simulation output to vtk files
 *
 * @param state    the state of the simulation
 *
 * @return non-zero if errors occurred
 */
int simulate(simulation_state *state){

    io::P8est_vtkWriter *l_writer = 
        new io::P8est_vtkWriter(state->output_name);

    // write initial state
    if(state->current_phase == 0) {
        l_writer->writeTimeStep(state->p8est);
        state->current_phase = 1;
    }
    l_writer->set_timestep(state->current_phase);

    unsigned int l_iterations = 0;

    // main simulation loop over phases
    for(; state->current_phase <= state->num_phases && 
                state->t <= state->t_end; 
            state->current_phase++) {

        // handle dynamic resources
        bool changed = state->is_new;
        if(!state->is_new) 
            changed = check_resource_change(state);
        if(changed) {
            if(state->is_new)
                state->p8est = NULL;
            state->p8est = p8est_dynres_replace(state->p8est, 
                    state->mpicomm);
        }

        if(state->p8est == NULL)
            return 0;

        // simulate until next timestep
        float new_time = simulate_interval(state);

        // update simulation time with time step width.
        state->t = new_time;
        l_iterations++;

        // write current simulation state
        l_writer->writeTimeStep(state->p8est);

        state->is_new = false;
    }

    //Finalize
    delete l_writer;
    return 0;
}


/**
 * @brief main function serving as init point for processes inside the
 *        MPIDYNRES_Sessions simulation layer
 */
int MPIDYNRES_main(int argc, char **argv) {

    simulation_state state = initialize(argc, argv);

    simulate(&state);

    return 0;
}


/**
 * @brief main program
 *
 * Starts the libmpidynres simulation layer, 
 * initializes new processes and starts the simulation
 */
int main(int argc, char **argv) {

    //initialize MPI
    sc_MPI_Comm mpicomm;
    int mpiret = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI(mpiret);
    mpicomm=sc_MPI_COMM_WORLD;
    int mpisize;
    MPI_Comm_size(mpicomm, &mpisize);

    //initialize mpidynres
    MPI_Info manager_config;
    MPI_Info_create(&manager_config);
    char buf[0x20];
    snprintf(buf, 0x20, "%d", 2);
    MPI_Info_set(manager_config, "manager_initial_number", buf);

    MPIDYNRES_SIM_config mpidynres_config = {
        .base_communicator = MPI_COMM_WORLD,
        .manager_config = manager_config,
    };

    MPIDYNRES_SIM_start(mpidynres_config, argc, argv, MPIDYNRES_main);

    mpiret = sc_MPI_Finalize();
    SC_CHECK_MPI (mpiret);

    return 0;
}

