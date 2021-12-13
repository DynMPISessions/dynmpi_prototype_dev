/**
 * @file examples/p4est_dynres.cpp
 * 
 * @brief example demonstrating the simulation of two-dimensional
 *        shallow water equations on a p4est with dynamic resource
 *        management by running the simulation on different subsets
 *        of a fixed number of processes.
 * 
 * This file contains the main simulation functionality, including
 * functions for initializing and manipulating the p4est, managing
 * dynamic resources, and the main simulation loop.
 *
 * The resource adaptation in this example does not react to the 
 * current workload, but steadily increases the number of processes.
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>

#include "../swe_solvers/src/solver/AugRie.hpp"

#include "scenarios/SWE_Scenario.hh"
#include "scenarios/SWE_simple_scenarios.hh"

#include "writer/p4est_vtk_writer.hh"

#include "tools/args.hh"

#include <../p4est/src/p4est_geometry.h>
#include <../p4est/src/p4est_extended.h>
#include <../p4est/src/p4est_communication.h>

/*
 * @brief state of the simulation
 */
struct simulation_state {
    float t;                        // the current simulation time
    float t_end;                    // the time to be simulated
    float *phases;             // an array containing all phases in time
    int num_phases;            // the number of phases in time
    int current_phase;         // the index of the current phase
    char scenario[100];             // the simulated scenario
    char bathymetry[100];           // the path to a file with bathymetry data
    char displacement[100];         // the path to a file with displacement data
    char output_name[100];          // the name of the output file written to
    p4est_t *p4est;                 // the p4est simulation object
    sc_MPI_Comm mpicomm;            // the MPI communicator
    sc_MPI_Group active;            // the group of all active processes
    bool is_new;                    // flag showing whether this process has been 
                                    // running for more than one iteration
};


/**
 * state of a single cell, stored inside every p4est quadrant
 */
struct cell_state {
    //cell state
    double h;
    double hu;
    double hv;
    double b;

    //derivatives / updates
    double update_h;
    double update_hu;
    double update_hv;
};


/**
 * helper struct for passing the ghost layer and the timestep 
 * to the p4est_iterate function simulatanously
 */
struct ghost_and_data {
   cell_state *ghost;
   double timestep;
};


// The riemann solver and the simulated scenario are static variables
// in order to make them accessible from the p4est callback functions below.
// An alternative way whould be storing then in the p4est->user_pointer.
// This is done for the scenario in the benchmarksin src/examples/benchmarks/

static solver::AugRie<double> riemann_solver;
static SWE_Scenario *scenario;

//////////////////////////////////////////////////////////////////////
////////////////////// simulation on the p4est ///////////////////////
//////////////////////////////////////////////////////////////////////


/************************************************
 *      Functions for manipulating the p4est    *
 ***********************************************/

/**
 * @brief Callback function for refining the mesh.
 */
static int refine_fn(p4est_t *p4est, p4est_topidx_t tree, 
        p4est_quadrant_t *quadrant){

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
static int coarsen_fn(p4est_t *p4est, p4est_topidx_t tree, 
        p4est_quadrant_t *quadrant[]){

    //coarsen, if the maximum h update is too low
    double update_h_max = 0;

    for(int i = 0; i < 4; i++){
        p4est_quadrant_t *quad = quadrant[i];
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
 * @brief Function for initializing cels after adapting the grid
 */
static void replace_fn(p4est_t * p4est, p4est_topidx_t which_tree, 
        int num_outgoing, p4est_quadrant_t * outgoing[],
        int num_incoming, p4est_quadrant_t * incoming[]) {

    if(num_outgoing > 1) {
        //coarsening
        cell_state *new_state = (cell_state *) (incoming[0]->p.user_data);

        new_state->h = 0;
        new_state->hu = 0;
        new_state->hv = 0;
        new_state->b = 0;

        new_state->update_h = 0;
        new_state->update_hu = 0;
        new_state->update_hv = 0;

        //average over all outgoing cells
        for(int i = 0; i < num_outgoing; i++){
            cell_state *old_state = (cell_state *) (outgoing[i]->p.user_data);

            new_state->h += old_state->h;
            new_state->hu += old_state->hu;
            new_state->hv += old_state->hv;
            new_state->b += old_state->b;

            new_state->update_h += old_state->update_h;
            new_state->update_hu += old_state->update_hu;
            new_state->update_hv += old_state->update_hv;
        }

        new_state->h /= num_outgoing; 
        new_state->hu /= num_outgoing;
        new_state->hv /= num_outgoing;
        new_state->b /= num_outgoing; 

        new_state->update_h /= num_outgoing; 
        new_state->update_hu /= num_outgoing;
        new_state->update_hv /= num_outgoing;

    } else {
        //refine
        
        //simply copy the state of the old cell to all new cells.
        cell_state *old_state = (cell_state *) (outgoing[0]->p.user_data);

        for(int i = 0; i < num_incoming; i++){
            cell_state *new_state = (cell_state *) (incoming[i]->p.user_data);

            new_state->h = old_state->h;
            new_state->hu = old_state->hu;
            new_state->hv = old_state->hv;
            new_state->b = old_state->b;

            new_state->update_h = old_state->update_h;
            new_state->update_hu = old_state->update_hu;
            new_state->update_hv = old_state->update_hv;
        }
    }
}

/**
 * @brief initialize a cell using the scenario data
 */
static void init_cell(p4est_t *p4est, p4est_topidx_t tree, 
        p4est_quadrant_t *quadrant){

    if(scenario == NULL)
        return;

    cell_state *state = (cell_state*) quadrant->p.user_data;
    float x = quadrant->x;
    float y = quadrant->y;

    double xyz[3];
    p4est_qcoord_to_vertex(p4est->connectivity, tree, x, y, xyz); 

    xyz[0] *= scenario->getBoundaryPos(BND_RIGHT) - 
                  scenario->getBoundaryPos(BND_LEFT);
    xyz[1] *= scenario->getBoundaryPos(BND_TOP) - 
                  scenario->getBoundaryPos(BND_BOTTOM);
    xyz[0] += scenario->getBoundaryPos(BND_LEFT);
    xyz[1] += scenario->getBoundaryPos(BND_BOTTOM);

    state->h = scenario->getWaterHeight(xyz[0], xyz[1]);
    state->hu = scenario->getVeloc_u(xyz[0], xyz[1]);
    state->hv = scenario->getVeloc_v(xyz[0], xyz[1]);
    state->b = scenario->getBathymetry(xyz[0], xyz[1]);

    state->update_h = 0;
    state->update_hu = 0;
    state->update_hv = 0;
}



/*****************************************
 *       Functions for simulation        *
 ****************************************/


/**
 * @brief reset all updates from previous timesteps
 */
static void resetUpdates(p4est_iter_volume_info_t * info, void *user_data){

    p4est_quadrant_t *quad = (p4est_quadrant_t *) info->quad;
    cell_state *state = (cell_state *) quad->p.user_data;
    
    state->update_h = 0;
    state->update_hu = 0;
    state->update_hv = 0;
}


/**
 * @brief compute updates for one face
 *
 * This function uses the AugRie solver to calculate updates 
 * at all faces of the tree.
 * For simplicity, the resulting updates are stored inside the quadrants.
 */
static void computeNumericalFluxes(p4est_iter_face_info_t * info, 
        void *user_data){

    cell_state *ghost_data = ((ghost_and_data *) user_data)->ghost;

    sc_array_t *sides = &(info->sides);
    p4est_iter_face_side_t *first_side = 
            p4est_iter_fside_array_index_int (sides, 0);
    p4est_iter_face_side_t *second_side = 
            p4est_iter_fside_array_index_int (sides, 1);

    cell_state *first;
    cell_state *second;

    //handle hanging faces
    if(first_side->is_hanging) {
        first = new cell_state {0, 0, 0, 0, 0, 0, 0};
        
        //average over small cells
        for(int i = 0; i < P4EST_HALF; i++){
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
            first->b += small_cell->b;
        }

        first->h /= P4EST_HALF;
        first->hu /= P4EST_HALF;
        first->hv /= P4EST_HALF;
        first->b /= P4EST_HALF;

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
        for(int i = 0; i < P4EST_HALF; i++){
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
            second->b += small_cell->b;
        }

        second->h /= P4EST_HALF;
        second->hu /= P4EST_HALF;
        second->hv /= P4EST_HALF;
        second->b /= P4EST_HALF;

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
                                 first->b, second->b,
                                 update_first_h, update_second_h,
                                 update_first_hu, update_second_hu,
                                 current_maxWaveSpeed);

        first->update_h += update_first_h;
        second->update_h += update_second_h;
        first->update_hu += update_first_hu;
        second->update_hu += update_second_hu;

        double min_size = P4EST_QUADRANT_LEN (lowestLevel);
        min_size *= (scenario->getBoundaryPos(BND_RIGHT) - 
                        scenario->getBoundaryPos(BND_LEFT)) / 
                    (double) P4EST_ROOT_LEN;
        double current_maxTimestep = 0.4 * min_size / current_maxWaveSpeed;

        if(current_maxTimestep < *maxTimestep)
            *maxTimestep = current_maxTimestep;

        // if one face is hanging, we have to write the updates back 
        // to all of the small cells
        if(first_side->is_hanging){
            for(int i = 0; i < P4EST_HALF; i++){
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
            for(int i = 0; i < P4EST_HALF; i++){
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

    } else {
        //y-direction

        double update_first_hv = 0;
        double update_second_hv = 0;

        riemann_solver.computeNetUpdates(first->h, second->h,
                                 first->hv, second->hv,
                                 first->b, second->b,
                                 update_first_h, update_second_h,
                                 update_first_hv, update_second_hv,
                                 current_maxWaveSpeed);

        first->update_h += update_first_h;
        second->update_h += update_second_h;
        first->update_hv += update_first_hv;
        second->update_hv += update_second_hv;


        double min_size = P4EST_QUADRANT_LEN (lowestLevel);
        min_size *= (scenario->getBoundaryPos(BND_TOP) - 
                            scenario->getBoundaryPos(BND_BOTTOM)) /
                        (double) P4EST_ROOT_LEN;
        double current_maxTimestep = 0.4 * min_size / current_maxWaveSpeed;

        if(current_maxTimestep < *maxTimestep)
            *maxTimestep = current_maxTimestep;

        // if one face is hanging, we have to write the updates back to 
        // the small cells
        if(first_side->is_hanging){
            for(int i = 0; i < P4EST_HALF; i++){
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
            for(int i = 0; i < P4EST_HALF; i++){
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
    }
}

/**
 * @brief apply the updates previously calculated, 
 *        and reset the updates for the next timestep
 */
static void updateUnknowns(p4est_iter_volume_info_t * info, void *user_data){

    double *maxTimestep = (double *) user_data;
    p4est_quadrant_t *quad = (p4est_quadrant_t *) info->quad;
    cell_state *state = (cell_state *) quad->p.user_data;
    double dx = P4EST_QUADRANT_LEN (quad->level);
    dx *= (scenario->getBoundaryPos(BND_TOP) - 
                   scenario->getBoundaryPos(BND_BOTTOM)) / 
               P4EST_ROOT_LEN;

    //apply stored updates
    state->h -= *maxTimestep / dx  * state->update_h;
    state->hu -= *maxTimestep / dx * state->update_hu;
    state->hv -= *maxTimestep / dx * state->update_hv;
}



/**
 * @brief initializes the p4est and the initial state of all 
 *        grid cells according to the given SWE_Scenario.
 * 
 * @param [in] mpicomm     MPI communicator the p4est should be running on
 * @param [in] scenario  scenario, which is used during the setup.
 *
 * @return the p4est structure
 */
p4est_t *initP4est(sc_MPI_Comm mpicomm, SWE_Scenario *scenario) {

    //create p4est
    p4est_connectivity_t *conn = p4est_connectivity_new_periodic();

    p4est_t *p4est = p4est_new(mpicomm, conn, sizeof(cell_state), 
            init_cell, NULL);

    //refine tree to minimum level
    for(int level = 0; level < 5; level++) {
        p4est_refine(p4est, 0, refine_fn, init_cell);
        p4est_balance(p4est, P4EST_CONNECT_FACE, init_cell);
        p4est_partition(p4est, 0, NULL);
    }

    //calculate one iteration to decide on further refinement,
    //refine by maximum three levels
    for(int i = 0; i < 3; i++) {
        p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);
        cell_state *ghost_data = P4EST_ALLOC(cell_state, 
                ghost->ghosts.elem_count);
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);
        double maxTimestep = 100;
        ghost_and_data *param = new ghost_and_data {ghost_data, maxTimestep};
        
        p4est_iterate (p4est, ghost,
                (void *) param,
                NULL,
                computeNumericalFluxes,
                NULL);

        p4est_refine(p4est, 0, refine_fn, init_cell);
        p4est_balance(p4est, P4EST_CONNECT_FACE, init_cell);
        p4est_partition(p4est, 0, NULL);

        p4est_ghost_destroy(ghost);
        P4EST_FREE(ghost_data);
        ghost = NULL;
        ghost_data = NULL;
    }

    return p4est;
}

/**
 * @brief main simulation loop between two phases
 *
 * @param state    the current simulation state
 * @return the actual end time reached
 */
float simulate_interval(simulation_state *state) {

    p4est_t *p4est = state->p4est;
    float t = state->t;

    // repartition the p4est
    p4est_partition(p4est, 0, NULL);

    // create the ghost quadrants
    p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);

    // create space for storing the ghost data
    cell_state *ghost_data = P4EST_ALLOC(cell_state, 
            ghost->ghosts.elem_count);

    int iterations = 0;

    do {
        // reset updates
        p4est_iterate (p4est, NULL,
                NULL,
                resetUpdates,
                NULL,
                NULL);

        // synchronize the ghost data
        p4est_ghost_exchange_data(p4est, ghost, ghost_data);

        // reset maxTimestep
        double maxTimestep = 100;

        ghost_and_data *param = new ghost_and_data {ghost_data, maxTimestep};

        // compute numerical fluxes for every edge
        // updates maxTimestep
        p4est_iterate (p4est, ghost,
                (void *) param,
                NULL,
                computeNumericalFluxes,
                NULL);
        
        //collect the timestep size from all processes
        int mpiret = sc_MPI_Allreduce (&(param->timestep), &maxTimestep, 1, 
                sc_MPI_DOUBLE, sc_MPI_MIN, state->mpicomm);
        SC_CHECK_MPI (mpiret);

        // update unknowns accordingly
        p4est_iterate (p4est, NULL,
                (void *) &maxTimestep,
                updateUnknowns,
                NULL,
                NULL);

        delete param;
        t += maxTimestep;

        iterations ++;

    } while (t <= state->phases[state->current_phase]);

    p4est_ghost_destroy(ghost);
    P4EST_FREE(ghost_data);
    ghost = NULL;
    ghost_data = NULL;

    // adapt the grid
    for(int i = 0; i < 2; i++) {
        p4est_refine_ext(p4est, 0, 10, refine_fn, NULL, replace_fn);
    }
    for(int i = 0; i < 5; i++) {
        p4est_coarsen_ext(p4est, 0, 0, coarsen_fn, NULL, replace_fn);
    }

    p4est_balance_ext(p4est, P4EST_CONNECT_FACE, NULL, replace_fn);

    return t;
}




/////////////////////////////////////////////////////////////////
////////////////////// Main functionality ///////////////////////
/////////////////////////////////////////////////////////////////



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

    // set up MPI environment for resource adaptation
    MPI_Comm world = MPI_COMM_WORLD;
    MPI_Comm mpicomm = state->mpicomm;
    MPI_Comm idle = MPI_COMM_NULL;

    MPI_Group world_group;
    MPI_Group active_group = state->active;
    MPI_Group idle_group;
    MPI_Comm_group(world, &world_group);
    MPI_Group_difference(world_group, active_group, &idle_group);

    MPI_Comm_create(world, idle_group, &idle);

    MPI_Group added_group = MPI_GROUP_NULL;

    io::P4est_vtkWriter *l_writer = 
        new io::P4est_vtkWriter(state->output_name);

    // write initial state
    if(state->mpicomm != MPI_COMM_NULL && state->current_phase == 0) {
        l_writer->writeTimeStep(state->p4est);
        state->current_phase = 1;
    }
    l_writer->set_timestep(state->current_phase);


    //set up initial MPI communicators and groups
    const int ranks[] = {0};
    int rank;
    MPI_Comm_rank(world, &rank);

    // main simulation loop over phases
    for(; state->current_phase <= state->num_phases 
            && state->t <= state->t_end; state->current_phase++) {

        if(mpicomm != MPI_COMM_NULL) {

            // simulate until next timestep
            float new_time = simulate_interval(state);

            // update simulation time with time step width.
            state->t = new_time;

            state->is_new = false;

            l_writer->writeTimeStep(state->p4est);
        }
        MPI_Barrier(world);

        // adapt resources
        if(idle_group != MPI_GROUP_EMPTY) {
            MPI_Group_incl(idle_group, 1, ranks, &added_group);
            MPI_Group tmp = MPI_GROUP_EMPTY;
            MPI_Group_union(added_group, active_group, &tmp);
            active_group = tmp;
            MPI_Group tmp2 = MPI_GROUP_EMPTY;
            MPI_Group_excl(idle_group, 1, ranks, &tmp2);
            idle_group = tmp2;
            MPI_Comm_create(world, active_group, &mpicomm);
            state->mpicomm = mpicomm;
            MPI_Comm_create(world, idle_group, &idle);
            MPI_Comm_rank(world, &rank);

            if(mpicomm != MPI_COMM_NULL) {
                int active_rank;
                MPI_Comm_rank(mpicomm, &active_rank);
                // synchronize dynamic state variables
                if(active_rank == 1) {
                    MPI_Send((void*) &state->t, 1, MPI_FLOAT, 0, 
                            0, mpicomm);

                    MPI_Send((void*) &state->current_phase, 1, 
                            MPI_INT, 0, 1, mpicomm);
                }
                if(active_rank == 0) {
                    MPI_Recv((void*) &state->t, 1, MPI_FLOAT, 1, 
                            0, mpicomm, MPI_STATUS_IGNORE);

                    MPI_Recv((void*) &state->current_phase, 1, 
                            MPI_INT, 1, 1, mpicomm, MPI_STATUS_IGNORE);
                }

                state->p4est = p4est_dynres_replace(state->p4est, mpicomm); 

                l_writer->set_timestep(state->current_phase);
                
            }
        }
    }

    //Finalize
    delete l_writer;
    return 0;
}



/**
 * @brief main program for the simulation on a p4est
 *
 * Initializes new processes and starts the simulation
 */
int main(int argc, char **argv) {

    //initialize MPI
    sc_MPI_Comm mpicomm;
    int mpiret = sc_MPI_Init(&argc, &argv);
    SC_CHECK_MPI(mpiret);
    mpicomm=sc_MPI_COMM_WORLD;

    simulation_state state;

    // Parse command line parameters
    tools::Args args;

    //parameters for regular startup
    args.addOption("scenario", 's', "Loads input conditions, (artificial | dambreak), default: dambreak", args.Required, false);
    args.addOption("bathymetry-path", 'b', "Loads bathymetry input from a file, overwrites -s option", args.Required, false);
    args.addOption("displacement-path", 'd', "Loads displacement data from a file, to be used in combination with -b", args.Required, false);
    args.addOption("output-basepath", 'o', "Output base file name, mandatory argument, no default.");
    args.addOption("simulated-time", 't', "Time in seconds to be simulated, default: 200", args.Required, false);
    args.addOption("num-phases", 'c', "Number of phases for visualization (at each phase in time, an output file is written), default: 200", args.Required, false);

    //additional parameters for processes spawned during simulation
    args.addOption("current-time", 'n', "Current simulation time, not to be set manyally", args.Required, false);
    args.addOption("current-phase", 'i', "Number of current phase, not to be set manually", args.Required, false);

    tools::Args::Result ret = args.parse(argc, argv);

    switch (ret) {
        case tools::Args::Error:
            exit(1);
        case tools::Args::Help:
            exit(0);
        default:
            break; 
    }

    //has this process been spawned during simulation?
    bool joining = args.isSet("current-time");

    // read command line parameters
    strcpy(state.output_name, 
            args.getArgument<std::string>("output-basepath").c_str());
    state.num_phases = args.getArgument<int>("num-phases", 200);
    state.t_end = args.getArgument<float>("simulated-time", 200);
    state.phases = new float[state.num_phases + 1];
    for (int cp = 0; cp <= state.num_phases; cp++) {
        state.phases[cp] = cp * (state.t_end / state.num_phases);
    }

    state.is_new = true;

    //initialize scenario
    std::string scenario_option = "";
    std::string bathymetry_path = "";
    std::string displacement_path = "";

    if(!args.isSet("bathymetry-path") && !args.isSet("displacement-path")) {
        scenario_option = args.getArgument<std::string>("scenario", "dambreak");
    } else if(args.isSet("bathymetry-path") && args.isSet("displacement-path")) {
        bathymetry_path = args.getArgument<std::string>("bathymetry-path");
        displacement_path = args.getArgument<std::string>("displacement-path");
    } else {
        cout << "The options -b and -d must be used together. Aborting.\n";
        exit(1);
    }

    strcpy(state.scenario, scenario_option.c_str());
    strcpy(state.bathymetry, bathymetry_path.c_str());
    strcpy(state.displacement, displacement_path.c_str());

    if(strcmp(state.scenario, "artificial") == 0)
        scenario = new SWE_ArtificialTsunamiScenario();
    else if(strcmp(state.scenario, "dambreak") == 0)
        scenario = new SWE_RadialDamBreakScenario();
    else
#ifdef NETCDF
        scenario = new SWE_TsunamiScenario(bathymetry_path.c_str(), 
                displacement_path.c_str(), state.t_end,
                OUTFLOW, OUTFLOW, OUTFLOW, OUTFLOW);
#else
        std::cout << "ERROR: In order to read the scenario from " <<
                     "netCDF files, compile with netCDF=true\n";
#endif

    // create and init the p4est simulation object on a single process

    const int ranks[] = {0};
    MPI_Group start_group;
    MPI_Group world_group;
    MPI_Comm start_comm;
    MPI_Comm_group(mpicomm, &world_group);
    MPI_Group_incl(world_group, 1, ranks, &start_group);
    MPI_Comm_create(mpicomm, start_group, &start_comm);
    if(start_comm != MPI_COMM_NULL)
        state.p4est = initP4est(start_comm, scenario);
    else
        state.p4est = NULL;
    
    state.mpicomm = start_comm;
    state.active = start_group;

    if(!joining) {
        //initialize other state variables
        state.t = 0;
        state.current_phase = 0;
    } else {
        // read command line parameters
        state.t = args.getArgument<float>("current-time");
        state.current_phase = args.getArgument<int>("current-phase");
    }

    simulate(&state);

    if(state.p4est != NULL) {
        p4est_connectivity_t *connectivity = state.p4est->connectivity;
        p4est_destroy(state.p4est);
        p4est_connectivity_destroy(connectivity);
    }
    delete[] state.phases;

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);

    return 0;
}

