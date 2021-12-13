/**
 * @file benchmarks/p4est_benchmark_constProcs_adapt.cpp
 * 
 * Benchmark for changing number of quadrants on a constant number
 * of processes. The grid is refined according to the wave steepness 
 * at every phase.
 * Initial number fo quadrants is 4^9, min 4^6, max 4^11.
 * The number of processes is 56 on two nodes.
 * Simulated with 10 phases.
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
    p4est_t *p4est;                 // the p4est simulation object
    sc_MPI_Comm mpicomm;            // the MPI communicator of all working processes
    bool is_new;                    // flag showing whether this process has been 
                                    // running for more than one iteration
    
    int iterations;
    clock_t time;
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


static solver::AugRie<double> riemann_solver;

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
    bool refine = false;

    // weighting factor depending on the quadrant's level, 
    // for iterative refining
    double factor = 1;
    factor = pow(2, ((int) quadrant->level - 10));

    return (quadrant->level < 7 || 
            (quadrant->level < 11 && state->update_h / (double) factor > 0.2));
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
    factor = pow(2, ((int) quadrant[0]->level - 10));

    return (quadrant[0]->level > 11 ||
            ((quadrant[0]->level > 7) && update_h_max / (double) factor < 0.1));
}

/*
 * @brief Function for initializing cells after adapting the grid
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

    if(p4est->user_pointer == NULL)
        return;

    SWE_Scenario *scenario = (SWE_Scenario *) p4est->user_pointer;

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

    SWE_Scenario *scenario = (SWE_Scenario *) info->p4est->user_pointer;

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
 * @brief initializes the unknowns and bathymetry 
 *        in all grid cells according to the given SWE_Scenario.
 * 
 * @param i_scenario scenario, which is used during the setup.
 */
p4est_t *initP4est(sc_MPI_Comm mpicomm, SWE_Scenario *scenario) {

    //create p4est
    p4est_connectivity_t *conn = p4est_connectivity_new_periodic();

    p4est_t *p4est = p4est_new(mpicomm, conn, sizeof(cell_state), 
            init_cell, scenario);

    //refine tree to minimum level
    for(int level = 0; level < 9; level++) {
        p4est_refine(p4est, 0, refine_fn, init_cell);
        p4est_balance(p4est, P4EST_CONNECT_FACE, init_cell);
        p4est_partition(p4est, 0, NULL);
    }

    return p4est;
}

/**
 * @brief apply the updates previously calculated, 
 *        and reset the updates for the next timestep
 */
static void updateUnknowns(p4est_iter_volume_info_t * info, void *user_data){

    SWE_Scenario *scenario = (SWE_Scenario *) info->p4est->user_pointer;

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
 * @brief main simulation loop between two phases
 *
 * @param	tStart	time where the simulation is started
 * @param	tEnd	time of the next phase 
 * @return	actual	end time reached
 */
float simulate_interval(simulation_state *state) {

    p4est_t *p4est = state->p4est;
    float t = state->t;
    state->time = 0;

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
    clock_t start = clock();

    for(int i = 0; i < 4; i++) {
        p4est_refine_ext(p4est, 0, 11, refine_fn, NULL, replace_fn);
    }
    for(int i = 0; i < 4; i++) {
        p4est_coarsen_ext(p4est, 0, 0, coarsen_fn, NULL, replace_fn);
    }

    p4est_balance_ext(p4est, P4EST_CONNECT_FACE, NULL, replace_fn);
    clock_t end = clock();
    state->time = end-start;
    state->iterations = iterations;

    return t;
}



/////////////////////////////////////////////////////////////////
////////////////////// Main functionality ///////////////////////
/////////////////////////////////////////////////////////////////



/**
 * @brief main simulation loop
 */
int simulate(simulation_state *state){

    state->current_phase = 1;

    unsigned int l_iterations = 0;


    ofstream output_file;
    if(state->p4est->mpirank == 0) {
        output_file.open("constProcs_adapt_log.csv");
        output_file << 
            "Checkpoint,Quadrants,Processes,Iterations,TotalTime,p4estTime,partitionTime,dynresTime\n";
    }


    // main simulation loop over phases
    for(; state->current_phase <= state->num_phases 
            && state->t <= state->t_end; state->current_phase++) {

        clock_t phase_start = clock();

        // simulate until next timestep
        float new_time = simulate_interval(state);

        clock_t partition_start = clock();
        p4est_partition(state->p4est, 0, NULL);
        clock_t partition_end = clock();

        // update simulation time with time step width.
        state->t = new_time;

        state->is_new = false;

        //log time mesurements
        clock_t phase_end = clock();

        if(state->p4est->mpirank == 0) {
            output_file << state->current_phase << ",";
            output_file << state->p4est->global_num_quadrants << ",";
            output_file << state->p4est->mpisize << ",";
            output_file << state->iterations << ",";
            output_file << ((double) (phase_end - phase_start)) / CLOCKS_PER_SEC << ",";
            output_file << ((double) state->time) / CLOCKS_PER_SEC << ",";
            output_file << ((double) (partition_end - partition_start)) / CLOCKS_PER_SEC << ",";
            output_file << "0\n";
        }
        MPI_Barrier(state->mpicomm);

    }

    //Finalize
    if(state->p4est->mpirank == 0) {
        output_file.close();
    }
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
    args.addOption("simulated-time", 't', "Time in seconds to be simulated, default: 200", args.Required, false);
    args.addOption("num-phases", 'c', "Number of phases for visualization (at each phase in time, an output file is written), default: 200", args.Required, false);

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
    state.num_phases = args.getArgument<int>("num-phases", 200);
    state.t_end = args.getArgument<float>("simulated-time", 200);
    state.phases = new float[state.num_phases + 1];
    for (int cp = 0; cp <= state.num_phases; cp++) {
        state.phases[cp] = cp * (state.t_end / state.num_phases);
    }

    state.is_new = true;

    SWE_Scenario *l_scenario = new SWE_RadialDamBreakScenario();

    state.mpicomm = mpicomm;

    // create and init the p4est simulation object
    state.p4est = initP4est(mpicomm, l_scenario);

    //initialize other state variables
    state.t = 0;
    state.current_phase = 0;

    simulate(&state);

    delete l_scenario;
    if(state.p4est != NULL) {
        p4est_connectivity_destroy(state.p4est->connectivity);
        p4est_destroy(state.p4est);
    }
    delete[] state.phases;

    mpiret = sc_MPI_Finalize ();
    SC_CHECK_MPI (mpiret);

    return 0;
}

