/**
 * @file benchmarks/p4est_benchmark_decMpidynres_fixed.cpp
 * 
 * Benchmark for p4est on libmpidynres using a constant 
 * number of quadrants.
 * Quadrants: 4^10
 * Processes: decreasing from 55 to 45, removing one process per timestep
 * Start with 56 processes, as one is required by libmpidynres!
 * Use 10 phases.
 * 
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../swe_solvers/src/solver/AugRie.hpp"

#include "scenarios/SWE_Scenario.hh"
#include "scenarios/SWE_simple_scenarios.hh"

#include "tools/args.hh"

#include <../p4est/src/p4est_geometry.h>
#include <../p4est/src/p4est_extended.h>
#include <../p4est/src/p4est_communication.h>

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
    float *phases;             // an array containing all phases in time
    int num_phases;            // the number of phases in time
    int current_phase;         // the index of the current phase
    SWE_Scenario *scenario;          // the simulated scenario
    p4est_t *p4est;                 // the p4est simulation object
    MPI_Session session;            // the MPI Session object
    char pset[MPI_MAX_PSET_NAME_LEN];   // the curent process set
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


    return 1;
}

/**
 * @brief initialize a cell using the scenario data
 */
static void init_cell(p4est_t *p4est, p4est_topidx_t tree, 
        p4est_quadrant_t *quadrant){

    SWE_Scenario *scenario = (SWE_Scenario *) p4est->user_pointer;

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
    for(int level = 0; level < 10; level++) {
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
    int iterations = 0;

    // create the ghost quadrants
    p4est_ghost_t *ghost = p4est_ghost_new(p4est, P4EST_CONNECT_FULL);

    // create space for storing the ghost data
    cell_state *ghost_data = P4EST_ALLOC(cell_state, 
            ghost->ghosts.elem_count);

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

        iterations++;

    } while (t <= state->phases[state->current_phase]);

    p4est_ghost_destroy(ghost);
    P4EST_FREE(ghost_data);
    ghost = NULL;
    ghost_data = NULL;

    // adapt the grid
    state->time = 0;
    state->iterations = iterations;

    return t;
}



////////////////////////////////////////////////////////////////
//////////// Functions for ressource management ////////////////
////////////////////////////////////////////////////////////////


/*
 * @brief spawn and new processes with correct input parameters
 *
 * @param state the current simulation state
 * @param n     the number of processes to be spawned
 *
 * @return an intercommunicator containing all processes
 */
MPI_Info create_info_from_state(simulation_state *state) {
   
    MPI_Info info;
    MPI_Info_create(&info);

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
    
    return info;
}


bool check_resource_change(simulation_state *state) {

    bool must_terminate = false;

    int mpirank = -1;
    if(state->p4est != NULL)
        mpirank = state->p4est->mpirank;

    MPI_Info rc_info;
    MPIDYNRES_RC_type rc_type = MPIDYNRES_RC_ADD;
    char delta_pset[MPI_MAX_PSET_NAME_LEN];
    int rc_tag;

    if(mpirank == 0) {
        // request a resource change: add one process
        DEBUG_MPIDYNRES_RC_user_request(state->session, -1, delta_pset, &rc_type, &rc_tag, &rc_info);
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
                MPI_Session_get_psets(state->session, MPI_INFO_NULL, 
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
        MPI_Group_from_session_pset(state->session, state->pset, &group);
        if(!must_terminate) {
            MPI_Comm_create_from_group(group, NULL, MPI_INFO_NULL,
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


simulation_state initialize(int argc, char* argv[]) {

    simulation_state state;

    int contains_key;
    int value_len;
    MPI_Info session_info;
    MPI_Info psets;

    //init MPI Session
    MPI_Session_init(MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &state.session);
    MPI_Session_get_psets(state.session, MPI_INFO_NULL, &psets);


    //check whether this process was spawned at startup time
    MPI_Info_get_valuelen(psets, "mpi://WORLD", &value_len, &contains_key);
    if(contains_key) {

        //initialize using cli parameters
        tools::Args args;

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
        state.num_phases = args.getArgument<int>("num-phases", 200);
        state.t_end = args.getArgument<float>("simulated-time", 200);

        state.t = 0;
        state.current_phase = 1;
        strcpy(state.pset, "mpi://WORLD");

    } else {

        // initialize using the session info object
        MPI_Info session_info;
        MPI_Session_get_info(state.session, &session_info);

        int exists;

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

    }

    state.phases = new float[state.num_phases + 1];
    for (int cp = 0; cp <= state.num_phases; cp++) {
        state.phases[cp] = cp * (state.t_end / state.num_phases);
    }

    state.is_new = true;

    //initialize scenario
    SWE_Scenario *scenario = new SWE_RadialDamBreakScenario();
    state.scenario = scenario;


    //create mpicomm
    int mpiret;
    MPI_Group group;
    mpiret = MPI_Group_from_session_pset(state.session, state.pset, &group);
    if(mpiret)
        MPIDYNRES_exit();

    mpiret = MPI_Comm_create_from_group(group, NULL, MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &state.mpicomm);
    if(mpiret)
        MPIDYNRES_exit();

    // create and init the p4est simulation object
    state.p4est = NULL;
    if(contains_key) {
        state.p4est = initP4est(state.mpicomm, scenario);
        state.is_new = false;
    }

    MPI_Group_free(&group);

    return state;
}


void finalize(simulation_state *state) {

    if(state->p4est != NULL) {
        p4est_connectivity_t *connectivity = state->p4est->connectivity;
        p4est_destroy(state->p4est);
        p4est_connectivity_destroy(connectivity);
    }
    delete[] state->phases;
}


/**
 * @brief main work loop
 *
 * calls simulate_interval
 */
int simulate(simulation_state *state){

    unsigned int l_iterations = 0;
    clock_t simulation_time = 0;
    clock_t dynres_time = 0;

    ofstream output_file;
    int rank;
    MPI_Comm_rank(state->mpicomm, &rank);

    if(rank == 0) {
        output_file.open("decMpidynres_fixed_log.csv");
        output_file << 
            "Checkpoint,Quadrants,Processes,Iterations,TotalTime,refinementTime,partitionTime,dynresTime\n";
    }

    // main simulation loop over phases
    for(; state->current_phase <= state->num_phases && 
                state->t <= state->t_end; 
            state->current_phase++) {

        dynres_time = clock();

        // handle dynamic resources
        bool changed = state->is_new;
        if(!state->is_new)
            changed = check_resource_change(state);
        if(changed) {
            state->p4est = p4est_dynres_replace_ext(state->p4est, state->mpicomm, 
                    -1, -1, -1, -1, (void *) state->scenario, NULL);
        }

        if(state->p4est == NULL)
            return 0;

        clock_t dynres_end = clock();
        dynres_time = dynres_end - dynres_time;
 
        clock_t partition_start = clock();
        p4est_partition(state->p4est, 0, NULL);
        clock_t partition_end = clock();
       
        simulation_time = clock();

        // simulate until next timestep
        float new_time = simulate_interval(state);

        // update simulation time with time step width.
        state->t = new_time;
        l_iterations++;

        state->is_new = false;

        clock_t simulation_end = clock();
        simulation_time = simulation_end - simulation_time;

        if(rank == 0) {
            output_file << state->current_phase << ",";
            output_file << state->p4est->global_num_quadrants << ",";
            output_file << state->p4est->mpisize << ",";
            output_file << state->iterations << ",";
            output_file << ((double) (simulation_time)) / CLOCKS_PER_SEC << ",";
            output_file << ((double) state->time) / CLOCKS_PER_SEC << ",";
            output_file << ((double) (partition_end - partition_start)) / CLOCKS_PER_SEC << ",";
            output_file << ((double) dynres_time) / CLOCKS_PER_SEC << "\n";
        }
    }

    return 0;
}


int MPIDYNRES_main(int argc, char **argv) {

    simulation_state state = initialize(argc, argv);


    simulate(&state);


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
    int mpisize;
    MPI_Comm_size(mpicomm, &mpisize);

    //initialize mpidynres
    MPI_Info manager_config;
    MPI_Info_create(&manager_config);
    char buf[0x20];
    snprintf(buf, 0x20, "%d", mpisize-1);
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

