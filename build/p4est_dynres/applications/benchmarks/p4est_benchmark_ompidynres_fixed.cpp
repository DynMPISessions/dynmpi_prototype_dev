/**
 * @file benchmarks/p4est_benchmark_mpidynres_fixed.cpp
 * 
 * Benchmark for p4est on libmpidynres using a constant 
 * number of quadrants.
 * Quadrants: 4^10
 * Processes: increasing from 45 to 55, adding one process per timestep
 * Start with 56 processes, as one is required by libmpidynres!
 * Use 10 phases.
 * 
 */

#include <cassert>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cstring>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include "../swe_solvers/src/solver/AugRie.hpp"

#include "scenarios/SWE_Scenario.hh"
#include "scenarios/SWE_simple_scenarios.hh"

#include "tools/args.hh"
extern "C" {
#include "tools/bench_timings.h"
}

#include <../p4est/src/p4est_geometry.h>
#include <../p4est/src/p4est_extended.h>
#include <../p4est/src/p4est_communication.h>

extern "C" {
#include "comm.h"
}
extern "C" {
#include <mpi.h>
}


#define BENCH_INCRIMENTAL 0
#define BENCH_SEQUENTIAL 1
#define BENCH_MIXED 2

static int rc_id = 0;
static int num_delta = 1;
static int cur_num_delta = 1;
static int proc_limit = -1;
static int proc_limit_down = -1;
static int rc_frequency = 0;
static int recvd_delta = -1;
static int blocking = 0;
static bool spawned = false;
static bool check_rc = true;

static char mode[3];
static char mode_str[64];
char cwd[PATH_MAX];
static int mode_num = BENCH_INCRIMENTAL;
static int mode_type = MPI_RC_ADD;
static int cur_type = MPI_RC_ADD;

int set_mode(const char *c_mode){
    if(strlen(c_mode) != 2){
        return -1;
    }
    strcpy(mode, c_mode);
    if(c_mode[0] == 'i'){
        mode_num = BENCH_INCRIMENTAL;
        cur_num_delta = num_delta;
    }else if(c_mode[0] == 's'){
        mode_num = BENCH_SEQUENTIAL;
        cur_num_delta = num_delta;
    }else if(c_mode[0] == 'm'){
        mode_num = BENCH_MIXED;
        cur_num_delta = proc_limit / num_delta * num_delta;
    }else{
        return -1;
    }
    strcpy(mode_str, mode);


    if(c_mode[1] == '+'){
        mode_type = MPI_RC_ADD;
        strcpy(mode_str + 1, ":plus");
    }else if(c_mode[1] == '_'){
        mode_type = MPI_RC_SUB;
        strcpy(mode_str + 1, ":minus");
        if(mode_num == BENCH_MIXED){
            return -1;
        }
    }else{
        return -1;
    }
    cur_type = mode_type;

    return 0;

}

void eval_mode(int current_size,int  mode_n){
    switch(mode_n){

        case BENCH_MIXED:
            cur_type = cur_type == MPI_RC_ADD ? MPI_RC_SUB : MPI_RC_ADD;
            cur_num_delta = cur_type != mode_type ? cur_num_delta - num_delta : cur_num_delta;           
            check_rc = cur_num_delta <= 0 ? false : cur_type == MPI_RC_ADD ? current_size + cur_num_delta <= proc_limit : (current_size - cur_num_delta >= 1);
            break;
        case BENCH_SEQUENTIAL:
            cur_type = cur_type == MPI_RC_ADD ? MPI_RC_SUB : MPI_RC_ADD;
            cur_num_delta = cur_type == mode_type ? cur_num_delta + num_delta: cur_num_delta;
            check_rc = cur_num_delta <= 0 ? false : cur_type == MPI_RC_ADD ? current_size + cur_num_delta <= proc_limit : (current_size - cur_num_delta >= proc_limit_down);
            break;
        case BENCH_INCRIMENTAL:
            check_rc = mode_type == MPI_RC_ADD ? current_size + cur_num_delta <= proc_limit : (current_size - cur_num_delta >= proc_limit && current_size - cur_num_delta >= 1);
            break;
        default:
            if(mode_num > BENCH_INCRIMENTAL){
                if(mode_type == MPI_RC_ADD && spawned){
                    cur_type = cur_type == MPI_RC_ADD ? MPI_RC_SUB : MPI_RC_ADD;
                }
                if(mode_num == BENCH_MIXED){
                    if(recvd_delta < 0)
                        cur_num_delta = (proc_limit - current_size) / num_delta * num_delta;
                    else 
                        cur_num_delta = recvd_delta;
                }
                check_rc = cur_num_delta <= 0 ? false : cur_type == MPI_RC_ADD ? current_size + cur_num_delta <= proc_limit : (current_size - cur_num_delta >= proc_limit_down);
            }else{
                check_rc = mode_type == MPI_RC_ADD ? current_size + cur_num_delta <= proc_limit : (current_size - cur_num_delta >= proc_limit && current_size - cur_num_delta >= 1);

            }
            break;            
    }
    //printf("cur num delta init %d, check_rc %d, cur type %d, current size %d\n", cur_num_delta, check_rc, cur_type, current_size);
}

node_t *base_timing_list;
timing_frame_base_t *cur_base_timing_frame;

node_t *p4est_timing_list;
timing_frame_p4est_t *cur_p4est_timing_frame;


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


int check_resource_change(simulation_state *state) {

    bool must_terminate = false;
    int terminate = 0;
    MPI_Comm new_comm;
    struct timespec ts;
    ts.tv_sec = 0;
    ts.tv_nsec = 10000;

    int mpirank = -1;
    if(state->p4est != NULL)
        mpirank = state->p4est->mpirank;

    MPI_Comm_rank(state->mpicomm, &mpirank);

    MPI_Info rc_info=MPI_INFO_NULL;
    int rc_type = MPI_RC_SUB;
    int rc_type_get= MPI_RC_NULL, status, rc = MPI_ERR_OTHER;
    char delta_pset[MPI_MAX_PSET_NAME_LEN];
    int incl_flag;

    //request and query the resource change
    if(mpirank == 0) {
    
        make_timestamp_root(&cur_p4est_timing_frame->get_rc_start_1);
        rc = MPI_Session_get_res_change(state->session, &rc_type_get, delta_pset, &incl_flag, &status, NULL);
        if(MPI_RC_NULL == rc_type_get){
            // if no rc and not requested: request a resource change: add one process
            make_timestamp_root(&cur_p4est_timing_frame->request_start);
            MPI_Session_request_res_change(state->session, cur_num_delta, delta_pset, cur_type, &rc_info);
            // As the request returns when notification is send, do an ugly wait until the rc is defined 
            make_timestamp_root(&cur_p4est_timing_frame->get_rc_start_2);
            while(MPI_RC_NULL == rc_type_get){
                rc = MPI_Session_get_res_change(state->session, &rc_type_get, delta_pset, &incl_flag, &status, NULL);
                if(MPI_RC_NULL == rc_type_get){
                    nanosleep(&ts, NULL);
                }
            }
        }

        if(MPI_RC_NULL != rc_type_get){
            timings_cur_rc_status = status;
            if(MPI_RC_ANNOUNCED == status){
                set_res_change_id(&cur_p4est_timing_frame->res_change_id, delta_pset);
            }
        }

        if(rc_type_get == MPI_RC_NULL || (status != MPI_RC_ANNOUNCED && status != MPI_RC_CONFIRMATION_PENDING)){
            rc = MPI_ERR_OTHER;
        }
        // ignore rc_info
        if(rc_info != MPI_INFO_NULL)
            MPI_Info_free(&rc_info);
    }

    MPI_Bcast(&rc, 1, MPI_INT, 0, state->mpicomm);
    MPI_Bcast(&rc_type_get, 1, MPI_INT, 0, state->mpicomm);
    if(rc == MPI_SUCCESS && rc_type_get == MPI_RC_SUB){
        MPI_Comm_disconnect(&state->p4est->mpicomm);
        MPI_Group group;
        MPI_Group_from_session_pset(state->session, state->pset, &group);
        MPI_Comm_create_from_group(group, "tag", MPI_INFO_NULL,
            MPI_ERRORS_ARE_FATAL, &state->p4est->mpicomm);
        MPI_Group_free(&group);
    }
    
    // handle changes
    if(mpirank == 0 && rc_type_get != MPI_RC_NULL && status == MPI_RC_ANNOUNCED) {
        make_timestamp_root(&cur_p4est_timing_frame->psetop_start);
        switch(rc_type_get) {
            case MPI_RC_ADD:
                //create new pset
                if(mpirank == 0) 
                    MPI_Session_pset_create_op(state->session, MPI_PSETOP_UNION,
                            state->pset, delta_pset, 
                            state->pset, NULL);
                break;
            case MPI_RC_SUB:
                //create new pset
                if(mpirank == 0) 
                    MPI_Session_pset_create_op(state->session, MPI_PSETOP_DIFFERENCE,
                            state->pset, delta_pset, 
                            state->pset, NULL);
                break;
            default:
                break;
        }
        //printf("out of pset_create_op\n");
    }

    make_timestamp_root(&cur_p4est_timing_frame->handle_rc_start);
    // bcast changes to all processes
    //MPI_Bcast(&rc, 1, MPI_INT, 0, state->mpicomm);
    if(MPI_SUCCESS != rc){
        return rc;
    }
    
    //MPI_Bcast(&rc_type_get, 1, MPI_INT, 0, state->mpicomm);
    
    // if we remove processes we need to replace the p4est communicator now
    if(rc_type_get == MPI_RC_SUB){
        // this is ugly but for now use this to check if we are included in the resource change
        if(mpirank != 0){
            rc_type_get = MPI_RC_NULL;
            while(MPI_RC_NULL == rc_type_get){
                rc = MPI_Session_get_res_change(state->session, &rc_type_get, delta_pset, &incl_flag, &status, NULL);
                if(rc_type_get == MPI_RC_NULL){
                    nanosleep(&ts, NULL);
                }
            }
        }
        //printf("Rank %d: bcast\n", mpirank);
        MPI_Bcast(state->pset,MPI_MAX_PSET_NAME_LEN, MPI_CHAR, 0, state->mpicomm);
        //MPI_Barrier(state->mpicomm);


        /* replace with MPI_COMM_NULL */
        if(incl_flag){
            state->p4est = p4est_dynres_replace_ext(state->p4est, MPI_COMM_NULL,
                        -1, -1, -1, -1, (void *) state->scenario, NULL);
            //state->mpicomm = MPI_COMM_NULL;
            //MPI_Comm_disconnect(&state->mpicomm);
        /* replace with new communicator */    
        }else{
            MPI_Group group;
            //printf("Rank %d: group_from_session_pset\n", mpirank);
            MPI_Group_from_session_pset(state->session, state->pset, &group);
            MPI_Comm_create_from_group(group, "tag", MPI_INFO_NULL,
                MPI_ERRORS_ARE_FATAL, &new_comm);
            MPI_Group_free(&group);
            state->p4est = p4est_dynres_replace_ext(state->p4est, new_comm, 
                        -1, -1, -1, -1, (void *) state->scenario, NULL);
            
        }

        make_timestamp_root(&cur_p4est_timing_frame->handle_rc_end);
    }

    
    
    //________________________________ACCEPT____________________//
    MPI_Info info = MPI_INFO_NULL;
    MPI_Info_create(&info);
    MPI_Info_set(info, "mpi_no_disconnect", "1");
    if(blocking){
        MPI_Info_set(info, "mpi_blocking", "1");
    }
    if(rc_type_get == MPI_RC_SUB){
        MPI_Info_set(info, "mpi_rc_type", "2");
    }
    make_timestamp_root(&cur_p4est_timing_frame->accept_start);
    //if(0 == mpirank)
    //    printf("calling accept with %s and %s\n", delta_pset, state->pset);
    if(rc_type_get == MPI_RC_SUB){
        MPI_Comm_disconnect(&state->mpicomm);
        if(!incl_flag){
            rc = MPI_Session_accept_res_change(&state->session, &info, delta_pset, state->pset, 0, &new_comm, &terminate);
            MPI_Comm_disconnect(&new_comm);
        }else{
            rc = MPI_Session_accept_res_change(&state->session, &info, delta_pset, state->pset, 0, NULL, &terminate);
            return MPI_SUCCESS;
        }
    }else{
        rc = MPI_Session_accept_res_change(&state->session, &info, delta_pset, state->pset, 0, &state->mpicomm, &terminate);
    }
    make_timestamp_root(&cur_p4est_timing_frame->accept_end);
    //________________________________ACCEPT_________________//


    if(MPI_SUCCESS != rc){
        //printf("rc not success: will terminate!\n");
        return rc;
    }

    //if we added processes we need to create a new communicator now
    if(rc_type_get == MPI_RC_ADD){
        make_timestamp_root(&cur_p4est_timing_frame->handle_rc_start);
        MPI_Group group;
        int retv = MPI_Group_from_session_pset(state->session, state->pset, &group);
        MPI_Comm_create_from_group(group, "tag", MPI_INFO_NULL,
            MPI_ERRORS_ARE_FATAL, &state->mpicomm);
        MPI_Info state_info = create_info_from_state(state);
        MPIDYNRES_Bcast_Send_MPI_Info(state_info, 0, state->mpicomm);
        if(mode_num == BENCH_MIXED || (mode_num == BENCH_SEQUENTIAL && mode_type == MPI_RC_SUB)){
            MPI_Bcast(&cur_num_delta, 1, MPI_INT, 0, state->mpicomm);
        }
        MPI_Group_free(&group);
        
    }else{

        MPI_Group group;
        MPI_Group_from_session_pset(state->session, state->pset, &group);
        MPI_Comm_create_from_group(group, "tag", MPI_INFO_NULL,
            MPI_ERRORS_ARE_FATAL, &state->mpicomm);
        MPI_Group_free(&group);
        int i = 0;
    }
    
    return MPI_SUCCESS;
  

}







/////////////////////////////////////////////////////////////////
////////////////////// Main functionality ///////////////////////
/////////////////////////////////////////////////////////////////


simulation_state initialize(int argc, char* argv[]) {

    simulation_state state;

    int contains_key;
    int value_len;
    int rc_type = MPI_RC_NULL;
    int rc_status = MPI_RC_INVALID;
    int incl;
    char delta_pset[MPI_MAX_PSET_NAME_LEN];
    MPI_Info session_info;
    MPI_Info psets;
    MPI_Group group;

    //initialize using cli parameters
    tools::Args args;

    args.addOption("simulated-time", 't', 
            "Time in seconds to be simulated, default: 200", 
            args.Required, false);
    args.addOption("num-phases", 'c', 
            "Number of phases for visualization (at each phase\
                in time, an output file is written), default: 200", 
            args.Required, false);
    args.addOption("rc-mode", 'm', 
            "The mode for resource changes: i (incremental), s (sequential), m (mixed)\
                appended by +/- for adding/substracting resources, default: i+", 
            args.Required, false);
    args.addOption("num-delta", 'n', 
            "Number of processes per resource change default: 1", 
            args.Required, false);
    args.addOption("proc-limit", 'l', 
            "Upper/Lower limit of processes. default: -1", 
            args.Required, false);
    args.addOption("rc-frequency", 'f', 
            "Number of constant phases after resource change. default: 0", 
            args.Required, false);
    args.addOption("blocking", 'b', 
            "Wether MPI_Session_accept_res_change should block. default: 0", 
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

    num_delta = args.getArgument<int>("num-delta", 1);
    proc_limit = args.getArgument<int>("proc-limit", -1);
    rc_frequency = args.getArgument<int>("rc-frequency", 0);
    blocking = args.getArgument<int>("blocking", 0);
    std::string str_mode = args.getArgument<std::string>("rc-mode", "i+");
    const char *c_mode = str_mode.c_str();
    if(0 != set_mode(c_mode)){
        printf("invalid arguments\n");
        _exit(-1);
    }

    //printf("Read paramters: mode=%s(%d:%d), num_delta=%d, proc_limit=%d, rc_freq=%d\n", mode, mode_num, mode_type, num_delta, proc_limit, rc_frequency);




    //init MPI Session
    make_timestamp_base(&cur_base_timing_frame->init_start);
    MPI_Session_init(MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &state.session);
    make_timestamp_base(&cur_base_timing_frame->init_end);

    //get res change
    MPI_Session_get_res_change(state.session, &rc_type, delta_pset, &incl, &rc_status, NULL);
    
    int initial_size;
    if(MPI_RC_NULL == rc_type) {

        state.t = 0;
        state.current_phase = 1;
        //change to other start pset
        strcpy(state.pset, "test1");
        int retval;
        retval = MPI_Group_from_session_pset(state.session, state.pset, &group);
        MPI_Group_size(group, &initial_size);
        retval = MPI_Comm_create_from_group(group, "tag", MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &state.mpicomm);

        

    } else {

        spawned = true;

        set_res_change_id(&cur_base_timing_frame->res_change_id, delta_pset);

        //confirm the res change
        
        int retval;
        char *result_pset = (char*)malloc(MPI_MAX_PSET_NAME_LEN);

        make_timestamp_base(&cur_base_timing_frame->confirm_start);
        retval = MPI_Session_confirm_res_change(&state.session, &session_info, delta_pset, &result_pset);
        make_timestamp_base(&cur_base_timing_frame->confirm_end);

        if(retval)
            die("confirm res change returned %d\n", retval);
        
        //create communicator
        retval = MPI_Group_from_session_pset(state.session, result_pset, &group);
        MPI_Group_size(group, &initial_size);
        retval = MPI_Comm_create_from_group(group, "tag", MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &state.mpicomm);
        free(result_pset);

        //receive state info
        MPI_Info session_info;
        retval = MPIDYNRES_Bcast_Recv_MPI_Info(&session_info, 0, state.mpicomm);
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

        MPI_Info_free(&session_info);

        if(mode_num == BENCH_MIXED || (mode_num == BENCH_SEQUENTIAL && mode_type == MPI_RC_SUB)){
            MPI_Bcast(&cur_num_delta, 1, MPI_INT, 0, state.mpicomm);
            cur_num_delta += num_delta;
        }
    }

    if(0 == strcmp(c_mode, "s_")){
        proc_limit_down = proc_limit;
        proc_limit = initial_size;
        //printf("setting limits for s_ mode: up: %d, down: %d\n", proc_limit, proc_limit_down);
    }

    state.phases = new float[state.num_phases + 1];
    for (int cp = 0; cp <= state.num_phases; cp++) {
        state.phases[cp] = cp * (state.t_end / state.num_phases);
    }

    state.is_new = true;

    //initialize scenario
    SWE_Scenario *scenario = new SWE_RadialDamBreakScenario();
    state.scenario = scenario;

    // create and init the p4est simulation object
    state.p4est = NULL;
    if(MPI_RC_NULL == rc_type) {
        MPI_Group dup_group;
        MPI_Comm dup_comm = MPI_COMM_NULL;
        MPI_Group_from_session_pset(state.session, state.pset, &dup_group);
        MPI_Comm_create_from_group(group, "tag2", MPI_INFO_NULL, MPI_ERRORS_ARE_FATAL, &dup_comm);
        MPI_Group_free(&dup_group);
        state.p4est = initP4est(dup_comm, scenario);
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
   
    int rank;
    MPI_Comm_rank(state->mpicomm, &rank);
    timings_my_rank = rank;
    printf("Rank %d starting simulation\n", rank);
    init_add_timing(p4est_timing_list, (void**) &cur_p4est_timing_frame, sizeof(timing_frame_p4est_t));

    make_timestamp_base(&cur_base_timing_frame->sim_start);
    

    // main simulation loop over phases
    for(; state->current_phase <= state->num_phases && state->t <= state->t_end; state->current_phase++) {

        make_timestamp_root(&cur_p4est_timing_frame->iteration_start);



        // handle dynamic resources
        bool changed = state->is_new;
        int size_before, size_after;
        MPI_Comm_size(state->mpicomm, &size_before);

        if(rank == 0)
            printf("start of phase %d ---------> %d processes\n", state->current_phase, size_before);

        int status;
        if(state->current_phase != 1 && !state->is_new && check_rc){
            //if(rank == 0)
            //    printf("start of check_resource change\n");
            make_timestamp_root(&cur_p4est_timing_frame->rc_start);
            status = check_resource_change(state);
            make_timestamp_root(&cur_p4est_timing_frame->rc_end);
            //if(rank == 0)
            //    printf("end of check_resource change: %d\n", status);

            changed = status == MPI_SUCCESS;

            if(state->p4est == NULL)
                break;
        }
        printf("Rank %d finished res change\n", rank);
        MPI_Comm_size(state->mpicomm, &size_after);
        if(changed) {
            if(state->is_new || size_after > size_before){
                //if(rank == 0)
                //    printf("replace communicator: add\n");
                MPI_Comm dup_comm = MPI_COMM_NULL;
                MPI_Comm_dup(state->mpicomm, &dup_comm);
                state->p4est = p4est_dynres_replace_ext(state->p4est, dup_comm, 
                        -1, -1, -1, -1, (void *) state->scenario, NULL);
            
                make_timestamp_root(&cur_p4est_timing_frame->handle_rc_end);
                //if(rank == 0)
                //    printf("replaced communicator: add\n");
            }

            eval_mode(state->p4est->mpisize, mode_num);
        }

        //if(rank == 0)
        //    printf("partition\n");
        make_timestamp_root(&cur_p4est_timing_frame->partition_start);
        p4est_partition(state->p4est, 0, NULL);


        // simulate until next timestep
        //if(rank == 0)
        //    printf("simlate interval start\n");

        make_timestamp_root(&cur_p4est_timing_frame->interval_start);
        float new_time = simulate_interval(state);
        make_timestamp_root(&cur_p4est_timing_frame->interval_end);
        //if(rank == 0)
        //   printf("interval end\n");
        // update simulation time with time step width.
        state->t = new_time;
        l_iterations++;

        state->is_new = false;

        make_timestamp_root(&cur_p4est_timing_frame->iteration_end);
        
        if(rank == 0){
            cur_p4est_timing_frame->res_change_id = timings_cur_res_change_id;
            cur_p4est_timing_frame->delta_procs = size_after - size_before; 
            cur_p4est_timing_frame->num_procs = state->p4est->mpisize;
            cur_p4est_timing_frame->rc_status = timings_cur_rc_status;
            cur_p4est_timing_frame->iteration = state->current_phase;
            
            if(state->current_phase <= state->num_phases && state->t <= state->t_end){
                init_add_timing(p4est_timing_list, (void**) &cur_p4est_timing_frame, sizeof(timing_frame_p4est_t));
            }
        }
        

    }

    make_timestamp_base(&cur_base_timing_frame->sim_end);

    if(rank == 0)
        printf("FINISHED SIMULATION.\nWriting timings to file...\n");

    char *cwd_res = getcwd(cwd, sizeof(cwd));

    if(rank == 0){

        char filename_root[256];
        char *path_base = getenv("DYNMPI_BASE");
        sprintf(filename_root, "%s/build/p4est_dynres/applications/output/p4est/timings_p4est_%s_root.csv",path_base, mode_str);
        print_list_to_file_p4est(p4est_timing_list, filename_root);
    }

    return 0;
}


int main(int argc, char **argv) {

    base_timing_list = (node_t *)calloc(1, sizeof(node_t));
    p4est_timing_list = (node_t *)calloc(1, sizeof(node_t));

    init_add_timing(base_timing_list, (void**) &cur_base_timing_frame, sizeof(timing_frame_base_t));
    make_timestamp_base(&cur_base_timing_frame->app_start);

    simulation_state state = initialize(argc, argv);
    
    simulate(&state);
    make_timestamp_base(&cur_base_timing_frame->finalize_start);

    MPI_Session_finalize(&state.session);
    make_timestamp_base(&cur_base_timing_frame->finalize_end);

    char filename_base[256];
    char *path_base = getenv("DYNMPI_BASE");
    sprintf(filename_base, "%s/build/p4est_dynres/applications/output/p4est/timings_p4est_%s_base_%d:%d.csv", path_base, mode_str, cur_base_timing_frame->res_change_id, timings_my_rank);    
    print_list_to_file_base(base_timing_list, filename_base);

    _exit(0);
}


/**
 * @brief main program for the simulation on a p4est
 *
 * Initializes new processes and starts the simulation
 *
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
    snprintf(buf, 0x20, "%d", 45);
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
*/
