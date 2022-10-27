/*
 * Copyright (c) 2004-2006 The Trustees of Indiana University and Indiana
 *                         University Research and Technology
 *                         Corporation.  All rights reserved.
 * Copyright (c) 2006      Cisco Systems, Inc.  All rights reserved.
 *
 * Sample MPI "hello world" application in C
 */

#include <stdio.h>
#include <stdlib.h>
#include <linux/limits.h>
#include <time.h>
#include <stdbool.h>
#include "mpi.h"
#include <unistd.h>
#include "tools/bench_timings.h"

#define noop

/* Simulation parameters */
const unsigned long long PROBLEM_SIZE = 100000000;
int ITER_MAX = 200;
int num_procs;
int my_work_rank;
unsigned long long start_index;
unsigned long long end_index;
int rank; 
bool spawned = false;
float result_last;



extern char *optarg;
extern int optind, opterr, optout, optopt;

/* Timing */
node_t *base_timing_list;
node_t *synth_timing_list;

timing_frame_base_t *cur_base_timing_frame;
timing_frame_synth_t *cur_synth_timing_frame;

#define BENCH_INCRIMENTAL 0
#define BENCH_SEQUENTIAL 1
#define BENCH_MIXED 2

static int n_start = 0;
static int rc_id = 0;
static int recvd_delta = -1;
static int num_delta = 1;
static int cur_num_delta = 1;
static int proc_limit = -1;
static int proc_limit_down = -1;
static int rc_frequency = 0;
static int iter_since_change = 0;
static int blocking = 0;
static bool check_rc = true;

static char mode[2];
static char mode_str[64];
static int mode_num = BENCH_INCRIMENTAL;
static int mode_type = MPI_RC_ADD;
static int cur_type = MPI_RC_ADD;

int set_mode(const char *c_mode){
    
    if(strlen(c_mode) != 2){
        return -1;
    }
    //printf("mode: %c%c\n", c_mode[0], c_mode[1]);
    mode[0] = c_mode[0];
    mode[1] = c_mode[1];
    if(c_mode[0] == 'i'){
        mode_num = BENCH_INCRIMENTAL;
        cur_num_delta = num_delta;
        //printf("cur_num_delta: %d\n", cur_num_delta);
    }else if(c_mode[0] == 's'){
        mode_num = BENCH_SEQUENTIAL;
        cur_num_delta = num_delta;
        //printf("cur_num_delta: %d\n", cur_num_delta);
    }else if(c_mode[0] == 'm'){
        mode_num = BENCH_MIXED;
        cur_num_delta = proc_limit / num_delta * num_delta;
        //printf("curnum delta : %d\n", cur_num_delta);
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
            check_rc = mode_type == MPI_RC_ADD ? current_size + cur_num_delta <= proc_limit : (current_size - cur_num_delta >= proc_limit);
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
                //printf("mode_num %d, cur_type %d, current_size %d, cur_num_delta %d\n", mode_num, cur_type, current_size, cur_num_delta);
                check_rc = cur_num_delta <= 0 ? false : cur_type == MPI_RC_ADD ? current_size + cur_num_delta <= proc_limit : (current_size - cur_num_delta >= proc_limit_down);
            }else{
                check_rc = mode_type == MPI_RC_ADD ? current_size + cur_num_delta <= proc_limit : (current_size - cur_num_delta >= proc_limit && current_size - cur_num_delta >= 1);

            }

            //printf("cur num delta init %d, check_rc %d, cur type %d, current size %d\n", cur_num_delta, check_rc, cur_type, current_size);
            break;            
    }
    iter_since_change = 0;
}

int parse_argument(int argc, char** argv){

    opterr = 0;

    int c;
    char *mode_s = NULL;

    while ((c = getopt (argc, (char**)argv, "c:m:l:n:f:b:")) != -1){
      switch (c)
        {
        case 'b':
            //blocking = 1;
            blocking = atoi(optarg);
            break;
        case 'c':
            ITER_MAX = atoi(optarg);
            break;
        case 'l':
            proc_limit = atoi(optarg);
            break;
        case 'n':
            num_delta = atoi(optarg);
            break;
        case 'f':
            rc_frequency = atoi(optarg);
            break;
        case 'm':
            mode_s = optarg;
            break;
        case '?':
          if (optopt == 'c' || optopt == 'm' || optopt == 'l' || optopt == 'n' || optopt == 'f' || optopt == 'b')
            fprintf (stderr, "Option -%c requires an argument.\n", optopt);
          else
            fprintf (stderr,
                     "Unknown option character `\\x%x'.\n",
                     optopt);
          return 1;
        default:
            printf("Bad Parameter: Abort!\n");
          abort ();
        }
    }
    if(NULL != mode_s){
        set_mode(mode_s);
    }
    return 0;
}


/*-------------------------------- HERE BEGINS THE MAIN FILE ------------------------------------------- */




/* create a communicator from the given process set */
#pragma region
int init(MPI_Session *session_handle, MPI_Comm *comm, char *main_pset_name){
    int rc;
    MPI_Group wgroup = MPI_GROUP_NULL;

    /* create a group from pset */
    rc= MPI_Group_from_session_pset (*session_handle, main_pset_name, &wgroup);
    MPI_Group_size(wgroup, &num_procs);
    /* create a communicator from group */
    //printf("creating communicator of size: %d\n", num_procs);
    if(MPI_SUCCESS != (rc = MPI_Comm_create_from_group(wgroup, "mpi.forum.example", MPI_INFO_NULL, MPI_ERRORS_RETURN, comm))){
        printf("MPI_comm_create_from_group failed\n");
        MPI_Session_finalize(session_handle);
        return -1;
    }

    /* new rank & size */
    MPI_Comm_size(*comm, &num_procs);
    MPI_Comm_rank(*comm, &my_work_rank);
    //printf("rank %d created communicator of size: %d\n", my_work_rank, num_procs);

    MPI_Group_free(&wgroup);
    return rc;
}

#pragma endregion
 

int resource_change_step(MPI_Session *session_handle, MPI_Comm *lib_comm, char *pset_name, char *delta_pset, char *pset_result, int *terminate){
    
    struct timespec ts;
    ts.tv_sec = 0;
    ts.tv_nsec = 1000000;

    int rc;
    int rc_type = MPI_RC_NULL;
    int new_rc_available;

    int incl_flag;
    int my_rank;
    MPI_Info info = MPI_INFO_NULL;
    int status = MPI_RC_INVALID;
    char dummy_pset_name[] = "test1";
    char prefname[MPI_MAX_PSET_NAME_LEN-1]; 

    //printf("_________________________resource change step\n");

    MPI_Comm_rank(*lib_comm, &my_rank);

    /*****************fetch resource change****************************/
    if(my_rank==0){
        MPI_Info_create(&info);
        MPI_Info_set(info, "MPI_RC_BOUND_PSET", pset_name);
        make_timestamp_root(&cur_synth_timing_frame->get_rc_start);
        rc = MPI_Session_get_res_change(*session_handle, &rc_type, delta_pset, &incl_flag, &status, &info);

        if(rc_type == MPI_RC_NULL || status == MPI_RC_FINALIZED){

            // if no rc and not requested: request a resource change: add one process
            MPI_Session_request_res_change(*session_handle, cur_num_delta, pset_name, cur_type, &info);
            // As the request returns when notification is send, do an ugly wait until the rc is defined 
            
            while(rc_type == MPI_RC_NULL){
                rc = MPI_Session_get_res_change(*session_handle, &rc_type, delta_pset, &incl_flag, &status, &info);
                if(rc_type == MPI_RC_NULL){
                    nanosleep(&ts, NULL);
                }
            }
        }

        if(rc_type == MPI_RC_NULL || (status != MPI_RC_ANNOUNCED && status != MPI_RC_CONFIRMATION_PENDING)){
            rc = MPI_ERR_OTHER;
        }
        MPI_Info_free(&info);
        cur_synth_timing_frame->rc_status = status;
        set_res_change_id(&cur_synth_timing_frame->res_change_id, delta_pset);
    }
    /*****************end of fetch resource change**********************/    

    /******************handle resource change***************************/
    if(my_rank == 0 && rc_type != MPI_RC_NULL && status == MPI_RC_ANNOUNCED){ 
        /* create a new pset name */ 
        strcpy(prefname,dummy_pset_name);
        strcat(prefname, delta_pset);
        
        MPI_Info_create(&info);
        MPI_Info_set(info, "MPI_PSETOP_PREF_NAME", prefname);

        make_timestamp_root(&cur_synth_timing_frame->psetop_start);
        if(rc_type ==  MPI_RC_ADD){
            
            MPI_Session_pset_create_op(*session_handle, MPI_PSETOP_UNION, pset_name, delta_pset, pset_result, &info);
            
        }else if(rc_type == MPI_RC_SUB){
            MPI_Session_pset_create_op(*session_handle, MPI_PSETOP_DIFFERENCE, pset_name, delta_pset, pset_result, &info);
        }

        MPI_Info_free(&info);
    }
    //if(my_rank != 0)
    //    sleep(1);
    /***************** end of handle resource change *******************/ 
    //printf("_________________________bcast\n");

    /* Root needs to tell other processes about the resource change */
    MPI_Bcast(&rc, 1, MPI_INT, 0, *lib_comm);
    if(MPI_ERR_OTHER == rc)return rc;
    //printf("init\n");
    //MPI_Comm_free(lib_comm);
    //init(session_handle, lib_comm, pset_name);
    //printf("init success\n");

    MPI_Bcast(pset_result, 511, MPI_CHAR, 0 , *lib_comm);

    
    
    /*********************** accept resource change *************************/

    if(my_rank == 0 && blocking){
        MPI_Info_create(&info);
        MPI_Info_set(info, "mpi_blocking", "1");
        //printf("blocking accept\n");
    }

    *terminate = 0;
    make_timestamp_root(&cur_synth_timing_frame->accept_start);
    rc = MPI_Session_accept_res_change(session_handle, &info, delta_pset, pset_result, 0, lib_comm, terminate);
    make_timestamp_root(&cur_synth_timing_frame->accept_end);
    /*********************** end of accept resource change *************************/

    if(my_rank == 0 && blocking)
        MPI_Info_free(&info);
    /* if successful create new communicator, else need to break (here indicated by MPI_COMM_NULL)*/
    if(MPI_SUCCESS == rc){
        if(0 == *terminate){
            make_timestamp_root(&cur_synth_timing_frame->reinit_start); 
            strcpy(pset_name, pset_result);
            rc = init(session_handle, lib_comm, pset_name);
        }else{
            return MPI_SUCCESS;
        }
    }
    
    return rc;

}


int work_step(){
    unsigned long long  n;
    int i = 0;
    double res = .0;
    //float result = static_cast <float> (rand())/(static_cast <float> (RAND_MAX));
    double a = 1.57, b = 20.1113, c = 5.4200102;
    double d = 5.223, e = 1.89, f = 22.223;
    double t, x, y;
    for(n = start_index; n < end_index; n++){
        t = (double) n;
        x = a*t*t + b*t + c;
        y = d*t*t + e*t + f;

        if(x > 42.42 && x < 42.4242 && y > 42.42 && y < 42.4242){
            //printf("You've found the answer to everything!\n");
            i++;
        }

        /*
        for(i=1; i<=100; i++){
            float factor = static_cast <float> (rand())/(static_cast <float> (RAND_MAX));
            if(i%2 == 0)result*=factor;
            else result/=factor;
        }
        */
    }
    return i;
    
}

void rebalance_step(){
    unsigned long long chunk_size=PROBLEM_SIZE/num_procs;
    start_index = chunk_size*my_work_rank;
    end_index   = (my_work_rank==num_procs-1)   ? PROBLEM_SIZE  : (my_work_rank+1)*chunk_size;
    //printf("    RANK %llu/%llu: rebalancing:  range[%llu,%llu]\n", my_work_rank,num_procs-1,start_index, end_index);
}

int main(int argc, char* argv[])
{
	char host[256];
	gethostname(host, 256);
    int  size, len, flag, npsets, counter=0;
    char pset_name[MPI_MAX_PSET_NAME_LEN-1];
    char app_pset_name[MPI_MAX_PSET_NAME_LEN-1];
    char delta_pset[MPI_MAX_PSET_NAME_LEN-1];
    char pset_result[MPI_MAX_PSET_NAME_LEN-1]={0};
    int incl_flag = 0;
    //system("rm -R /dev/shm/*");
    MPI_Session session_handle;
    MPI_Info info = MPI_INFO_NULL;
    MPI_Info info2 = MPI_INFO_NULL;
    MPI_Comm lib_comm = MPI_COMM_NULL;
 
    int rc_type;
    int rc_status;
    int rc = -16;

    base_timing_list = (node_t *)calloc(1, sizeof(node_t));
    synth_timing_list = (node_t *)calloc(1, sizeof(node_t));

    
    init_add_timing(base_timing_list, (void**) &cur_base_timing_frame, sizeof(timing_frame_base_t));
    make_timestamp_base(&cur_base_timing_frame->app_start);
    parse_argument(argc, argv);
    //printf("Read paramters: mode=%s(%d:%d), num_delta=%d, proc_limit=%d, rc_freq=%d\n", mode, mode_num, mode_type, num_delta, proc_limit, rc_frequency);

    strcpy(pset_name, "test1");

    /* initialize the session */
    MPI_Info_create(&info2);
    MPI_Info_set(info2, "MPI_RC_BOUND_PSET", "mpi://SELF");

    make_timestamp_base(&cur_base_timing_frame->init_start);
    int init_ret = MPI_Session_init(MPI_INFO_NULL, MPI_ERRORS_RETURN, &session_handle);
    make_timestamp_base(&cur_base_timing_frame->init_end);
    rank = 0; 
    //MPI_Session_get_num_psets(session_handle, info, &npsets);
    strcpy(pset_result,pset_name);
    /*
    if(ompi_proc_local_proc->super.proc_name.vpid==0){
        printf("pset_op\n");
        MPI_Session_pset_create_op(session_handle, MPI_PSETOP_UNION, pset_name, pset_name, pset_result);
    }
    */
    /* check if there is a resource change right at the beginning */
    int ret = MPI_Session_get_res_change(session_handle, &rc_type, delta_pset, &incl_flag, &rc_status, &info2);
    
    /* if we are included in the delta_pset we are a dynamically added process, so we need to confirm the resource change */
    if(rc_type != MPI_RC_NULL && rc_status == MPI_RC_ANNOUNCED && incl_flag){
        //printf("    DELTA PSET RANK %d: I was added dynamically. Need to confirm \n", rank);
        spawned = true;
        
        char *pset_name_fresh = (char*) malloc(MPI_MAX_PSET_NAME_LEN);

        set_res_change_id(&cur_base_timing_frame->res_change_id, delta_pset);
        make_timestamp_base(&cur_base_timing_frame->confirm_start);
        //printf("Starting CONFIRM\n");
        rc = MPI_Session_confirm_res_change(&session_handle, &info, delta_pset, &pset_name_fresh);

        make_timestamp_base(&cur_base_timing_frame->confirm_end);
        
        strcpy(pset_name, pset_name_fresh);
        if(MPI_SUCCESS == rc){
            printf("    DELTA PSET RANK %d: Confirmation succeeded. Communication will happen via pset: %s\n", rank, pset_name);
        }else{
            printf("    DELTA PSET RANK %d: Confirmation failed. ERROR: %d\n", rank, rc);
        }
        init(&session_handle, &lib_comm, pset_name);
        //printf("Rank %d: start recv counter\n", my_work_rank);

        /* get the current iteration from the main rank */
        MPI_Bcast(&counter, 1, MPI_INT, 0, lib_comm);
        if(mode_num == BENCH_MIXED){
            MPI_Bcast(&recvd_delta, 1, MPI_INT, 0, lib_comm);
        }
        //printf("    RANK %d: received counter value %d from root in new communicator\n", my_work_rank, counter);
        if(mode_num == BENCH_SEQUENTIAL){
            MPI_Bcast(&cur_num_delta, 1, MPI_INT, 0, lib_comm);
        }
        

    }else{
        
        /* initialize communication */
        init(&session_handle, &lib_comm, pset_name);
        //printf("init end\n");
        rc = MPI_SUCCESS;
        if(my_work_rank == 0){
        	n_start = num_procs;
        }
    }
    timings_my_rank = my_work_rank;

    if(0 == strcmp(mode, "s_")){
        proc_limit_down = proc_limit;
        proc_limit = num_procs;
        //printf("setting limits for s_ mode: up: %d, down: %d\n", proc_limit, proc_limit_down);
    }

    eval_mode(num_procs, -1);
    
    int res;
    make_timestamp_base(&cur_base_timing_frame->sim_start);
    /********************** START OF MAIN LOOP *******************************/
    while(counter++< ITER_MAX){
        if(my_work_rank == 0){
            printf("\n START OF ITERATION: %d with %d processes\n", counter, num_procs);
            init_add_timing(synth_timing_list, (void**) &cur_synth_timing_frame, sizeof(timing_frame_synth_t));
            cur_synth_timing_frame->iteration = counter;
            cur_synth_timing_frame->num_procs = num_procs;
        }
        make_timestamp_root(&cur_synth_timing_frame->rebalance_start);
        
        /* Rebalance */
        if(rc == MPI_SUCCESS){
            rebalance_step();
        }
        /* Work step */
        //printf("rank %d: first barrier\n", my_work_rank);
        MPI_Barrier(lib_comm);
        make_timestamp_root(&cur_synth_timing_frame->work_start);
        res += work_step();
        make_timestamp_root(&cur_synth_timing_frame->work_end);
        //printf("rank %d: second barrier\n", my_work_rank);
        MPI_Barrier(lib_comm);
        MPI_Bcast(&cur_num_delta, 1, MPI_INT, 0, lib_comm);

        /* Resource Change step */
        make_timestamp_root(&cur_synth_timing_frame->rc_start);

        if(check_rc && iter_since_change++ >= rc_frequency){
            //printf("%d, %d, num_procs - cur_num_delta = %d\n", num_procs, cur_num_delta, num_procs - cur_num_delta);
            //if(cur_type == MPI_RC_ADD)printf("**************ADD");
            //else printf("*****************SUB\n");
            int size_before = num_procs;
            int terminate = 0;
            rc = resource_change_step(&session_handle, &lib_comm, pset_name, delta_pset, pset_result, &terminate);

            if(my_work_rank == 0)
                cur_synth_timing_frame->delta_procs = num_procs - size_before;


            make_timestamp_root(&cur_synth_timing_frame->rc_end);

            /* Data redistribution */
            if(rc == MPI_ERR_PENDING){
                noop;
                //printf("    RANK %d: Resource change pending. Communication will still happen via pset: %s\n", my_work_rank, pset_name);
            }else if (rc == MPI_SUCCESS) {
                /* terminate substracted processes */
                if(terminate){
                    printf("    Old RANK %d: Resource change succeeded. I am not needed anymore. Goodbye!\n", my_work_rank);
                    break;
                }

                //printf("    RANK %d: Resource change succeeded. Communication will now happen via pset: %s\n", my_work_rank, pset_name);
                eval_mode(num_procs, mode_num);
                
                MPI_Bcast(&counter, 1, MPI_INT, 0, lib_comm);
                if(mode_num == BENCH_MIXED || mode_num == BENCH_SEQUENTIAL){
                    MPI_Bcast(&cur_num_delta, 1, MPI_INT, 0, lib_comm);
                }
            }
            

        }
        make_timestamp_root(&cur_synth_timing_frame->iteration_end);

        if(my_work_rank == 0 && counter == ITER_MAX){
            printf("finished benchmark successfully!\n"); 
        }
    }
    /************************ END OF MAIN LOOP ************************************/
    make_timestamp_base(&cur_base_timing_frame->sim_end);

    if(MPI_INFO_NULL != info){
        MPI_Info_free(&info);
    }
    if(res > ITER_MAX*100){
        //printf("stay positive!\n");
        iter_since_change = 0;
    }
    make_timestamp_base(&cur_base_timing_frame->finalize_start);
    MPI_Session_finalize(&session_handle);
    make_timestamp_base(&cur_base_timing_frame->finalize_end);
    //printf("rank%d: session_finalized\n", my_work_rank);

    /* print the timings to file */
    char *path_base = getenv("DYNMPI_BASE");


    char filename_base[256];
    sprintf(filename_base, "%s/build/p4est_dynres/applications/output/synth/timings_p4est_%s_l%d_n%d_c%d_f%d_base_%d:%d.csv", path_base, mode_str, proc_limit, num_delta, ITER_MAX, rc_frequency, cur_base_timing_frame->res_change_id, timings_my_rank);    

    print_list_to_file_base(base_timing_list, filename_base);
    free(base_timing_list);
    
    if(my_work_rank == 0){
        char filename_root[256];
        sprintf(filename_root, "%s/build/p4est_dynres/applications/output/synth/timings_p4est_%s_np%d_l%d_n%d_c%d_f%d_root.csv",path_base, mode_str, n_start, proc_limit, num_delta, ITER_MAX, rc_frequency);
        print_list_to_file_synth(synth_timing_list, filename_root);
        free(synth_timing_list);
    }

    //free(path_base);
    printf("rank%d terminating\n", my_work_rank);
    _exit(0);
    //return 0;
}
