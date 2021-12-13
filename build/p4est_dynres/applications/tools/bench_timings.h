#include <unistd.h>
#include <stdio.h>
#include <sys/time.h>
#include <stddef.h>
#include <string.h>
#include <errno.h>


static int timings_my_rank = -1;
static int timings_cur_res_change_id = 0;
static int timings_cur_rc_status = 0;


/*
 * LIST FUNCTIONS
 */
typedef struct node {
    void *val;
    struct node * next;
} node_t;

void push(node_t * head, void *val) {
    node_t * current = head;
    while (current->next != 0) {
        current = current->next;
    }

    /* now we can add a new variable */
    current->next = (node_t *) malloc(sizeof(node_t));
    current->next->val = val;
    current->next->next = 0;
}

int pop(node_t ** head) {
    int retval = -1;
    node_t * next_node = 0;

    if (*head == 0) {
        return -1;
    }

    next_node = (*head)->next;
    free((*head)->val);
    free(*head);
    *head = next_node;
    return retval;
}

/*
 ************* INIT & TIMESTAMP ****************      
 */

void init_add_timing(node_t * list, void ** timing, size_t frame_size){
    *timing = calloc(1, frame_size);
    push(list, (void*) *timing);
}

void make_timestamp_base(long * ptr_to_dest){
    struct timeval timestamp;
    gettimeofday(&timestamp, 0);

    *ptr_to_dest = timestamp.tv_sec * 1000000 + timestamp.tv_usec;
}

void make_timestamp_root(long * ptr_to_dest){
    if(0 != timings_my_rank)return;

    struct timeval timestamp;
    gettimeofday(&timestamp, 0);

    *ptr_to_dest = timestamp.tv_sec * 1000000 + timestamp.tv_usec;
}

void make_timestamp_rootonly(long * ptr_to_dest, int rank){

    struct timeval timestamp;
    gettimeofday(&timestamp, 0);

    *ptr_to_dest = timestamp.tv_sec * 1000000 + timestamp.tv_usec;
}

void set_res_change_id(int *dest, const char *delta_pset){
    *dest = atoi(&delta_pset[strlen(delta_pset)-1]) + 1;
    timings_cur_res_change_id = *dest;
}

/*
 * *********** WRITE OUPUT ********************
 */


/* Base frame */
typedef struct{
    long app_start;
    long init_start;
    long init_end;
    long confirm_start;
    long confirm_end;
    long sim_start;
    long sim_end;
    long finalize_start;
    long finalize_end;
    int res_change_id;
}timing_frame_base_t;

void print_list_to_file_base(node_t * head,  const char *filename) {

    
    struct timeval start, end;
    gettimeofday(&start, 0);
    long l_start = start.tv_sec * 1000000 + start.tv_usec;
 
    FILE *f = fopen(filename, "w+");
    if(NULL == f){
        printf("Error opening file: %s\n", strerror(errno));
        return;
    }

    char header[] = "res_change_id, app_start, init_start, init_end, confirm_start, confirm_end, sim_start, sim_end, finalize_start, finalize_end, start_writing, end_writing";
    fprintf(f, "%s\n", header);

    node_t * current = head->next;
    while (current != 0) {
        timing_frame_base_t *timing = (timing_frame_base_t *) current->val;

        fprintf(f, "%d,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,",   timing->res_change_id, timing->app_start, timing->init_start, timing->init_end, timing->confirm_start, 
                                                timing->confirm_end, timing->sim_start, timing->sim_end, timing->finalize_start, timing->finalize_end);
        pop(&current);
    }

    
    gettimeofday(&end, 0);
    long l_end = end.tv_sec * 1000000 + end.tv_usec;

    fprintf(f, "%ld,%ld\n", l_start, l_end);

    fclose(f);


}



/* P4est frame */
typedef struct{
    long iteration_start;
    long rc_start;
    long get_rc_start_1;
    long request_start;
    long get_rc_start_2;
    long psetop_start;
    long handle_rc_start;
    long handle_rc_end;
    long accept_start;
    long accept_end;
    long rc_end;
    long partition_start;
    long interval_start;
    long interval_end;
    long iteration_end;
    int iteration;
    int res_change_id;
    int num_procs;
    int delta_procs;
    int rc_status;
}timing_frame_p4est_t;


void print_list_to_file_p4est(node_t * head, char *filename) {


    FILE *f = fopen(filename, "w+");
    


    char header[] = "iteration, res_change_id, num_procs, delta_procs, rc_status, rc_start, get_rc_start_1, request_start, get_rc_start_2, psetop_start, handle_rc_start, handle_rc_end, accept_start, accept_end, rc_end, partition_start, interval_start, interval_end, iteration_end";
    fprintf(f, "%s\n", header);
    
    node_t * current = head->next;
    while (current != 0) {
        timing_frame_p4est_t *timing = (timing_frame_p4est_t *) current->val;
        fprintf(f, "%d,%d,%d,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n",   timing->iteration, timing->res_change_id, timing->num_procs, timing->delta_procs, timing->rc_status, 
                                                timing->rc_start, timing->get_rc_start_1, timing->request_start, timing->get_rc_start_2, 
                                                timing->psetop_start, timing->handle_rc_start, timing->handle_rc_end, timing->accept_start, 
                                                timing->accept_end, timing->rc_end, timing->partition_start, timing->interval_start, 
                                                timing->interval_end, timing->iteration_end);
        pop(&current);
    }

    fclose(f);

}

/* Synth frame */
typedef struct{
    long iteration_start;
    long rc_start;
    long get_rc_start;
    long psetop_start;
    long accept_start;
    long accept_end;
    long reinit_start;
    long rc_end;
    long rebalance_start;
    long work_start;
    long work_end;
    long iteration_end;
    int iteration;
    int res_change_id;
    int num_procs;
    int delta_procs;
    int rc_status;
}timing_frame_synth_t;


void print_list_to_file_synth(node_t * head, char *filename) {


    FILE *f = fopen(filename, "w+");
    char header[] = "iteration, res_change_id, num_procs, delta_procs, rc_status, rc_start, get_rc_start, psetop_start, accept_start, accept_end, reinit_start, rc_end, rebalance_start, work_start, work_end, iteration_end";
    fprintf(f, "%s\n", header);
    
    node_t * current = head->next;
    while (current != 0) {
        timing_frame_synth_t *timing = (timing_frame_synth_t *) current->val;
        fprintf(f, "%d,%d,%d,%d,%d,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld\n",   timing->iteration, timing->res_change_id, timing->num_procs, timing->delta_procs, timing->rc_status, 
                                                timing->rc_start, timing->get_rc_start, timing->psetop_start, timing->accept_start, timing->accept_end, timing->reinit_start, 
                                                timing->rc_end, timing->rebalance_start, timing->work_start, timing->work_end, timing->iteration_end);
        pop(&current);
    }

    fclose(f);
}
