#include "logging.h"
#include "math.h"
#include "scheduler_datatypes.h"
#include "scheduler_mgmt.h"
#include "util.h"

#define DEFAULT_CHANGE_PROB 1.0

struct user_request_manager {
  MPIDYNRES_scheduler *scheduler;
  int num_processes;
};
typedef struct user_request_manager user_request_manager;


/**
 * @brief      Initialize the manager
 *
 * @param      scheduler The scheduler that is using the management interface
 *
 * @return     The new manager object
 */
MPIDYNRES_manager MPIDYNRES_manager_init(MPIDYNRES_scheduler *scheduler) {
  user_request_manager *res;
  res = calloc(1, sizeof(user_request_manager));
  if (res == NULL) {
    die("Memory Error\n");
  }
  res->scheduler = scheduler;
  res->num_processes = scheduler->num_scheduling_processes;

  return res;
}



/**
 * @brief      Free a manger
 *
 * @param      manager The manager to be freed
 *
 * @return     if != 0, an error occured
 */
int MPIDYNRES_manager_free(MPIDYNRES_manager manager) {
  free(manager);
  return 0;
}


/*
 * No support yet
 */
int MPIDYNRES_manager_register_scheduling_hints(MPIDYNRES_manager manager,
                                                int src_process_id,
                                                MPI_Info scheduling_hints,
                                                MPI_Info *o_answer) {
  (void) manager;
  (void) src_process_id;
  (void) scheduling_hints;
  *o_answer = MPI_INFO_NULL;
  return 0;
}


/**
 * @brief      Get initial process set
 *
 * @details    Check for config or default to cr id 1
 *
 * @param      manager The manager used
 *
 * @param      o_initial_pset The initial pset is returned here
 *
 * @return     if != 0, an error occured
 */
int MPIDYNRES_manager_get_initial_pset(MPIDYNRES_manager manager,
                                       set_int *o_initial_pset) {
  user_request_manager *mgr = (user_request_manager *)manager;
  int vlen;
  int in_there;
  int num_init = 1;
  MPI_Info config = mgr->scheduler->config->manager_config;

  if (config == MPI_INFO_NULL) {
    in_there = false;
  } else {
    MPI_Info_get_valuelen(config, "manager_initial_number", &vlen, &in_there);
  }
  if (in_there) {
    char *value = calloc(vlen+1, sizeof(char));
    MPI_Info_get(config, "manager_initial_number", vlen, value, &in_there);
    num_init = atoi(value);
    free(value);
    if (num_init < 1 || num_init > mgr->num_processes) {
      die("Key manager_initial_number is invalid.\n");
    }
  } else {
      num_init = 1;
  }

  // start processes 1...num_init
  *o_initial_pset = set_int_init(int_compare);
  for (int i = 0; i < num_init; i++) {
    set_int_insert(o_initial_pset, i+1);
  }

  return 0;
}

/**
 * @brief      Generate a cr set of specific size either part or not part
 * of current running crs
 *
 * @param      scheduler the scheduler
 *
 * @param      pointer where the new set will be returned
 *
 * @param      size the number of elements in the set
 *
 */
static void gen_set(MPIDYNRES_scheduler *scheduler, set_int *set, size_t delta, MPIDYNRES_RC_type type) {
  *set = set_int_init(int_compare);

  if(type == MPIDYNRES_RC_ADD) {
      for (int i = 0; i < delta; i++) {
          set_int_insert(set, i+scheduler->running_crs.size+1);
      }
  } else {
      for(int i = 0; i < delta; i++) {
          set_int_insert(set, scheduler->running_crs.size-i);
      }
  }
}


/**
 * @brief      Handle a resource change query, don't use with this manager
 *
 * @param      manager The manager used
 *
 * @param      src_process_id The cr id of the calling computing resource
 *
 * @param      o_rc_info The resource change info that should be returned to the application
 *
 * @param      o_rc_type The type of resource change is returned here
 *
 * @param      o_new_pset The new process set is returned here
 *
 * @return     if != 0, an error occured
 */
int MPIDYNRES_manager_handle_rc_msg(MPIDYNRES_manager manager,
                                    int src_process_id, MPI_Info *o_rc_info,
                                    MPIDYNRES_RC_type *o_rc_type,
                                    set_int *o_new_pset) {
  return 1;
}


/**
 * @brief      Handle a resource change query posed by the user
 *
 * @param      manager The manager used
 *
 * @param      src_process_id The cr id of the calling computing resource
 *
 * @param      rc_type The type of resource change
 *
 * @param      delta The size of the delta pset
 *
 * @param      o_rc_info The resource change info that should be returned to the application
 *
 * @param      o_new_pset The new process set is returned here
 *
 * @return     if != 0, an error occured
 */
int DEBUG_MPIDYNRES_manager_handle_user_rc_msg(MPIDYNRES_manager manager,
                                    int src_process_id, 
                                    int delta,
                                    MPI_Info *o_rc_info,
                                    set_int *o_new_pset) {
  (void)src_process_id;
  user_request_manager *mgr = (user_request_manager *)manager;
  int num_running = mgr->scheduler->running_crs.size;
  int num_free = mgr->scheduler->num_scheduling_processes - num_running;
  printf("user_request_manager: available for this manager\n");
  if (delta < 0) {
      assert(num_running >= -delta);
      printf("num running: %d but delta is %d\n", num_running, delta);
      gen_set(mgr->scheduler, o_new_pset, -delta, MPIDYNRES_RC_SUB);
  } else if (delta > 0) {
      assert(num_free >= delta);
      printf("num free: %d but delta is %d\n", num_free, delta);
      gen_set(mgr->scheduler, o_new_pset, delta, MPIDYNRES_RC_ADD);
  }

  return 0;
}


