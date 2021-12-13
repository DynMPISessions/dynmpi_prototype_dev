/*
 * Functions defined in this file can be called from the simulated application
 * They are used for interaction with the simulated runtime (the resource
 * manager)
 */
#ifndef MPIDYNRES_H
#define MPIDYNRES_H

#include <limits.h>
#include <mpi.h>
#include <stdbool.h>
#include <stddef.h>

#define MPI_MAX_PSET_NAME_LEN MPI_MAX_INFO_KEY

#define MPIDYNRES_INVALID_SESSION_ID INT_MAX

#define MPIDYNRES_NO_ORIGIN_RC_TAG -1


struct internal_MPIDYNRES_Session {
  int session_id;
  MPI_Info info;
};
typedef struct internal_MPIDYNRES_Session internal_MPIDYNRES_Session;

typedef internal_MPIDYNRES_Session *MPIDYNRES_Session;

extern MPIDYNRES_Session MPIDYNRES_Session_NULL;



/*
 * MPI Draft API
 */
int MPIDYNRES_Session_init(MPI_Info info, MPI_Errhandler errhandler,
                     MPIDYNRES_Session *session);

int MPIDYNRES_Session_finalize(MPIDYNRES_Session *session);

int MPIDYNRES_Session_get_info(MPIDYNRES_Session session, MPI_Info *info_used);

/*
 * Modified Draft API
 */
int MPIDYNRES_Session_get_psets(MPIDYNRES_Session session, MPI_Info info, MPI_Info *psets);

int MPIDYNRES_Session_get_pset_info(MPIDYNRES_Session session, char const *pset_name,
                              MPI_Info *info);

int MPIDYNRES_Group_from_session_pset(MPIDYNRES_Session session, const char *pset_name,
                                MPI_Group *newgroup);

int MPIDYNRES_Comm_create_from_group(MPI_Group group, const char *stringtag, MPI_Info info, MPI_Errhandler errhandler, MPI_Comm *newcomm);

/*
 * pset Management
 */
/**
 * @brief      different set operations possible on pset
 */
enum MPIDYNRES_pset_op {
  MPIDYNRES_PSET_UNION,
  MPIDYNRES_PSET_INTERSECT,
  MPIDYNRES_PSET_DIFFERENCE,
};
typedef enum MPIDYNRES_pset_op MPIDYNRES_pset_op;

/*
 * Create a new URI based on (valid) URIs and a set operation
 */
int MPIDYNRES_pset_create_op(MPIDYNRES_Session session,
                             MPI_Info hints,
                             char const pset1[],
                             char const pset2[],
                             MPIDYNRES_pset_op op,
                             char pset_result[MPI_MAX_PSET_NAME_LEN]);

/*
 * Mark pset as free, if all processes in the pset have marked it as free or have exited, it will be deleted
 */
int MPIDYNRES_pset_free(MPIDYNRES_Session session, char i_pset_name[]);




/*
 * Query Runtime (Resource Manager) for Resource Changes (RCs)
 */
int MPIDYNRES_add_scheduling_hints(MPIDYNRES_Session session, MPI_Info hints, MPI_Info *answer);

/**
 * @brief      Different types of resource changes
 */
enum MPIDYNRES_RC_type {
  MPIDYNRES_RC_NONE,
  MPIDYNRES_RC_ADD,
  MPIDYNRES_RC_SUB,
};
typedef enum MPIDYNRES_RC_type MPIDYNRES_RC_type;

typedef int MPIDYNRES_RC_tag;

/*
 * debug function for requesting specific RCs from the user application
 */
int DEBUG_MPIDYNRES_RC_user_request(MPIDYNRES_Session session, 
                     int delta,
                     char delta_pset[MPI_MAX_PSET_NAME_LEN],
                     MPIDYNRES_RC_type *rc_type,
                     MPIDYNRES_RC_tag *tag, MPI_Info *info);


int MPIDYNRES_RC_get(MPIDYNRES_Session session,
                     MPIDYNRES_RC_type *rc_type,
                     char delta_pset[MPI_MAX_PSET_NAME_LEN],
                     MPIDYNRES_RC_tag *tag,
                     MPI_Info *info);

/*
 * Accept runtime change and provide info that will be added to the new pset
 * processes
 */
int MPIDYNRES_RC_accept(MPIDYNRES_Session session, MPIDYNRES_RC_tag tag,
                        MPI_Info info);


/*
 * Exit from simulated process (real process will idle until all exit or its
 * started again)
 */
void MPIDYNRES_exit();



/*
 * Utility function to convert an array of strings (key, value alternating)
 * into an MPI_Info object
 * kvlist_size is the size of the kvlist array, NOT the number of keys (which is kvlist_size / 2)
 */
int MPIDYNRES_Info_create_strings(size_t kvlist_size, char const * const kvlist[], MPI_Info *info);

#endif
