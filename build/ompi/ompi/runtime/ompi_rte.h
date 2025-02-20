/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil -*- */
/*
 * Copyright (c) 2012-2015 Los Alamos National Security, LLC.  All rights reserved.
 * Copyright (c) 2013      Mellanox Technologies, Inc.
 *                         All rights reserved.
 * Copyright (c) 2014-2020 Intel, Inc.  All rights reserved.
 * Copyright (c) 2019      Research Organization for Information Science
 *                         and Technology (RIST).  All rights reserved.
 * Copyright (c) 2020      Triad National Security, LLC. All rights
 *                         reserved.
 * Copyright (c) 2020      Amazon.com, Inc. or its affiliates.  All Rights
 *                         reserved.
 *
 * Copyright (c) 2021      Nanook Consulting.  All rights reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 */

/* This is the public RTE interface to the OMPI layer. Any RTE can be
 * connected to the OMPI layer by creating a new static component in
 * this framework, assigning it a priority and including a configure.m4
 * to define when it should be built.
 *
 * Each component must provide a number of types and functions that mimic
 * those provided by PRTE. These include (where flexibility exists, the
 * pmix data type is shown, but any compatible type is allowed. For example,
 * the jobid field in ompi_process_name_t could be any type of integer, but
 * cannot be a string):
 *
 * (a) Process name objects and operations
 *     1. Definitions for integral types ompi_jobid_t and ompi_vpid_t.
 *        The jobid must be unique for a given MPI_COMM_WORLD capable of
 *        connecting to another OMPI_COMM_WORLD and the vpid will be the
 *        process's rank in MPI_COMM_WORLD.
 *     2. ompi_process_name_t - a struct that must contain at least two integer-typed fields:
 *           a. ompi_jobid_t jobid
 *           b. ompi_vpid_t vpid
 *        Note that the structure can contain any number of fields beyond these
 *        two, so the process name struct for any particular RTE can be whatever
 *        is desired.
 *     3. OMPI_NAME_PRINT - a macro that prints a process name when given
 *        a pointer to ompi_process_name_t. The output format is to be
 *        a single string representing the name.  This function should
 *        be thread-safe for multiple threads to call simultaneously.
 *     4. OMPI_PROC_MY_NAME - a pointer to a global variable containing
 *        the ompi_process_name_t for this process. Typically, this is
 *        stored as a field in the ompi_process_info_t struct, but that
 *        is not a requirement.
 *     5. OMPI_NAME_WIlDCARD - a wildcard name.
 *     6. ompi_rte_compare_name_fields - a function used to compare fields
 *        in the ompi_process_name_t struct. The function prototype must be
 *        of the form:
 *        int ompi_rte_compare_name_fields(ompi_rte_cmp_bitmask_t mask,
 *                                         ompi_process_name_t *name1,
 *                                         ompi_process_name_t *name2);
 *        The bitmask must be defined to indicate the fields to be used
 *        in the comparison. Fields not included in the mask must be ignored.
 *        Supported bitmask values must include:
 *           b. OMPI_RTE_CMP_JOBID
 *           c. OMPI_RTE_CMP_VPID
 *           d. OMPI_RTE_CMP_ALL
 *      7. uint64_t ompi_rte_hash_name(name) - return a string hash uniquely
 *         representing the ompi_process_name passed in.
 *      8. OMPI_NAME - an Opal DSS constant for a handler already registered
 *         to serialize/deserialize an ompi_process_name_t structure.
 *
 * (b) Collective objects and operations
 *     1. ompi_rte_collective_t - an OPAL object used during RTE collective operations
 *        such as modex and barrier. It must be an opal_list_item_t and contain the
 *        following fields:
 *           a. id (pmix type: int32_t)
 *           b. bool active
 *              flag that user can poll on to know when collective
 *              has completed - set to false just prior to
 *              calling user callback function, if provided
 *     2. ompi_rte_modex - a function that performs an exchange of endpoint information
 *        to wireup the MPI transports. The function prototype must be of the form:
 *        int ompi_rte_modex(ompi_rte_collective_t *coll);
 *        At the completion of the modex operation, the coll->active flag must be set
 *        to false, and the endpoint information must be stored in the modex database.
 *        This function must have barrier semantics across the MPI_COMM_WORLD of the
 *        calling process.
 *     3. ompi_rte_barrier - a function that performs a barrier operation within the
 *        RTE. The function prototype must be of the form:
 *        int ompi_rte_barrier(ompi_rte_collective_t *coll);
 *        At the completion of the barrier operation, the coll->active flag must be set
 *        to false
 *
 * (c) Process info struct
 *     1. ompi_process_info_t - a struct containing info about the current process.
 *        The struct must contain at least the following fields:
 *           a. app_num -
 *           b. pid - this process's pid.  Should be same as getpid().
 *           c. num_procs - Number of processes in this job (ie, MCW)
 *           d. my_node_rank - relative rank on local node to other peers this run-time
 *                    instance knows about.  If doing dynamics, this may be something
 *                    different than my_local_rank, but will be my_local_rank in a
 *                    static job.
 *           d. my_local_rank - relative rank on local node with other peers in this job (ie, MCW)
 *           e. num_local_peers - Number of local peers (peers in MCW on your node)
 *           f. my_hnp_uri -
 *           g. peer_modex - a collective id for the modex operation
 *           h. peer_init_barrier - a collective id for the barrier during MPI_Init
 *           i. peer_fini_barrier - a collective id for the barrier during MPI_Finalize
 *           j. job_session_dir -
 *           k. proc_session_dir -
 *           l. nodename - a string representation for the name of the node this
 *              process is on
 *           m. cpuset -
 *     2. ompi_process_info - a global instance of the ompi_process_t structure.
 *     3. ompi_rte_proc_is_bound - global boolean that will be true if the runtime bound
 *        the process to a particular core or set of cores and is false otherwise.
 *
 * (d) Error handling objects and operations
 *     1. void ompi_rte_abort(int err_code, char *fmt, ...) - Abort the current
 *        process with the specified error code and message.
 *     2. int ompi_rte_abort_peers(ompi_process_name_t *procs, size_t nprocs) -
 *        Abort the specified list of peers
 *     3. OMPI_ERROR_LOG(rc) - print error message regarding the given return code
 *
 * (e) Init and finalize objects and operations
 *     1. ompi_rte_init - a function to initialize the RTE. The function
 *        prototype must be of the form:
 *        int ompi_rte_init(int *argc, char ***argv);
 *     2. ompi_rte_finalize - a function to finalize the RTE. The function
 *        prototype must be of the form:
 *        int ompi_rte_finalize(void);
 *     3. void ompi_rte_wait_for_debugger(void) - Called during MPI_Init, this
 *        function is used to wait for debuggers to do their pre-MPI attach.
 *        If there is no attached debugger, this function will not block.
 *
 * (f) Database operations
 *     1. ompi_rte_db_store - a function to store modex and other data in
 *        a local database. The function is primarily used for storing modex
 *        data, but can be used for general purposes. The prototype must be
 *        of the form:
 *        int ompi_rte_db_store(const ompi_process_name_t *proc,
 *                              const char *key, const void *data,
 *                              opal_data_type_t type);
 *        The implementation of this function must store a COPY of the data
 *        provided - the data is NOT guaranteed to be valid after return
 *        from the call.
 *     3. ompi_rte_db_fetch -
 *        NOTE: Fetch accepts an 'ompi_proc_t'.
 *        int ompi_rte_db_fetch(const struct ompi_proc_t *proc,
 *                              const char *key,
 *                              void **data,
 *                              opal_data_type_t type);
 *     4. ompi_rte_db_fetch_pointer -
 *        NOTE: Fetch accepts an 'ompi_proc_t'.
 *        int ompi_rte_db_fetch_pointer(const struct ompi_proc_t *proc,
 *                                      const char *key,
 *                                      void **data,
 *                                      opal_data_type_t type);
 *     5. Pre-defined db keys (with associated values after rte_init)
 *        a. OMPI_DB_HOSTNAME
 *        b. OMPI_DB_LOCALITY
 *
 * (g) Communication support
 *
 */

#ifndef OMPI_MCA_RTE_H
#define OMPI_MCA_RTE_H

#include "ompi_config.h"
#include "ompi/constants.h"

#include <stdint.h>
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

struct opal_proc_t;

#include "opal/mca/threads/threads.h"
#include "opal/util/proc.h"
#include "opal/mca/hwloc/hwloc-internal.h"
#include "opal/mca/pmix/pmix-internal.h"

struct ompi_proc_t;
struct ompi_communicator_t;

BEGIN_C_DECLS

/* Process name objects and operations */
typedef opal_process_name_t ompi_process_name_t;
typedef uint32_t ompi_jobid_t;
typedef uint32_t ompi_vpid_t;

/* some local storage */
OMPI_DECLSPEC extern hwloc_cpuset_t ompi_proc_applied_binding;

#define OMPI_PROC_MY_NAME (&opal_process_info.my_name)
#define OMPI_NAME_WILDCARD  (&opal_name_wildcard)
#define OMPI_PROC_MYID (&opal_process_info.myprocid)
#define OMPI_PRINT_ID(a) ompi_pmix_print_id(a)
OMPI_DECLSPEC char* ompi_pmix_print_id(const pmix_proc_t *procid);

typedef uint8_t ompi_rte_cmp_bitmask_t;
#define OMPI_RTE_CMP_NONE   0x00
#define OMPI_RTE_CMP_JOBID  0x02
#define OMPI_RTE_CMP_VPID   0x04
#define OMPI_RTE_CMP_ALL    0x0f
#define OMPI_RTE_CMP_WILD   0x10

OMPI_DECLSPEC char* ompi_pmix_print_name(const ompi_process_name_t *name);

#define OMPI_NAME_PRINT(a) ompi_pmix_print_name(a)
OMPI_DECLSPEC int ompi_rte_compare_name_fields(ompi_rte_cmp_bitmask_t mask,
                                               const opal_process_name_t* name1,
                                               const opal_process_name_t* name2);
OMPI_DECLSPEC int ompi_rte_convert_string_to_process_name(opal_process_name_t *name,
                                                          const char* name_string);
OMPI_DECLSPEC int ompi_rte_convert_process_name_to_string(char** name_string,
                                                          const opal_process_name_t *name);

OMPI_DECLSPEC void ompi_rte_breakpoint(char *name);

#define OMPI_LOCAL_JOBID(n) \
    ( (n) & 0x0000ffff)
#define OMPI_JOB_FAMILY(n)  \
    (((n) >> 16) & 0x0000ffff)
#define OMPI_CONSTRUCT_LOCAL_JOBID(local, job) \
    ( ((local) & 0xffff0000) | ((job) & 0x0000ffff) )
#define OMPI_CONSTRUCT_JOB_FAMILY(n) \
    ( ((n) << 16) & 0xffff0000)

#define OMPI_CONSTRUCT_JOBID(family, local) \
    OMPI_CONSTRUCT_LOCAL_JOBID(OMPI_CONSTRUCT_JOB_FAMILY(family), local)

/* This is the DSS tag to serialize a proc name */
#define OMPI_NAME OPAL_NAME
#define OMPI_PROCESS_NAME_HTON OPAL_PROCESS_NAME_HTON
#define OMPI_PROCESS_NAME_NTOH OPAL_PROCESS_NAME_NTOH

#if OPAL_ENABLE_DEBUG
static inline opal_process_name_t * OMPI_CAST_RTE_NAME(opal_process_name_t * name) {
    return (opal_process_name_t *)name;
}
#else
#define OMPI_CAST_RTE_NAME(a) ((opal_process_name_t*)(a))
#endif

/* Process info struct and values */
#define ompi_process_info opal_process_info
#define ompi_rte_proc_is_bound opal_process_info.proc_is_bound

/* Error handling objects and operations */
OMPI_DECLSPEC void __opal_attribute_noreturn__
  ompi_rte_abort(int error_code, char *fmt, ...);
OMPI_DECLSPEC void ompi_rte_abort_peers(opal_process_name_t *procs,
                                        int32_t num_procs,
                                        int error_code);
#define OMPI_ERROR_LOG OPAL_ERROR_LOG

/* Init, update and finalize operations */
OMPI_DECLSPEC int ompi_rte_init(int *argc, char ***argv);
OMPI_DECLSPEC int ompi_rte_refresh_job_size(void);
OMPI_DECLSPEC int ompi_rte_refresh_peers(bool self);
OMPI_DECLSPEC int ompi_rte_finalize(void);
OMPI_DECLSPEC void ompi_rte_wait_for_debugger(void);
OMPI_DECLSPEC extern bool ompi_singleton;

/* In a few places, we need to barrier until something happens
 * that changes a flag to indicate we can release - e.g., waiting
 * for a specific RTE message to arrive. We don't want to block MPI
 * progress while waiting, so we loop over opal_progress, letting
 * the RTE progress thread move the RTE along
 */
#define OMPI_WAIT_FOR_COMPLETION(flg)                                       \
    do {                                                                    \
        while ((flg)) {                                                     \
            opal_progress();                                                \
        }                                                                   \
    }while(0);

#define OMPI_LAZY_WAIT_FOR_COMPLETION(flg)                                  \
    do {                                                                    \
        while ((flg)) {                                                     \
            opal_progress();                                                \
            usleep(100);                                                    \
        }                                                                   \
    }while(0);

END_C_DECLS

#endif /* OMPI_RTE_H_ */
