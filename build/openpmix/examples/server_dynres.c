/*
 * Copyright (c) 2004-2010 The Trustees of Indiana University and Indiana
 *                         University Research and Technology
 *                         Corporation.  All rights reserved.
 * Copyright (c) 2004-2011 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2004-2005 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 * Copyright (c) 2004-2005 The Regents of the University of California.
 *                         All rights reserved.
 * Copyright (c) 2006-2013 Los Alamos National Security, LLC.
 *                         All rights reserved.
 * Copyright (c) 2009-2012 Cisco Systems, Inc.  All rights reserved.
 * Copyright (c) 2011      Oak Ridge National Labs.  All rights reserved.
 * Copyright (c) 2013-2020 Intel, Inc.  All rights reserved.
 * Copyright (c) 2015      Research Organization for Information Science
 *                         and Technology (RIST). All rights reserved.
 * Copyright (c) 2016      IBM Corporation.  All rights reserved.
 * Copyright (c) 2021      Nanook Consulting.  All rights reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 *
 */

#include "src/include/pmix_config.h"
#include "../include/pmix_server.h"
#include "src/include/pmix_globals.h"
#include "src/include/types.h"
#include "src/server/pmix_server_ops.h"


#include <dirent.h>
#include <errno.h>
#include <pwd.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <time.h>
#include <unistd.h>
#include <sched.h>

#include "src/class/pmix_list.h"
#include "src/util/argv.h"
#include "src/util/output.h"
#include "src/util/pmix_environ.h"
#include "src/util/printf.h"



static pmix_status_t connected(const pmix_proc_t *proc, void *server_object,
                               pmix_op_cbfunc_t cbfunc, void *cbdata);
static pmix_status_t finalized(const pmix_proc_t *proc, void *server_object,
                               pmix_op_cbfunc_t cbfunc, void *cbdata);
static pmix_status_t abort_fn(const pmix_proc_t *proc, void *server_object, int status,
                              const char msg[], pmix_proc_t procs[], size_t nprocs,
                              pmix_op_cbfunc_t cbfunc, void *cbdata);
static pmix_status_t fencenb_fn(const pmix_proc_t procs[], size_t nprocs, const pmix_info_t info[],
                                size_t ninfo, char *data, size_t ndata, pmix_modex_cbfunc_t cbfunc,
                                void *cbdata);
static pmix_status_t dmodex_fn(const pmix_proc_t *proc, const pmix_info_t info[], size_t ninfo,
                               pmix_modex_cbfunc_t cbfunc, void *cbdata);
static pmix_status_t publish_fn(const pmix_proc_t *proc, const pmix_info_t info[], size_t ninfo,
                                pmix_op_cbfunc_t cbfunc, void *cbdata);
static pmix_status_t lookup_fn(const pmix_proc_t *proc, char **keys, const pmix_info_t info[],
                               size_t ninfo, pmix_lookup_cbfunc_t cbfunc, void *cbdata);
static pmix_status_t unpublish_fn(const pmix_proc_t *proc, char **keys, const pmix_info_t info[],
                                  size_t ninfo, pmix_op_cbfunc_t cbfunc, void *cbdata);
static pmix_status_t spawn_fn(const pmix_proc_t *proc, const pmix_info_t job_info[], size_t ninfo,
                              const pmix_app_t apps[], size_t napps, pmix_spawn_cbfunc_t cbfunc,
                              void *cbdata);
static pmix_status_t connect_fn(const pmix_proc_t procs[], size_t nprocs, const pmix_info_t info[],
                                size_t ninfo, pmix_op_cbfunc_t cbfunc, void *cbdata);
static pmix_status_t disconnect_fn(const pmix_proc_t procs[], size_t nprocs,
                                   const pmix_info_t info[], size_t ninfo, pmix_op_cbfunc_t cbfunc,
                                   void *cbdata);
static pmix_status_t register_event_fn(pmix_status_t *codes, size_t ncodes,
                                       const pmix_info_t info[], size_t ninfo,
                                       pmix_op_cbfunc_t cbfunc, void *cbdata);
static pmix_status_t deregister_events(pmix_status_t *codes, size_t ncodes, pmix_op_cbfunc_t cbfunc,
                                       void *cbdata);
static pmix_status_t notify_event(pmix_status_t code, const pmix_proc_t *source,
                                  pmix_data_range_t range, pmix_info_t info[], size_t ninfo,
                                  pmix_op_cbfunc_t cbfunc, void *cbdata);
static pmix_status_t query_fn(pmix_proc_t *proct, pmix_query_t *queries, size_t nqueries,
                              pmix_info_cbfunc_t cbfunc, void *cbdata);
static void tool_connect_fn(pmix_info_t *info, size_t ninfo, pmix_tool_connection_cbfunc_t cbfunc,
                            void *cbdata);
static void log_fn(const pmix_proc_t *client, const pmix_info_t data[], size_t ndata,
                   const pmix_info_t directives[], size_t ndirs, pmix_op_cbfunc_t cbfunc,
                   void *cbdata);
static pmix_status_t pset_operation_fn( const pmix_proc_t *client,
                                        pmix_psetop_directive_t directive,
                                        const pmix_info_t data[], size_t ndata,
                                        pmix_psetop_cbfunc_t cbfunc, void *cbdata);
static pmix_status_t group_fn(pmix_group_operation_t op, char *gpid, const pmix_proc_t procs[],
                                   size_t nprocs, const pmix_info_t directives[], size_t ndirs,
                                   pmix_info_cbfunc_t cbfunc, void *cbdata);

static pmix_server_module_t mymodule = {.client_connected = connected,
                                        .client_finalized = finalized,
                                        .abort = abort_fn,
                                        .fence_nb = fencenb_fn,
                                        .direct_modex = dmodex_fn,
                                        .publish = publish_fn,
                                        .lookup = lookup_fn,
                                        .unpublish = unpublish_fn,
                                        .spawn = spawn_fn,
                                        .connect = connect_fn,
                                        .disconnect = disconnect_fn,
                                        .register_events = register_event_fn,
                                        .deregister_events = deregister_events,
                                        .notify_event = notify_event,
                                        .query = query_fn,
                                        .tool_connected = tool_connect_fn,
                                        .log = log_fn,
                                        .pset_operation= pset_operation_fn,
                                        .group=group_fn};

typedef struct {
    pmix_list_item_t super;
    pmix_pdata_t pdata;
} pmix_locdat_t;
PMIX_CLASS_INSTANCE(pmix_locdat_t, pmix_list_item_t, NULL, NULL);

#define PMIX_WAIT_FOR_COMPLETION(a) \
    do {                            \
        while ((a)) {               \
            usleep(10);             \
        }                           \
        PMIX_ACQUIRE_OBJECT((a));   \
    } while (0)

typedef struct {
    pmix_object_t super;
    volatile bool active;
    pmix_proc_t caller;
    pmix_info_t *info;
    size_t ninfo;
    pmix_op_cbfunc_t cbfunc;
    pmix_spawn_cbfunc_t spcbfunc;
    void *cbdata;
} myxfer_t;
static void xfcon(myxfer_t *p)
{
    p->info = NULL;
    p->ninfo = 0;
    p->active = true;
    p->cbfunc = NULL;
    p->spcbfunc = NULL;
    p->cbdata = NULL;
}
static void xfdes(myxfer_t *p)
{
    if (NULL != p->info) {
        PMIX_INFO_FREE(p->info, p->ninfo);
    }
}
PMIX_CLASS_INSTANCE(myxfer_t, pmix_object_t, xfcon, xfdes);

typedef struct {
    pmix_list_item_t super;
    pid_t pid;
} wait_tracker_t;
PMIX_CLASS_INSTANCE(wait_tracker_t, pmix_list_item_t, NULL, NULL);

static volatile int wakeup;
static pmix_list_t pubdata;
static pmix_event_t handler;
static pmix_list_t children;
static pmix_list_t procs_universe;
static uint32_t num_psets=0;
static int num_ranks_universe=0;
static int rc_counter=0;
pthread_mutex_t rc_mutex;
char rc_pset[PMIX_MAX_KEYLEN];
char rc_tag[PMIX_MAX_KEYLEN];
char rc_op[PMIX_MAX_KEYLEN];

static volatile bool rc_in_chain=false;
static volatile int cur_rc_type;
static volatile int cur_rc_nprocs;   

pmix_status_t rcexec(int _rc_type, int nprocs);

typedef struct {
    pmix_list_item_t super;
    pmix_proc_t proc;
} pmix_proc_t_item;
PMIX_CLASS_INSTANCE(pmix_proc_t_item, pmix_list_item_t, NULL, NULL);

static cpu_set_t mask;

static pthread_mutex_t pset_buf_mutex;

typedef struct _pset_args{
    size_t nmembers;
    pmix_proc_t *procs;
    char name[32];
}pset_args;

typedef struct _pset_buf{
    volatile pset_args *buffer;
    size_t npsets;
    size_t size;
}pset_buf;

static volatile pset_buf pset_buffer;

static void assign_to_core(int core_id){
    CPU_ZERO(&mask);
    CPU_SET(core_id, &mask);
    sched_setaffinity(0, sizeof(mask), &mask);
}

static void release_pset_args(pset_args *args){
    //free(args->name);
    free(args->procs);
}

static void create_pset_args(pset_args *dest, pmix_proc_t *procs, size_t nprocs, char* name){
    int n;
    
    dest->nmembers=nprocs;
    strcpy(dest->name, name);
    PMIX_PROC_CREATE(dest->procs, nprocs);
    for(n=0; n<nprocs; n++){
        PMIX_PROC_LOAD(&dest->procs[n], procs[n].nspace, procs[n].rank);
    }
    
}

static void produce_pset_def(pmix_proc_t *procs, size_t nprocs, char* name){
    
    pthread_mutex_lock(&pset_buf_mutex);
    if(++(pset_buffer.npsets)>=pset_buffer.size){
        pset_buffer.size= pset_buffer.size==0 ? 2 : pset_buffer.size*2;
        pset_buffer.buffer=realloc(pset_buffer.buffer, pset_buffer.size*sizeof(pset_args));
    }
    create_pset_args(&pset_buffer.buffer[pset_buffer.npsets-1],procs, nprocs, name);
    pthread_mutex_unlock(&pset_buf_mutex);
}

static void consume_pset_def(){

    pthread_mutex_lock(&pset_buf_mutex);
    if(pset_buffer.npsets<=0){
        pthread_mutex_unlock(&pset_buf_mutex);
        return 0;
    }

    PMIx_server_define_process_set(pset_buffer.buffer[0].procs, pset_buffer.buffer[0].nmembers, pset_buffer.buffer[0].name);
    num_psets++;
    release_pset_args(pset_buffer.buffer);
    if(--pset_buffer.npsets>0){
        memmove(&pset_buffer.buffer[0], &(pset_buffer.buffer[1]), pset_buffer.npsets*sizeof(pset_args));
    }
    
    if(pset_buffer.npsets<pset_buffer.size/2){
        pset_buffer.buffer=realloc(pset_buffer.buffer, pset_buffer.size/2);
        pset_buffer.size/=2;
    }
    
    pthread_mutex_unlock(&pset_buf_mutex);
}


static void set_namespace(int nprocs, char *ranks, char *nspace, pmix_op_cbfunc_t cbfunc,
                          myxfer_t *x, bool keep_nlocalprocs);
static void errhandler(size_t evhdlr_registration_id, pmix_status_t status,
                       const pmix_proc_t *source, pmix_info_t info[], size_t ninfo,
                       pmix_info_t results[], size_t nresults,
                       pmix_event_notification_cbfunc_fn_t cbfunc, void *cbdata);
static void rchandler(size_t evhdlr_registration_id, pmix_status_t status,
                       const pmix_proc_t *source, pmix_info_t info[], size_t ninfo,
                       pmix_info_t results[], size_t nresults,
                       pmix_event_notification_cbfunc_fn_t cbfunc, void *cbdata);
static void rc_finalize_handler(size_t evhdlr_registration_id, pmix_status_t status,
                       const pmix_proc_t *source, pmix_info_t info[], size_t ninfo,
                       pmix_info_t results[], size_t nresults,
                       pmix_event_notification_cbfunc_fn_t cbfunc, void *cbdata);
static void wait_signal_callback(int fd, short event, void *arg);
static void errhandler_reg_callbk(pmix_status_t status, size_t errhandler_ref, void *cbdata);
//static void rchandler_reg_callbk(pmix_status_t status, size_t errhandler_ref, void *cbdata);

static void opcbfunc(pmix_status_t status, void *cbdata)
{
    myxfer_t *x = (myxfer_t *) cbdata;

    /* release the caller, if necessary */
    if (NULL != x->cbfunc) {
        x->cbfunc(PMIX_SUCCESS, x->cbdata);
    }
    x->active = false;
}


int main(int argc, char **argv)
{
    char **client_env = NULL;
    char **client_argv = NULL;
    char *tmp, **atmp, *executable = NULL, *tmpdir, *cleanup;
    int rc, nprocs = 1, n, k;
    uid_t myuid;
    gid_t mygid;
    pid_t pid;
    myxfer_t *x;
    pmix_proc_t proc;
    wait_tracker_t *child;
    pmix_proc_t_item *proc_item;
    char *tdir;
    uid_t uid = geteuid();
    pmix_info_t *info;
    struct stat buf;
    pset_buffer.npsets=0;
    pset_buffer.size=0;
    rc_tag[0]=rc_pset[0]=rc_op[0]='\0';

    /* define and pass a personal tmpdir to protect the system */
    if (NULL == (tdir = getenv("TMPDIR"))) {
        if (NULL == (tdir = getenv("TEMP"))) {
            if (NULL == (tdir = getenv("TMP"))) {
                tdir = "/tmp";
            }
        }
    }
    if (0 > asprintf(&tmpdir, "%s/pmix.%lu", tdir, (long unsigned) uid)) {
        fprintf(stderr, "Out of memory\n");
        exit(1);
    }
    /* create the directory */
    if (0 != stat(tmpdir, &buf)) {
        /* try to make directory */
        if (0 != mkdir(tmpdir, S_IRWXU)) {
            fprintf(stderr, "Cannot make tmpdir %s", tmpdir);
            exit(1);
        }
    }
    asprintf(&cleanup, "rm -rf %s", tmpdir);
    PMIX_INFO_CREATE(info, 2);
    PMIX_INFO_LOAD(&info[0], PMIX_SERVER_TMPDIR, tmpdir, PMIX_STRING);
    bool tool_support=true;
    PMIX_INFO_LOAD(&info[1], PMIX_SERVER_TOOL_SUPPORT, &tool_support, PMIX_BOOL);

    /* setup the server library */
    if (PMIX_SUCCESS != (rc = PMIx_server_init(&mymodule, info, 2))) {
        fprintf(stderr, "Init failed with error %d\n", rc);
        return rc;
    }
    PMIX_INFO_FREE(info, 1);

    /* register the errhandler */
    //PMIx_Register_event_handler(NULL, 0, NULL, 0, errhandler, errhandler_reg_callbk, NULL);

    pmix_status_t rc_finalized=PMIX_RC_FINALIZED;
    pmix_status_t rc_define=PMIX_RC_DEFINE;
    /* register the resource change cmd handler */
    PMIx_Register_event_handler(&rc_define, 1, NULL, 0, rchandler, errhandler_reg_callbk, NULL);

    /* register the resource change finalized handler */
    PMIx_Register_event_handler(&rc_finalized, 1, NULL, 0, rc_finalize_handler, errhandler_reg_callbk, NULL);

    /* setup the proc list */
    PMIX_CONSTRUCT(&procs_universe, pmix_list_t);

    /* setup the pub data, in case it is used */
    PMIX_CONSTRUCT(&pubdata, pmix_list_t);

    /* setup to see sigchld on the forked tests */
    PMIX_CONSTRUCT(&children, pmix_list_t);
    pmix_event_assign(&handler, pmix_globals.evbase, SIGCHLD, EV_SIGNAL | EV_PERSIST,
                      wait_signal_callback, &handler);
    pmix_event_add(&handler, NULL);

    /* see if we were passed the number of procs to run or
     * the executable to use */
    for (n = 1; n < (argc - 1); n++) {
        if (0 == strcmp("-n", argv[n]) && NULL != argv[n + 1]) {
            nprocs = strtol(argv[n + 1], NULL, 10);
            ++n; // step over the argument
        } else if (0 == strcmp("-e", argv[n]) && NULL != argv[n + 1]) {
            executable = strdup(argv[n + 1]);
            for (k = n + 2; NULL != argv[k]; k++) {
                pmix_argv_append_nosize(&client_argv, argv[k]);
            }
            n += k;
        }
    }
    if (NULL == executable) {
        //executable = strdup("./simpclient");
        executable = strdup("./dynres_c");
    }

    /* we have a single namespace for all clients */
    atmp = NULL;
    for (n = 0; n < nprocs; n++) {
        asprintf(&tmp, "%d", n);
        pmix_argv_append_nosize(&atmp, tmp);
        free(tmp);
    }
    tmp = pmix_argv_join(atmp, ',');
    pmix_argv_free(atmp);
    /* register the nspace */
    x = PMIX_NEW(myxfer_t);
    set_namespace(nprocs, tmp, "foobar", opcbfunc, x, false);

    /* set common argv and env */
    client_env = pmix_argv_copy(environ);
    pmix_argv_prepend_nosize(&client_argv, executable);

    wakeup = nprocs;
    printf("Server will start the test with %d client processes\n", nprocs);
    myuid = getuid();
    mygid = getgid();

    /* if the nspace registration hasn't completed yet,
     * wait for it here */
    PMIX_WAIT_FOR_COMPLETION(x->active);
    free(tmp);
    PMIX_RELEASE(x);

    /* prep the local node for launch */
    x = PMIX_NEW(myxfer_t);
    if (PMIX_SUCCESS != (rc = PMIx_server_setup_local_support("foobar", NULL, 0, opcbfunc, x))) {
        fprintf(stderr, "Setup local support failed: %d\n", rc);
        PMIx_server_finalize();
        system(cleanup);
        return rc;
    }
    PMIX_WAIT_FOR_COMPLETION(x->active);
    PMIX_RELEASE(x);
    assign_to_core(0);
    char psetname[6];
    char delta_pset_name[6];
    sprintf(delta_pset_name, "test1");
    pmix_proc_t *delta_pset_procs=malloc((nprocs)*sizeof(pmix_proc_t));
    /* fork/exec the test */
    (void) strncpy(proc.nspace, "foobar", PMIX_MAX_NSLEN);
    for (n = 0; n < nprocs; n++) {
        proc.rank = n;
        delta_pset_procs[n]=proc;
        /*
        if(n<1){
            sprintf(psetname, "test%d", n+1);
            PMIx_server_define_process_set(&proc,1,psetname);
            num_psets++;
        }else{
            sprintf(psetname, "test%d", n+2);
            //PMIx_server_define_process_set(&proc,1,psetname);
            delta_pset_procs[n-1]=proc;
        }
        */
        
        if (PMIX_SUCCESS != (rc = PMIx_server_setup_fork(&proc, &client_env))) { // n
            fprintf(stderr, "Server fork setup failed with error %d\n", rc);
            PMIx_server_finalize();
            system(cleanup);
            return rc;
        }
        
        x = PMIX_NEW(myxfer_t);
        if (PMIX_SUCCESS
            != (rc = PMIx_server_register_client(&proc, myuid, mygid, NULL, opcbfunc, x))) {
            fprintf(stderr, "Server fork setup failed with error %d\n", rc);
            PMIx_server_finalize();
            system(cleanup);
            return rc;
        }
        /* don't fork/exec the client until we know it is registered
         * so we avoid a potential race condition in the server */
        PMIX_WAIT_FOR_COMPLETION(x->active);
        PMIX_RELEASE(x);
        pid = fork();
        if (pid < 0) {
            fprintf(stderr, "Fork failed\n");
            PMIx_server_finalize();
            system(cleanup);
            return -1;
        }
        child = PMIX_NEW(wait_tracker_t);
        child->pid = pid;
        pmix_list_append(&children, &child->super);

        if (pid == 0) {
            assign_to_core(n%7+1);
            execve(executable, client_argv, client_env);
            /* Does not return */
            exit(0);
        }
        num_ranks_universe++;
        proc_item=PMIX_NEW(pmix_proc_t_item);
        PMIX_PROC_CONSTRUCT(&proc_item->proc);
        PMIX_PROC_LOAD(&proc_item->proc, proc.nspace, proc.rank);
        pmix_list_append(&procs_universe, &proc_item->super);

    }
    
    PMIx_server_define_process_set(delta_pset_procs,nprocs,delta_pset_name);
    num_psets++;

    free(executable);
    pmix_argv_free(client_argv);
    pmix_argv_free(client_env);

    /* hang around until the client(s) finalize */
    while (0 < wakeup) {
        
        struct timespec ts;
        ts.tv_sec = 0;
        ts.tv_nsec = 100000;
        nanosleep(&ts, NULL);
        int rc_nprocs;
        int rc_type=-1;
        pthread_mutex_lock(&rc_mutex);
        if(rc_in_chain){
            rc_type=cur_rc_type;
            rc_nprocs=cur_rc_nprocs;
            rc_in_chain=false;
        }
        pthread_mutex_unlock(&rc_mutex);
        if(rc_type>-1){
            rcexec(rc_type,rc_nprocs);
        }
        consume_pset_def();
    }

    /* deregister the errhandler */
    PMIx_Deregister_event_handler(0, NULL, NULL);

    /* release any pub data */
    PMIX_LIST_DESTRUCT(&pubdata);

    PMIX_LIST_DESTRUCT(&procs_universe);

    /* finalize the server library */
    if (PMIX_SUCCESS != (rc = PMIx_server_finalize())) {
        fprintf(stderr, "Finalize failed with error %d\n", rc);
    }
    free(pset_buffer.buffer);
    free(delta_pset_procs);
    fprintf(stderr, "Test finished OK!\n");
    system(cleanup);

    return rc;
}

static void setup_cbfunc(pmix_status_t status, pmix_info_t info[], size_t ninfo,
                         void *provided_cbdata, pmix_op_cbfunc_t cbfunc, void *cbdata)
{
    myxfer_t *myxfer = (myxfer_t *) provided_cbdata;
    size_t i;

    if (PMIX_SUCCESS == status && 0 < ninfo) {
        myxfer->ninfo = ninfo;
        PMIX_INFO_CREATE(myxfer->info, ninfo);
        for (i = 0; i < ninfo; i++) {
            PMIX_INFO_XFER(&myxfer->info[i], &info[i]);
        }
    }
    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, cbdata);
    }
    myxfer->active = false;
}

static void pmix_relfn(void *cbdata){
    pmix_info_caddy_t *cd = (pmix_info_caddy_t *)cbdata;
    if(NULL != cd->info){
        PMIX_INFO_FREE(cd->info, cd->ninfo);
    }
    PMIX_RELEASE(cd);
}

static void set_namespace(int nprocs, char *ranks, char *nspace, pmix_op_cbfunc_t cbfunc,
                          myxfer_t *x, bool keep_nlocalprocs)
{
    char *regex, *ppn;
    char hostname[PMIX_MAXHOSTNAMELEN];
    pmix_status_t rc;
    myxfer_t myxfer;
    size_t i = 0;

    gethostname(hostname, sizeof(hostname));

    /* request application setup information - e.g., network
     * security keys or endpoint info */
    PMIX_CONSTRUCT(&myxfer, myxfer_t);
    myxfer.active = true;
    if (PMIX_SUCCESS
        != (rc = PMIx_server_setup_application(nspace, NULL, 0, setup_cbfunc, &myxfer))) {
        PMIX_DESTRUCT(&myxfer);
        fprintf(stderr, "Failed to setup application: %d\n", rc);
        exit(1);
    }
    PMIX_WAIT_FOR_COMPLETION(myxfer.active);
    x->ninfo = myxfer.ninfo + 7;

    PMIX_INFO_CREATE(x->info, x->ninfo);
    if (0 < myxfer.ninfo) {
        for (i = 0; i < myxfer.ninfo; i++) {
            PMIX_INFO_XFER(&x->info[i], &myxfer.info[i]);
        }
    }
    PMIX_DESTRUCT(&myxfer);

    (void) strncpy(x->info[i].key, PMIX_UNIV_SIZE, PMIX_MAX_KEYLEN);
    x->info[i].value.type = PMIX_UINT32;
    x->info[i].value.data.uint32 = nprocs;

    ++i;
    (void) strncpy(x->info[i].key, PMIX_SPAWNED, PMIX_MAX_KEYLEN);
    x->info[i].value.type = PMIX_UINT32;
    x->info[i].value.data.uint32 = 0;

    ++i;
    (void) strncpy(x->info[i].key, PMIX_LOCAL_SIZE, PMIX_MAX_KEYLEN);
    x->info[i].value.type = PMIX_UINT32;
    x->info[i].value.data.uint32 = nprocs;

    ++i;
    (void) strncpy(x->info[i].key, PMIX_LOCAL_PEERS, PMIX_MAX_KEYLEN);
    x->info[i].value.type = PMIX_STRING;
    x->info[i].value.data.string = strdup(ranks);

    ++i;
    PMIx_generate_regex(hostname, &regex);
    (void) strncpy(x->info[i].key, PMIX_NODE_MAP, PMIX_MAX_KEYLEN);
    x->info[i].value.type = PMIX_STRING;
    x->info[i].value.data.string = regex;

    ++i;
    PMIx_generate_ppn(ranks, &ppn);
    (void) strncpy(x->info[i].key, PMIX_PROC_MAP, PMIX_MAX_KEYLEN);
    x->info[i].value.type = PMIX_STRING;
    x->info[i].value.data.string = ppn;

    ++i;
    (void) strncpy(x->info[i].key, PMIX_JOB_SIZE, PMIX_MAX_KEYLEN);
    x->info[i].value.type = PMIX_UINT32;
    x->info[i].value.data.uint32 = nprocs;

    if(keep_nlocalprocs){
        nprocs=INT_MIN;
    }

    PMIx_server_register_nspace(nspace, nprocs, x->info, x->ninfo, cbfunc, x);
}

static void errhandler(size_t evhdlr_registration_id, pmix_status_t status,
                       const pmix_proc_t *source, pmix_info_t info[], size_t ninfo,
                       pmix_info_t results[], size_t nresults,
                       pmix_event_notification_cbfunc_fn_t cbfunc, void *cbdata)
{
    pmix_output(0, "SERVER: ERRHANDLER CALLED WITH STATUS %d", status);
}

static void errhandler_reg_callbk(pmix_status_t status, size_t errhandler_ref, void *cbdata)
{
    return;
}

pmix_status_t rcexec(int _rc_type, int nprocs){
    // need a lock
    char **client_env = NULL;
    char **client_argv = NULL;
    char *tmp, **atmp, *executable = NULL, *tmpdir, *cleanup;
    int rc= 1, n, k;
    uid_t myuid;
    gid_t mygid;
    pid_t pid;
    myxfer_t *x;
    pmix_proc_t proc;
    wait_tracker_t *child;
    pmix_proc_t_item *proc_item;
    char *tdir;
    uid_t uid = geteuid();
    pmix_info_t *info;
    struct stat buf;

    if (NULL == executable) {
        //executable = strdup("./simpclient");
        executable = strdup("./dynres_c");
    }
    /* TODO: choose new nspace */
    char nspace[PMIX_MAX_NSLEN];
    //sprintf(nspace, "foobar%d", rc_counter);
    sprintf(nspace, "foobar");
    /* set common argv and env */
    client_env = pmix_argv_copy(environ);
    pmix_argv_prepend_nosize(&client_argv, executable);
    /* we have a single namespace for all clients */
    atmp = NULL;



    pthread_mutex_lock(&rc_mutex);

    if(_rc_type==0)strcpy(rc_op,"ADD");
    if(_rc_type==1)strcpy(rc_op,"SUB");
    sprintf(rc_pset, "rc%d", rc_counter);
    strcpy(rc_tag,rc_pset);
    rc_counter++;
    pmix_proc_t *delta_pset_procs=malloc((nprocs)*sizeof(pmix_proc_t));

    /* fork/exec the test */
    
    if(_rc_type==0){


        /* if new nspace 
        * for (n = 0; n < nprocs; n++) {
        *     asprintf(&tmp, "%d", n);
        *     pmix_argv_append_nosize(&atmp, tmp);
        *     free(tmp);
        * }
        */
        pmix_proc_t_item *proc_item;
        PMIX_LIST_FOREACH(proc_item, &procs_universe, pmix_proc_t_item){
            asprintf(&tmp, "%d", proc_item->proc.rank);
            pmix_argv_append_nosize(&atmp, tmp);
            free(tmp);
        }
        int start_rank=num_ranks_universe;
        int end_rank=start_rank+nprocs;
        for (n = start_rank; n < end_rank; n++) {
            asprintf(&tmp, "%d", n);
            pmix_argv_append_nosize(&atmp, tmp);
            free(tmp);
        }

        tmp = pmix_argv_join(atmp, ',');
        pmix_argv_free(atmp);
        /* register the nspace */
        x = PMIX_NEW(myxfer_t);

        set_namespace(nprocs+pmix_list_get_size(&procs_universe), tmp, nspace, opcbfunc, x, false);
        //x = PMIX_NEW(myxfer_t);
        //if (PMIX_SUCCESS != (rc = PMIx_server_setup_local_support("foobar", NULL, 0, opcbfunc, x))) {
        //    fprintf(stderr, "Setup local support failed: %d\n", rc);
        //    PMIx_server_finalize();
        //    system(cleanup);
        //    return rc;
        //}
        //PMIX_WAIT_FOR_COMPLETION(x->active);
        //PMIX_RELEASE(x);

        /* set common argv and env */
        client_env = pmix_argv_copy(environ);
        pmix_argv_prepend_nosize(&client_argv, executable);

        /* if the nspace registration hasn't completed yet,
         * wait for it here */
        PMIX_WAIT_FOR_COMPLETION(x->active);
        free(tmp);
        PMIX_RELEASE(x);

        wakeup +=nprocs;
        myuid = getuid();
        mygid = getgid();

        (void) strncpy(proc.nspace, nspace, PMIX_MAX_NSLEN);

        for (n = start_rank; n < end_rank; n++) {
            proc.rank = n;
            delta_pset_procs[n-start_rank]=proc;

            if (PMIX_SUCCESS != (rc = PMIx_server_setup_fork(&proc, &client_env))) { // n
                fprintf(stderr, "Server fork setup failed with error %d\n", rc);
                PMIx_server_finalize();
                system(cleanup);
                return rc;
            }
            x = PMIX_NEW(myxfer_t);
            if (PMIX_SUCCESS
                != (rc = PMIx_server_register_client(&proc, myuid, mygid, NULL, opcbfunc, x))) {
                fprintf(stderr, "Server fork setup failed with error %d\n", rc);
                PMIx_server_finalize();
                system(cleanup);
                return rc;
            }
            /* don't fork/exec the client until we know it is registered
             * so we avoid a potential race condition in the server */
            PMIX_WAIT_FOR_COMPLETION(x->active);
            PMIX_RELEASE(x);

            pid = fork();
            if (pid < 0) {
                fprintf(stderr, "Fork failed\n");
                PMIx_server_finalize();
                system(cleanup);
                return -1;
            }
            child = PMIX_NEW(wait_tracker_t);
            child->pid = pid;
            pmix_list_append(&children, &child->super);

            if (pid == 0) {
                assign_to_core(n%7+1);
                execve(executable, client_argv, client_env);
                /* Does not return */
                exit(0);
            }
            num_ranks_universe++;
            proc_item=PMIX_NEW(pmix_proc_t_item);
            PMIX_PROC_CONSTRUCT(&proc_item->proc);
            PMIX_PROC_LOAD(&proc_item->proc, proc.nspace, proc.rank);
            pmix_list_append(&procs_universe, &proc_item->super);
        }
        //assign_to_core(n);
    }else if(_rc_type==1){
        n=0;
        pmix_proc_t_item *proc_to_remove, *prev; 
        PMIX_LIST_FOREACH_SAFE_REV(proc_to_remove, prev, &procs_universe, pmix_proc_t_item){
            proc.rank = proc_to_remove->proc.rank;
            (void)strncpy(proc.nspace, proc_to_remove->proc.nspace, PMIX_MAX_NSLEN);
            delta_pset_procs[n]=proc;
            pmix_list_remove_item (&procs_universe, (pmix_list_item_t *) proc_to_remove);
            if(++n>=nprocs)break;
        }
        pmix_proc_t_item *proc_item;
        PMIX_LIST_FOREACH(proc_item, &procs_universe, pmix_proc_t_item){
            asprintf(&tmp, "%d", proc_item->proc.rank);
            pmix_argv_append_nosize(&atmp, tmp);
            free(tmp);
        }
        tmp = pmix_argv_join(atmp, ',');
        pmix_argv_free(atmp);
        /* register the nspace */
        x = PMIX_NEW(myxfer_t);

        set_namespace(pmix_list_get_size(&procs_universe), tmp, nspace, opcbfunc, x, true);
        free(tmp);
    }
    rc=PMIx_server_define_process_set(delta_pset_procs,nprocs,rc_pset);
    free(delta_pset_procs);
    pmix_argv:free(client_env);
    pthread_mutex_unlock(&rc_mutex);
    return PMIX_SUCCESS;
}

static pmix_status_t parse_rc_cmd(char *cmd, int *_rc_type, int *nprocs){

    char * token= strtok(cmd, " ");
    if(token==NULL || 0!=strncmp(token, "pmix_session", 12)){
        return PMIX_ERR_BAD_PARAM;
    }

    token= strtok(NULL, " ");
    if(token== NULL)return PMIX_ERR_BAD_PARAM;
    if(0 == strncmp(token, "add", 3))*_rc_type=0;
    else if(0 == strncmp(token, "sub", 3))*_rc_type=1;
    else return PMIX_ERR_BAD_PARAM;

    token= strtok(NULL, " ");
    return ((*nprocs=atoi(token))<=0) ? PMIX_ERR_BAD_PARAM : PMIX_SUCCESS;



}

static void rchandler(size_t evhdlr_registration_id, pmix_status_t status,
                       const pmix_proc_t *source, pmix_info_t info[], size_t ninfo,
                       pmix_info_t results[], size_t nresults,
                       pmix_event_notification_cbfunc_fn_t cbfunc, void *cbdata)
{
    
    size_t n;
    pmix_status_t rc=PMIX_SUCCESS;
    size_t sz;
    char *cmd;
    for(n=0; n<ninfo; n++){
        if(0 == strcmp(info[n].key, "PMIX_RC_CMD")){
            PMIX_VALUE_UNLOAD(rc, &info[n].value, (void**)&cmd, &sz);
            pmix_output(0, "SERVER: RECEIVED RC_CMD %s\n", cmd);
            int rc_type, nprocs;
            rc=parse_rc_cmd(cmd, &rc_type, &nprocs);
            if(rc!=PMIX_SUCCESS){
                printf("Error parsing rc command\n");
                return rc;
            }
            if(rc_type==1 && nprocs>=pmix_list_get_size(&procs_universe)){
                return PMIX_ERR_BAD_PARAM;
            }
            /* push rc into list */
            pthread_mutex_lock(&rc_mutex);
            cur_rc_type=rc_type;
            cur_rc_nprocs=nprocs;
            rc_in_chain=true;
            pthread_mutex_unlock(&rc_mutex);
            //rcexec(rc_type, nprocs);
            
        }
    }
    cbfunc(PMIX_SUCCESS, NULL, 0, NULL, NULL, cbdata);
    return rc;
}

static void rc_finalize_handler(size_t evhdlr_registration_id, pmix_status_t status,
                       const pmix_proc_t *source, pmix_info_t info[], size_t ninfo,
                       pmix_info_t results[], size_t nresults,
                       pmix_event_notification_cbfunc_fn_t cbfunc, void *cbdata)
{
    pmix_output(0, "SERVER: RC_FINALIZE_HANDLER CALLED WITH STATUS %d", status);
    size_t n;
    /*
    pmix_status_t rc;
    pmix_data_type_t sz;
    char cmd[PMIX_MAX_KEYLEN];
    for(n=0; n<ninfo; n++){
        if(0 == strcmp(info[n].key, "PMIX_RC_CMD")){
            rc=PMIX_VALUE_UNLOAD(rc, &info[n].value, (void**)&cmd, &sz);
            int rc_type, nprocs;
            rc=parse_rc_cmd(cmd, &rc_type, &nprocs);
            rcexec(rc_type, nprocs);
            return rc;
        }
    }
    */
    pthread_mutex_lock(&rc_mutex);
    /*
    if(0==strcmp(rc_op, "SUB")){
    
        pmix_pset_t *pset_list_iter;
        pmix_pset_t *pset1=NULL;
        PMIX_LIST_FOREACH(pset_list_iter, &pmix_server_globals.psets, pmix_pset_t){
            if(0 == strcmp(pset_list_iter->name, rc_pset)){
                pset1=pset_list_iter;
                break;
            }
        }
        for(n=0; n<pset1->nmembers; n++){
            (void*)pmix_list_remove_last(&procs_universe);
        }

    }*/
    rc_tag[0]=rc_pset[0]=rc_op[0]='\0';
    pthread_mutex_unlock(&rc_mutex);
    if(NULL != cbfunc){
        cbfunc(PMIX_SUCCESS, NULL, 0, NULL, NULL, cbdata);
    }
    
    return PMIX_SUCCESS;

}


static pmix_status_t connected(const pmix_proc_t *proc, void *server_object,
                               pmix_op_cbfunc_t cbfunc, void *cbdata)
{
    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, cbdata);
    }
    return PMIX_SUCCESS;
}
static pmix_status_t finalized(const pmix_proc_t *proc, void *server_object,
                               pmix_op_cbfunc_t cbfunc, void *cbdata)
{
    pmix_output(0, "SERVER: FINALIZED %s:%d", proc->nspace, proc->rank);
    --wakeup;
    /* ensure we call the cbfunc so the proc can exit! */
    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, cbdata);
    }
    return PMIX_SUCCESS;
}

static void abcbfunc(pmix_status_t status, void *cbdata)
{
    myxfer_t *x = (myxfer_t *) cbdata;

    /* be sure to release the caller */
    if (NULL != x->cbfunc) {
        x->cbfunc(status, x->cbdata);
    }
    PMIX_RELEASE(x);
}

static pmix_status_t abort_fn(const pmix_proc_t *proc, void *server_object, int status,
                              const char msg[], pmix_proc_t procs[], size_t nprocs,
                              pmix_op_cbfunc_t cbfunc, void *cbdata)
{
    pmix_status_t rc;
    myxfer_t *x;

    if (NULL != procs) {
        pmix_output(0, "SERVER: ABORT on %s:%d", procs[0].nspace, procs[0].rank);
    } else {
        pmix_output(0, "SERVER: ABORT OF ALL PROCS IN NSPACE %s", proc->nspace);
    }

    /* instead of aborting the specified procs, notify them
     * (if they have registered their errhandler) */

    /* use the myxfer_t object to ensure we release
     * the caller when notification has been queued */
    x = PMIX_NEW(myxfer_t);
    (void) strncpy(x->caller.nspace, proc->nspace, PMIX_MAX_NSLEN);
    x->caller.rank = proc->rank;

    PMIX_INFO_CREATE(x->info, 2);
    (void) strncpy(x->info[0].key, "DARTH", PMIX_MAX_KEYLEN);
    x->info[0].value.type = PMIX_INT8;
    x->info[0].value.data.int8 = 12;
    (void) strncpy(x->info[1].key, "VADER", PMIX_MAX_KEYLEN);
    x->info[1].value.type = PMIX_DOUBLE;
    x->info[1].value.data.dval = 12.34;
    x->cbfunc = cbfunc;
    x->cbdata = cbdata;

    if (PMIX_SUCCESS
        != (rc = PMIx_Notify_event(status, &x->caller, PMIX_RANGE_NAMESPACE, x->info, 2, abcbfunc,
                                   x))) {
        pmix_output(0, "SERVER: FAILED NOTIFY ERROR %d", (int) rc);
    }

    return PMIX_SUCCESS;
}

static pmix_status_t fencenb_fn(const pmix_proc_t procs[], size_t nprocs, const pmix_info_t info[],
                                size_t ninfo, char *data, size_t ndata, pmix_modex_cbfunc_t cbfunc,
                                void *cbdata)
{
    pmix_output(0, "SERVER: FENCENB");
    /* pass the provided data back to each participating proc */
    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, data, ndata, cbdata, NULL, NULL);
    }
    return PMIX_SUCCESS;
}

static pmix_status_t dmodex_fn(const pmix_proc_t *proc, const pmix_info_t info[], size_t ninfo,
                               pmix_modex_cbfunc_t cbfunc, void *cbdata)
{
    pmix_output(0, "SERVER: DMODEX");

    /* we don't have any data for remote procs as this
     * test only runs one server - so report accordingly */
    if (NULL != cbfunc) {
        cbfunc(PMIX_ERR_NOT_FOUND, NULL, 0, cbdata, NULL, NULL);
    }
    return PMIX_SUCCESS;
}

static pmix_status_t publish_fn(const pmix_proc_t *proc, const pmix_info_t info[], size_t ninfo,
                                pmix_op_cbfunc_t cbfunc, void *cbdata)
{
    pmix_locdat_t *p;
    size_t n;

    pmix_output(2, "SERVER: PUBLISH");

    for (n = 0; n < ninfo; n++) {
        p = PMIX_NEW(pmix_locdat_t);
        (void) strncpy(p->pdata.proc.nspace, proc->nspace, PMIX_MAX_NSLEN);
        p->pdata.proc.rank = proc->rank;
        (void) strncpy(p->pdata.key, info[n].key, PMIX_MAX_KEYLEN);
        pmix_value_xfer(&p->pdata.value, (pmix_value_t *) &info[n].value);
        pmix_list_append(&pubdata, &p->super);
    }
    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, cbdata);
    }
    return PMIX_SUCCESS;
}

static pmix_status_t lookup_fn(const pmix_proc_t *proc, char **keys, const pmix_info_t info[],
                               size_t ninfo, pmix_lookup_cbfunc_t cbfunc, void *cbdata)
{
    pmix_locdat_t *p, *p2;
    pmix_list_t results;
    size_t i, n;
    pmix_pdata_t *pd = NULL;
    pmix_status_t ret = PMIX_ERR_NOT_FOUND;

    pmix_output(2, "SERVER: LOOKUP");

    PMIX_CONSTRUCT(&results, pmix_list_t);

    for (n = 0; NULL != keys[n]; n++) {
        PMIX_LIST_FOREACH (p, &pubdata, pmix_locdat_t) {
            if (0 == strncmp(keys[n], p->pdata.key, PMIX_MAX_KEYLEN)) {
                p2 = PMIX_NEW(pmix_locdat_t);
                (void) strncpy(p2->pdata.proc.nspace, p->pdata.proc.nspace, PMIX_MAX_NSLEN);
                p2->pdata.proc.rank = p->pdata.proc.rank;
                (void) strncpy(p2->pdata.key, p->pdata.key, PMIX_MAX_KEYLEN);
                pmix_value_xfer(&p2->pdata.value, &p->pdata.value);
                pmix_list_append(&results, &p2->super);
                break;
            }
        }
    }
    if (0 < (n = pmix_list_get_size(&results))) {
        ret = PMIX_SUCCESS;
        PMIX_PDATA_CREATE(pd, n);
        for (i = 0; i < n; i++) {
            p = (pmix_locdat_t *) pmix_list_remove_first(&results);
            if (p) {
                (void) strncpy(pd[i].proc.nspace, p->pdata.proc.nspace, PMIX_MAX_NSLEN);
                pd[i].proc.rank = p->pdata.proc.rank;
                (void) strncpy(pd[i].key, p->pdata.key, PMIX_MAX_KEYLEN);
                pmix_value_xfer(&pd[i].value, &p->pdata.value);
            }
        }
    }
    PMIX_LIST_DESTRUCT(&results);
    if (NULL != cbfunc) {
        cbfunc(ret, pd, n, cbdata);
    }
    if (0 < n) {
        PMIX_PDATA_FREE(pd, n);
    }
    return PMIX_SUCCESS;
}

static pmix_status_t unpublish_fn(const pmix_proc_t *proc, char **keys, const pmix_info_t info[],
                                  size_t ninfo, pmix_op_cbfunc_t cbfunc, void *cbdata)
{
    pmix_locdat_t *p, *p2;
    size_t n;

    pmix_output(0, "SERVER: UNPUBLISH");

    for (n = 0; NULL != keys[n]; n++) {
        PMIX_LIST_FOREACH_SAFE (p, p2, &pubdata, pmix_locdat_t) {
            if (0 == strncmp(keys[n], p->pdata.key, PMIX_MAX_KEYLEN)) {
                pmix_list_remove_item(&pubdata, &p->super);
                PMIX_RELEASE(p);
                break;
            }
        }
    }
    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, cbdata);
    }
    return PMIX_SUCCESS;
}

static void spcbfunc(pmix_status_t status, void *cbdata)
{
    myxfer_t *x = (myxfer_t *) cbdata;

    if (NULL != x->spcbfunc) {
        x->spcbfunc(PMIX_SUCCESS, "DYNSPACE", x->cbdata);
    }
}

static pmix_status_t spawn_fn(const pmix_proc_t *proc, const pmix_info_t job_info[], size_t ninfo,
                              const pmix_app_t apps[], size_t napps, pmix_spawn_cbfunc_t cbfunc,
                              void *cbdata)
{
    myxfer_t *x;

    pmix_output(0, "SERVER: SPAWN");

    /* in practice, we would pass this request to the local
     * resource manager for launch, and then have that server
     * execute our callback function. For now, we will fake
     * the spawn and just pretend */

    /* must register the nspace for the new procs before
     * we return to the caller */
    x = PMIX_NEW(myxfer_t);
    x->spcbfunc = cbfunc;
    x->cbdata = cbdata;

    set_namespace(2, "0,1", "DYNSPACE", spcbfunc, x, false);

    return PMIX_SUCCESS;
}

static pmix_status_t connect_fn(const pmix_proc_t procs[], size_t nprocs, const pmix_info_t info[],
                                size_t ninfo, pmix_op_cbfunc_t cbfunc, void *cbdata)
{
    pmix_output(0, "SERVER: CONNECT");

    /* in practice, we would pass this request to the local
     * resource manager for handling */

    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, cbdata);
    }

    return PMIX_SUCCESS;
}

static pmix_status_t disconnect_fn(const pmix_proc_t procs[], size_t nprocs,
                                   const pmix_info_t info[], size_t ninfo, pmix_op_cbfunc_t cbfunc,
                                   void *cbdata)
{
    pmix_output(0, "SERVER: DISCONNECT");

    /* in practice, we would pass this request to the local
     * resource manager for handling */

    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, cbdata);
    }

    return PMIX_SUCCESS;
}

static pmix_status_t register_event_fn(pmix_status_t *codes, size_t ncodes,
                                       const pmix_info_t info[], size_t ninfo,
                                       pmix_op_cbfunc_t cbfunc, void *cbdata)
{
    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, cbdata);
    }
    return PMIX_SUCCESS;
}

static pmix_status_t deregister_events(pmix_status_t *codes, size_t ncodes, pmix_op_cbfunc_t cbfunc,
                                       void *cbdata)
{
    return PMIX_SUCCESS;
}

static pmix_status_t notify_event(pmix_status_t code, const pmix_proc_t *source,
                                  pmix_data_range_t range, pmix_info_t info[], size_t ninfo,
                                  pmix_op_cbfunc_t cbfunc, void *cbdata)
{
    return PMIX_SUCCESS;
}

typedef struct query_data_t {
    pmix_info_t *data;
    size_t ndata;
} query_data_t;

static int proc_cmp(pmix_proc_t p1, pmix_proc_t p2){
    return (0==strcmp(p1.nspace, p2.nspace) && p1.rank==p2.rank);
}

static pmix_status_t pset_intersection(pmix_pset_t p1, pmix_pset_t p2, pmix_proc_t *result, size_t *nmembers){
    size_t n, k;
    size_t res_ptr=0;
    size_t nprocs_max = p1.nmembers +p2.nmembers;

    /* Greedily fill in all procs from p2 which are not in p1 */
    for(n=0; n<p2.nmembers; n++){
        int found=0;
        for(k=0; k<p1.nmembers; k++){
            found+= proc_cmp(p2.members[n], p1.members[k]);
            if(found){
                PMIX_PROC_LOAD(&result[res_ptr],p2.members[n].nspace, p2.members[n].rank);
                res_ptr++;
                break;
            }
        }
    }
    *nmembers=res_ptr;
    return PMIX_SUCCESS;

}

static pmix_status_t pset_difference(pmix_pset_t p1, pmix_pset_t p2, pmix_proc_t *result, size_t *nmembers){
    size_t n, k;
    size_t res_ptr=0;
    size_t nprocs_max = p1.nmembers +p2.nmembers;

    /* Greedily fill in all procs from p2 which are not in p1 */
    for(n=0; n<p1.nmembers; n++){
        int found=0;
        for(k=0; k<p2.nmembers; k++){

            found+= proc_cmp(p1.members[n], p2.members[k]);
        }
        if(!found){
            PMIX_PROC_LOAD(&result[res_ptr],p1.members[n].nspace, p1.members[n].rank);
            res_ptr++;
        }
    }
    *nmembers=res_ptr;
    return PMIX_SUCCESS;

}

static pmix_status_t pset_union(pmix_pset_t p1, pmix_pset_t p2, pmix_proc_t *result, size_t *nmembers){
    size_t n, k;
    size_t res_ptr=0;
    size_t nprocs_max = p1.nmembers +p2.nmembers;

    /* fill in all procs from p1 */
    for(n=0; n<p1.nmembers; n++){
        PMIX_PROC_LOAD(&result[res_ptr],p1.members[n].nspace, p1.members[n].rank);
        res_ptr++;
    }

    /* Greedily fill in all procs from p2 which are not in p1 (b.c. procs from p1 were already added) */
    for(n=0; n<p2.nmembers; n++){
        int found=0;
        for(k=0; k<p1.nmembers; k++){

            found+= proc_cmp(p2.members[n], p1.members[k]);
        }
        if(!found){
            PMIX_PROC_LOAD(&result[res_ptr],p2.members[n].nspace, p2.members[n].rank);
            res_ptr++;
        }
    }
    *nmembers=res_ptr;

    return PMIX_SUCCESS;

}




static pmix_status_t pset_operation_fn( const pmix_proc_t *client,
                                        pmix_psetop_directive_t directive,
                                        const pmix_info_t data[], size_t ndata,
                                        pmix_psetop_cbfunc_t cbfunc, void *cbdata)
{
    size_t n, k, sz;
    size_t nresp=0;
    pmix_status_t rc;
    pmix_info_t *info;
    char *pset1_name, *pset2_name, *pset_result_name;



    if (NULL == cbfunc) {
        return PMIX_ERROR;
    }

    for (n = 0; n < ndata; n++) { 
        if(0 == strcmp(data[n].key, PMIX_PSETOP_P1)){
            rc=PMIX_VALUE_UNLOAD(rc, &data[n].value, (void**)&pset1_name, &sz);
        }else if(0 == strcmp(data[n].key, PMIX_PSETOP_P2)){
            rc=PMIX_VALUE_UNLOAD(rc, &data[n].value, (void**)&pset2_name, &sz);
        }else if(0 == strcmp(data[n].key, PMIX_PSETOP_PREF_NAME)){
            rc=PMIX_VALUE_UNLOAD(rc, &data[n].value, (void**)&pset_result_name, &sz);
        }
    }
    if(pset1_name==NULL || pset2_name==NULL){
        printf("pstop failed: no such psets\n");
        free(pset1_name);
        free(pset2_name);
        free(pset_result_name);
        return PMIX_ERR_BAD_PARAM;
    }
    /* TODO: check if preferred name is okay */
    char *pset_result_name_used=pset_result_name;
    

    /* pass name of new set to callback */
    PMIX_INFO_CREATE(info, 3);
    PMIX_INFO_LOAD(&info[0], PMIX_PSETOP_PRESULT, pset_result_name_used, PMIX_STRING);;
    PMIX_INFO_LOAD(&info[1], PMIX_PSETOP_P1, pset1_name, PMIX_STRING);
    PMIX_INFO_LOAD(&info[2], PMIX_PSETOP_P2, pset2_name, PMIX_STRING);

    pmix_info_caddy_t *cd = PMIX_NEW(pmix_info_caddy_t);
    cd->info = info;
    cd->ninfo = 3;
    cbfunc(PMIX_SUCCESS, directive, info, 3, cbdata, pmix_relfn, (void*)cd);
    ///* push pset define command into buffer */
    //produce_pset_def(result_pset_members, nmembers, pset_result_name_used);
    
    free(pset_result_name_used);
    free(pset1_name);
    free(pset2_name);
    return PMIX_SUCCESS;
}

static pmix_status_t group_fn(pmix_group_operation_t op, char *gpid, const pmix_proc_t procs[],
                                   size_t nprocs, const pmix_info_t directives[], size_t ndirs,
                                   pmix_info_cbfunc_t cbfunc, void *cbdata)
{
    
    //prte_pmix_mdx_caddy_t *cd;
    int rc;
    size_t i, mode = 0;
    bool fence = false;
    pmix_byte_object_t *bo = NULL;
    pmix_info_t *results;
    PMIX_INFO_CREATE(results, 1);


    /* they are required to pass us an id */
    if (NULL == gpid) {
        return PMIX_ERR_BAD_PARAM;
    }

    /* check the directives */
    for (i = 0; i < ndirs; i++) {
        /* see if they want a context id assigned */
        if (PMIX_CHECK_KEY(&directives[i], PMIX_GROUP_ASSIGN_CONTEXT_ID)) {
            if (PMIX_INFO_TRUE(&directives[i])) {
                mode = 1;
            }
        } else if (PMIX_CHECK_KEY(&directives[i], PMIX_EMBED_BARRIER)) {
            fence = PMIX_INFO_TRUE(&directives[i]);
        } else if (PMIX_CHECK_KEY(&directives[i], PMIX_GROUP_ENDPT_DATA)) {
            bo = (pmix_byte_object_t *) &directives[i].value.data.bo;
        }
    }

    /* internally track the groups (PRRTE uses psets?)
    if (PMIX_GROUP_CONSTRUCT == op) {
        pset = PRTE_NEW(pmix_server_pset_t);
        pset->name = strdup(gpid);
        pset->num_members = nprocs;
        PMIX_PROC_CREATE(pset->members, pset->num_members);
        memcpy(pset->members, procs, nprocs * sizeof(pmix_proc_t));
        prte_list_append(&prte_pmix_server_globals.psets, &pset->super);
    } else if (PMIX_GROUP_DESTRUCT == op) {
        PRTE_LIST_FOREACH(pset, &prte_pmix_server_globals.psets, pmix_server_pset_t)
        {
            if (0 == strcmp(pset->name, gpid)) {
                prte_list_remove_item(&prte_pmix_server_globals.psets, &pset->super);
                PRTE_RELEASE(pset);
                break;
            }
        }
    }
    */

    /* if they don't want us to do a fence and they don't want a
     * context id assigned, then we are done */
    if (!fence && 0 == mode) {
        return PMIX_OPERATION_SUCCEEDED;
    }

    /* PRTE WAY, We do it the PMIX way
    cd = PRTE_NEW(prte_pmix_mdx_caddy_t);
    cd->infocbfunc = cbfunc;
    cd->cbdata = cbdata;
    cd->mode = mode;
    */
    size_t pmix_group_context_id;
    /* compute the signature of this collective */
    if (NULL != procs) {
        pmix_group_context_id=nprocs * sizeof(pmix_proc_t);
        PMIX_INFO_LOAD(&results[0], PMIX_GROUP_CONTEXT_ID, &pmix_group_context_id, PMIX_SIZE);
        //PMIX_VALUE_LOAD(&results[0].value, &pmix_group_context_id, PMIX_SIZE);
        /*
        cd->sig = PRTE_NEW(prte_grpcomm_signature_t);
        cd->sig->sz = nprocs;
        cd->sig->signature = (pmix_proc_t *) malloc(cd->sig->sz * sizeof(pmix_proc_t));
        memcpy(cd->sig->signature, procs, cd->sig->sz * sizeof(pmix_proc_t));
        */
    }
    //PMIX_DATA_BUFFER_CREATE(cd->buf);
    
    /* For now only PMIX_GROUP_CONTEXT_ID */
    //bo=NULL; 
    /* if they provided us with a data blob, send it along */
    //if (NULL != bo) {
        /* We don't own the byte_object and so we have to
         * copy it here */
    /*    rc = PMIx_Data_embed(cd->buf, bo);
        if (PMIX_SUCCESS != rc) {
            PMIX_ERROR_LOG(rc);
        }
    }*/
    pmix_info_caddy_t *cd = PMIX_NEW(pmix_info_caddy_t);
    cd->info = results;
    cd->ninfo = 1;
    cbfunc(PMIX_SUCCESS, results, 1, cbdata, pmix_relfn, (void*)cd);   

    /* NO! We don't do this here. Go home PRRTE*/
    /*
    if (PRTE_SUCCESS != (rc = prte_grpcomm.allgather(cd->sig, cd->buf, mode, group_release, cd))) {
        PRTE_ERROR_LOG(rc);
        PRTE_RELEASE(cd);
        return PMIX_ERROR;
    }
    */
    return PMIX_SUCCESS;
}




static pmix_status_t query_fn(pmix_proc_t *proct, pmix_query_t *queries, size_t nqueries,
                              pmix_info_cbfunc_t cbfunc, void *cbdata)
{
    size_t n, k, q, l;
    pmix_info_t *info;
    pmix_status_t rc=PMIX_SUCCESS;

    pmix_output(2, "SERVER: QUERY");

    if (NULL == cbfunc) {
        return PMIX_ERROR;
    }
    /* keep this simple */
    
    for (n = 0; n < nqueries; n++) {
        
        int num_keys, num_keys_not_found =0;
        int nqual=queries->nqual;
        PMIX_ARGV_COUNT(num_keys, queries[n].keys);
        PMIX_INFO_CREATE(info, num_keys + nqual);
        for(q=0; q<queries[n].nqual; q++){
            PMIX_INFO_XFER(&info[q], &queries[n].qualifiers[q]);
        }
        for(k = nqual, l=0; l < num_keys; k++, l++){
            (void) strncpy(info[k].key, queries[n].keys[l], PMIX_MAX_KEYLEN);
            if( 0 == strcmp(queries[n].keys[l], PMIX_QUERY_NUM_PSETS)){
               
                PMIX_VALUE_LOAD(&info[k].value, (void**)&num_psets, PMIX_UINT32);
                
            }else if( 0 == strcmp(queries[n].keys[l], PMIX_QUERY_PSET_MEMBERSHIP) ){
                char pset_name[PMIX_MAX_KEYLEN];
                pmix_pset_t *pset;
                pmix_data_array_t data;
                pmix_proc_t *ptr;
                int flag=0;
                for(q=0; q<queries[n].nqual; q++){
                    if(0 == strcmp(queries[n].qualifiers[q].key, PMIX_PSET_NAME) ){
                        PMIX_LIST_FOREACH(pset, &pmix_server_globals.psets, pmix_pset_t){
                            
                            if(0 == strcmp(pset->name, queries[n].qualifiers[q].value.data.string)){
                                flag=1;
                                break;
                            }
                        }
                        if(flag==1){
                            PMIX_DATA_ARRAY_CONSTRUCT(&data, pset->nmembers, PMIX_PROC);
                            ptr=(pmix_proc_t*)data.array;
                            size_t proc;
                            for(proc=0;proc<pset->nmembers;proc++){
                                PMIX_PROC_LOAD(&ptr[proc],pset->members[proc].nspace, pset->members[proc].rank);
                            }
                        }
                    }
                }
                if(flag==0){
                    rc=PMIX_ERR_NOT_FOUND;
                    num_keys_not_found++;
                    continue;
                }else{
                    PMIX_VALUE_LOAD(&info[k].value, (void**)&data, PMIX_DATA_ARRAY);
                    PMIX_DATA_ARRAY_DESTRUCT(&data);
                }
                
            }else if( 0 == strcmp(queries[n].keys[l], "PMIX_RC_TAG")){
                pthread_mutex_lock(&rc_mutex);
                if(strlen(rc_tag)==0){
                    pthread_mutex_unlock(&rc_mutex);
                    rc=PMIX_ERR_NOT_FOUND;
                    num_keys_not_found++;
                    continue;
                }
                info[k].value.type = PMIX_STRING;
                info[k].value.data.string=rc_tag;
                pthread_mutex_unlock(&rc_mutex);

            }else if( 0 == strcmp(queries[n].keys[l], "PMIX_RC_TYPE")){
                pthread_mutex_lock(&rc_mutex);
                if(strlen(rc_op)==0){
                    pthread_mutex_unlock(&rc_mutex);
                    rc=PMIX_ERR_NOT_FOUND;
                    num_keys_not_found++;
                    continue;
                }
                uint8_t rc_type=(0 == strcmp(rc_op, "ADD")) ? PMIX_RES_CHANGE_ADD : PMIX_RES_CHANGE_SUB;
                PMIX_VALUE_LOAD(&info[k].value, (void**)&rc_type, PMIX_UINT8);
                pthread_mutex_unlock(&rc_mutex);

            }else if( 0 == strcmp(queries[n].keys[l], "PMIX_RC_PSET")){
                pthread_mutex_lock(&rc_mutex);
                if(strlen(rc_pset)==0){
                    pthread_mutex_unlock(&rc_mutex);
                    rc=PMIX_ERR_NOT_FOUND;
                    num_keys_not_found++;
                    continue;
                }
                PMIX_VALUE_LOAD(&info[k].value, rc_pset, PMIX_STRING);
                //info[k].value.type = PMIX_STRING;
                //info[k].value.data.string=rc_pset;
                pthread_mutex_unlock(&rc_mutex);

            }else if( 0 == strcmp(queries[n].keys[l], "PMIX_RC_PORT")){
                pthread_mutex_lock(&rc_mutex);
                info[k].value.type = PMIX_STRING;
                info[k].value.data.string="test port";
                pthread_mutex_unlock(&rc_mutex);

            }else if( 0 == strcmp(queries[n].keys[l], PMIX_QUERY_PSET_NAMES)){
                size_t npsets= pmix_list_get_size(&pmix_server_globals.psets);
                char *names;
                size_t old_size=0;
                size_t pset_index=0;
                pmix_pset_t *pset;
                PMIX_LIST_FOREACH(pset, &pmix_server_globals.psets, pmix_pset_t){
                    if(old_size==0){
                        names=strdup(pset->name);
                        
                    }else{
                        names=realloc(names, strlen(names)+strlen(pset->name)+2);
                        names[old_size]=',';
                        strcpy(&names[old_size+1], pset->name);
                    }
                    old_size=strlen(names);
                }
                
                PMIX_VALUE_LOAD(&info[k].value, names, PMIX_STRING);
                if(NULL != names){
                    free(names);
                }
            }
            
            else{
                info[k].value.type = PMIX_STRING;
                if (0 > asprintf(&info[n].value.data.string, "%d", (int) n)) {
                    return PMIX_ERROR;
                }
            }
        }
        pmix_info_caddy_t *cd = PMIX_NEW(pmix_info_caddy_t);
        cd->info = info;
        cd->ninfo = num_keys + nqual;
        cbfunc(rc, info, num_keys + nqual - num_keys_not_found, cbdata, pmix_relfn, (void*)cd);
        //PMIX_INFO_FREE(info, num_keys-num_keys_not_found);
        
    }
    
    return PMIX_SUCCESS;
}

static void tool_connect_fn(pmix_info_t *info, size_t ninfo, pmix_tool_connection_cbfunc_t cbfunc,
                            void *cbdata)
{
    pmix_proc_t proc;

    pmix_output(0, "SERVER: TOOL CONNECT");

    /* just pass back an arbitrary nspace */
    (void) strncpy(proc.nspace, "TOOL", PMIX_MAX_NSLEN);
    proc.rank = 0;

    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, &proc, cbdata);
    }
}

static void log_fn(const pmix_proc_t *client, const pmix_info_t data[], size_t ndata,
                   const pmix_info_t directives[], size_t ndirs, pmix_op_cbfunc_t cbfunc,
                   void *cbdata)
{
    pmix_output(0, "SERVER: LOG");

    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, cbdata);
    }
}

static void wait_signal_callback(int fd, short event, void *arg)
{
    pmix_event_t *sig = (pmix_event_t *) arg;
    int status;
    pid_t pid;
    wait_tracker_t *t2;

    if (SIGCHLD != pmix_event_get_signal(sig)) {
        return;
    }

    /* we can have multiple children leave but only get one
     * sigchild callback, so reap all the waitpids until we
     * don't get anything valid back */
    while (1) {
        pid = waitpid(-1, &status, WNOHANG);
        if (-1 == pid && EINTR == errno) {
            /* try it again */
            continue;
        }
        /* if we got garbage, then nothing we can do */
        if (pid <= 0) {
            return;
        }

        /* we are already in an event, so it is safe to access the list */
        PMIX_LIST_FOREACH (t2, &children, wait_tracker_t) {
            if (pid == t2->pid) {
                /* found it! */
                //--wakeup;
                printf("signal");
                break;
            }
        }
    }
}
