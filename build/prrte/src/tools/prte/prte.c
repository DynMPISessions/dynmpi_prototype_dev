/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil -*- */
/*
 * Copyright (c) 2004-2010 The Trustees of Indiana University and Indiana
 *                         University Research and Technology
 *                         Corporation.  All rights reserved.
 * Copyright (c) 2004-2021 The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 * Copyright (c) 2004-2005 High Performance Computing Center Stuttgart,
 *                         University of Stuttgart.  All rights reserved.
 * Copyright (c) 2004-2005 The Regents of the University of California.
 *                         All rights reserved.
 * Copyright (c) 2006-2020 Cisco Systems, Inc.  All rights reserved
 * Copyright (c) 2007-2009 Sun Microsystems, Inc. All rights reserved.
 * Copyright (c) 2007-2017 Los Alamos National Security, LLC.  All rights
 *                         reserved.
 * Copyright (c) 2013-2020 Intel, Inc.  All rights reserved.
 * Copyright (c) 2015-2019 Research Organization for Information Science
 *                         and Technology (RIST).  All rights reserved.
 * Copyright (c) 2020      Geoffroy Vallee. All rights reserved.
 * Copyright (c) 2020      IBM Corporation.  All rights reserved.
 * Copyright (c) 2021      Nanook Consulting.  All rights reserved.
 * Copyright (c) 2021      Amazon.com, Inc. or its affiliates.  All Rights
 *                         reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "prte_config.h"
#include "src/include/constants.h"
#include "src/include/version.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef HAVE_STRINGS_H
#    include <strings.h>
#endif /* HAVE_STRINGS_H */
#ifdef HAVE_UNISTD_H
#    include <unistd.h>
#endif
#ifdef HAVE_SYS_PARAM_H
#    include <sys/param.h>
#endif
#include <ctype.h>
#include <errno.h>
#include <signal.h>
#ifdef HAVE_SYS_TYPES_H
#    include <sys/types.h>
#endif /* HAVE_SYS_TYPES_H */
#ifdef HAVE_SYS_WAIT_H
#    include <sys/wait.h>
#endif /* HAVE_SYS_WAIT_H */
#ifdef HAVE_SYS_TIME_H
#    include <sys/time.h>
#endif /* HAVE_SYS_TIME_H */
#include <fcntl.h>
#ifdef HAVE_SYS_STAT_H
#    include <sys/stat.h>
#endif
#ifdef HAVE_POLL_H
#    include <poll.h>
#endif

#include "src/event/event-internal.h"
#include "src/mca/base/base.h"
#include "src/mca/prteinstalldirs/prteinstalldirs.h"
#include "src/pmix/pmix-internal.h"
#include "src/threads/mutex.h"
#include "src/util/argv.h"
#include "src/util/basename.h"
#include "src/util/cmd_line.h"
#include "src/util/daemon_init.h"
#include "src/util/fd.h"
#include "src/util/os_path.h"
#include "src/util/output.h"
#include "src/util/path.h"
#include "src/util/printf.h"
#include "src/util/prte_environ.h"
#include "src/util/prte_getcwd.h"
#include "src/util/show_help.h"

#include "src/class/prte_pointer_array.h"
#include "src/runtime/prte_progress_threads.h"

#include "src/mca/errmgr/errmgr.h"
#include "src/mca/ess/base/base.h"
#include "src/mca/plm/plm.h"
#include "src/mca/plm/base/base.h"
#include "src/mca/prteif/prteif.h"
#include "src/mca/rmaps/rmaps_types.h"
#include "src/mca/rml/rml.h"
#include "src/mca/schizo/base/base.h"
#include "src/mca/state/base/base.h"
#include "src/runtime/prte_globals.h"
#include "src/runtime/runtime.h"

#include "prte.h"
#include "src/prted/pmix/pmix_server_internal.h"
#include "src/prted/prted.h"

typedef struct {
    prte_pmix_lock_t lock;
    pmix_status_t status;
    pmix_info_t *info;
    size_t ninfo;
} mylock_t;

pmix_rank_t highest_rank_global = 6; 


static int res_change_cnt=0;
static pmix_nspace_t spawnednspace;
static pmix_proc_t myproc;
static bool signals_set = false;
static bool forcibly_die = false;
static prte_event_t term_handler;
static prte_event_t die_handler;
static prte_event_t epipe_handler;
static int term_pipe[2];
static prte_mutex_t prun_abort_inprogress_lock = PRTE_MUTEX_STATIC_INIT;
static prte_event_t *forward_signals_events = NULL;
static char *mypidfile = NULL;
static bool verbose = false;
static prte_cmd_line_t *prte_cmd_line = NULL;
static bool want_prefix_by_default = (bool) PRTE_WANT_PRTE_PREFIX_BY_DEFAULT;
static void abort_signal_callback(int signal);
static void clean_abort(int fd, short flags, void *arg);
static void signal_forward_callback(int fd, short args, void *cbdata);
static void epipe_signal_callback(int fd, short args, void *cbdata);
static int prep_singleton(const char *name);

static void rchandler(size_t evhdlr_registration_id, pmix_status_t status,
                       const pmix_proc_t *source, pmix_info_t info[], size_t ninfo,
                       pmix_info_t results[], size_t nresults,
                       pmix_event_notification_cbfunc_fn_t cbfunc, void *cbdata);

static void opcbfunc(pmix_status_t status, void *cbdata)
{
    prte_pmix_lock_t *lock = (prte_pmix_lock_t *) cbdata;
    PRTE_ACQUIRE_OBJECT(lock);
    PRTE_PMIX_WAKEUP_THREAD(lock);
}

static void setupcbfunc(pmix_status_t status, pmix_info_t info[], size_t ninfo,
                        void *provided_cbdata, pmix_op_cbfunc_t cbfunc, void *cbdata)
{
    mylock_t *mylock = (mylock_t *) provided_cbdata;
    size_t n;

    if (NULL != info) {
        mylock->ninfo = ninfo;
        PMIX_INFO_CREATE(mylock->info, mylock->ninfo);
        /* cycle across the provided info */
        for (n = 0; n < ninfo; n++) {
            PMIX_INFO_XFER(&mylock->info[n], &info[n]);
        }
    } else {
        mylock->info = NULL;
        mylock->ninfo = 0;
    }
    mylock->status = status;

    /* release the caller */
    if (NULL != cbfunc) {
        cbfunc(PMIX_SUCCESS, cbdata);
    }

    PRTE_PMIX_WAKEUP_THREAD(&mylock->lock);
}

static void toolcbfunc(pmix_status_t status, pmix_proc_t *name, void *cbdata)
{
    mylock_t *mylock = (mylock_t *) cbdata;

    mylock->status = status;
    PRTE_PMIX_WAKEUP_THREAD(&mylock->lock);
}

static void spcbfunc(pmix_status_t status, char nspace[], void *cbdata)
{
    prte_pmix_lock_t *lock = (prte_pmix_lock_t *) cbdata;

    PRTE_ACQUIRE_OBJECT(lock);
    lock->status = status;
    if (PMIX_SUCCESS == status) {
        lock->msg = strdup(nspace);
    }
    PRTE_PMIX_WAKEUP_THREAD(lock);
}

static int wait_pipe[2];

static int wait_dvm(pid_t pid)
{
    char reply;
    int rc;
    int status;

    close(wait_pipe[1]);
    do {
        rc = read(wait_pipe[0], &reply, 1);
    } while (0 > rc && EINTR == errno);

    if (1 == rc && 'K' == reply) {
        return 0;
    } else if (0 == rc) {
        waitpid(pid, &status, 0);
        if (WIFEXITED(status)) {
            return WEXITSTATUS(status);
        }
    }
    return 255;
}

static void setup_sighandler(int signal, prte_event_t *ev, prte_event_cbfunc_t cbfunc)
{
    prte_event_signal_set(prte_event_base, ev, signal, cbfunc, ev);
    prte_event_signal_add(ev, NULL);
}

static prte_cmd_line_init_t cmd_line_init[] = {
    /* override personality */
    {'\0', "personality", 1, PRTE_CMD_LINE_TYPE_STRING, "Specify the personality to be used",
     PRTE_CMD_LINE_OTYPE_DVM},

    /* End of list */
    {'\0', NULL, 0, PRTE_CMD_LINE_TYPE_NULL, NULL}};

int prte(int argc, char *argv[])
{
    int rc = 1, i, j;
    char *param, *ptr, *tpath, *fullpath;
    prte_pmix_lock_t lock;
    prte_list_t apps;
    prte_pmix_app_t *app;
    pmix_info_t *iptr, info;
    pmix_status_t ret;
    bool flag;
    size_t n, ninfo, param_len;
    pmix_app_t *papps;
    size_t napps;
    mylock_t mylock;
    prte_value_t *pval;
    uint32_t ui32;
    char **pargv;
    int pargc;
    prte_job_t *jdata;
    prte_app_context_t *dapp;
    bool proxyrun = false;
    void *jinfo;
    pmix_proc_t pname, parent;
    pmix_value_t *val;
    pmix_data_array_t darray;
    char **hostfiles = NULL;
    char **hosts = NULL;
    bool donotlaunch = false;
    prte_schizo_base_module_t *schizo;
    prte_ess_base_signal_t *sig;
    char **targv;

    /* init the globals */
    PRTE_CONSTRUCT(&apps, prte_list_t);
    /* init the tiny part of PRTE we use */
    prte_init_util(PRTE_PROC_MASTER);

    fullpath = prte_find_absolute_path(argv[0]);
    prte_tool_basename = prte_basename(argv[0]);
    pargc = argc;
    pargv = prte_argv_copy(argv);

    /** setup callbacks for abort signals - from this point
     * forward, we need to abort in a manner that allows us
     * to cleanup. However, we cannot directly use libevent
     * to trap these signals as otherwise we cannot respond
     * to them if we are stuck in an event! So instead use
     * the basic POSIX trap functions to handle the signal,
     * and then let that signal handler do some magic to
     * avoid the hang
     *
     * NOTE: posix traps don't allow us to do anything major
     * in them, so use a pipe tied to a libevent event to
     * reach a "safe" place where the termination event can
     * be created
     */
    if (0 != (rc = pipe(term_pipe))) {
        exit(1);
    }
    /* setup an event to attempt normal termination on signal */
    rc = prte_event_base_open();
    if (PRTE_SUCCESS != rc) {
        fprintf(stderr, "Unable to initialize event library\n");
        exit(1);
    }
    prte_event_set(prte_event_base, &term_handler, term_pipe[0], PRTE_EV_READ, clean_abort, NULL);
    prte_event_add(&term_handler, NULL);

    /* Set both ends of this pipe to be close-on-exec so that no
     children inherit it */
    if (prte_fd_set_cloexec(term_pipe[0]) != PRTE_SUCCESS
        || prte_fd_set_cloexec(term_pipe[1]) != PRTE_SUCCESS) {
        fprintf(stderr, "unable to set the pipe to CLOEXEC\n");
        prte_progress_thread_finalize(NULL);
        exit(1);
    }

    /* setup callback for SIGPIPE */
    setup_sighandler(SIGPIPE, &epipe_handler, epipe_signal_callback);

    /* point the signal trap to a function that will activate that event */
    signal(SIGTERM, abort_signal_callback);
    signal(SIGINT, abort_signal_callback);
    signal(SIGHUP, abort_signal_callback);

    /* because we have to use the schizo framework prior to parsing the
     * incoming argv for cmd line options, do a hacky search to support
     * passing of options (e.g., verbosity) for schizo */
    for (i = 1; NULL != argv[i]; i++) {
        if (0 == strcmp(argv[i], "--prtemca") || 0 == strcmp(argv[i], "--mca")) {
            if (0 == strncmp(argv[i + 1], "schizo", 6)) {
                prte_asprintf(&param, "PRTE_MCA_%s", argv[i + 1]);
                prte_setenv(param, argv[i + 2], true, &environ);
                free(param);
                i += 2;
            }
        }
    }

    /* open the SCHIZO framework */
    if (PRTE_SUCCESS
        != (rc = prte_mca_base_framework_open(&prte_schizo_base_framework,
                                              PRTE_MCA_BASE_OPEN_DEFAULT))) {
        PRTE_ERROR_LOG(rc);
        return rc;
    }

    if (PRTE_SUCCESS != (rc = prte_schizo_base_select())) {
        PRTE_ERROR_LOG(rc);
        return rc;
    }

    /* look for any personality specification */
    ptr = NULL;
    for (i = 0; NULL != argv[i]; i++) {
        if (0 == strcmp(argv[i], "--personality")) {
            ptr = argv[i + 1];
            break;
        }
    }
    if (NULL == ptr) {
        ptr = fullpath;
    }

    /* detect if we are running as a proxy and select the active
     * schizo module for this tool */
    schizo = prte_schizo.detect_proxy(ptr);
    if (NULL == schizo) {
        prte_show_help("help-schizo-base.txt", "no-proxy", true, prte_tool_basename, ptr);
        return 1;
    }
    if (0 != strcmp(schizo->name, "prte")) {
        proxyrun = true;
    } else {
        /* if we are using the "prte" personality, but we
         * are not actually running as "prte" or are actively
         * testing the proxy capability , then we are acting
         * as a proxy */
        if (0 != strcmp(prte_tool_basename, "prte") || prte_schizo_base.test_proxy_launch) {
            proxyrun = true;
        }
    }

    /* setup the cmd line - this is specific to the proxy */
    prte_cmd_line = PRTE_NEW(prte_cmd_line_t);
    if (PRTE_SUCCESS != (rc = schizo->define_cli(prte_cmd_line))) {
        PRTE_ERROR_LOG(rc);
        return rc;
    }
    /* add any prte-specific options */
    if (PRTE_SUCCESS != (rc = prte_cmd_line_add(prte_cmd_line, cmd_line_init))) {
        PRTE_ERROR_LOG(rc);
        return rc;
    }

    /* handle deprecated options */
    if (PRTE_SUCCESS != (rc = schizo->parse_deprecated_cli(prte_cmd_line, &pargc, &pargv))) {
        if (PRTE_OPERATION_SUCCEEDED == rc) {
            /* the cmd line was restructured - show them the end result */
            param = prte_argv_join(pargv, ' ');
            fprintf(stderr, "\n******* Corrected cmd line: %s\n\n\n", param);
            free(param);
        } else {
            return rc;
        }
    }

    /* parse the result to get values - this will not include MCA params */
    if (PRTE_SUCCESS != (rc = prte_cmd_line_parse(prte_cmd_line, true, false, pargc, pargv))) {
        if (PRTE_ERR_SILENT != rc) {
            fprintf(stderr, "%s: command line error (%s)\n", prte_tool_basename, prte_strerror(rc));
        }
        return rc;
    }

    /* check if we are running as root - if we are, then only allow
     * us to proceed if the allow-run-as-root flag was given. Otherwise,
     * exit with a giant warning message
     */
    if (0 == geteuid()) {
        schizo->allow_run_as_root(prte_cmd_line); // will exit us if not allowed
    }

    /* if we were given a keepalive pipe, set up to monitor it now */
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "keepalive", 0, 0))) {
        prte_event_set(prte_event_base, &die_handler, pval->value.data.integer, PRTE_EV_READ,
                       clean_abort, NULL);
        prte_event_add(&die_handler, NULL);
        prte_fd_set_cloexec(pval->value.data.integer); // don't let children inherit this
    }

    /* let the schizo components take a pass at it to get the MCA params */
    if (PRTE_SUCCESS != (rc = schizo->parse_cli(pargc, 0, pargv, NULL))) {
        if (PRTE_ERR_SILENT != rc) {
            fprintf(stderr, "%s: command line error (%s)\n", prte_tool_basename, prte_strerror(rc));
        }
        return rc;
    }

    /* check command line sanity - ensure there aren't multiple instances of
     * options where there should be only one */
    rc = schizo->check_sanity(prte_cmd_line);
    if (PRTE_SUCCESS != rc) {
        if (PRTE_ERR_SILENT != rc) {
            fprintf(stderr, "%s: command line error (%s)\n", prte_tool_basename, prte_strerror(rc));
        }
        param = prte_argv_join(pargv, ' ');
        fprintf(stderr, "\n******* Cmd line: %s\n\n\n", param);
        free(param);
        return rc;
    }

    if (prte_cmd_line_is_taken(prte_cmd_line, "verbose")) {
        verbose = true;
    }

    /* see if print version is requested. Do this before
     * check for help so that --version --help works as
     * one might expect. */
    if (prte_cmd_line_is_taken(prte_cmd_line, "version")) {
        if (proxyrun) {
            fprintf(stdout, "%s (%s) %s\n\nReport bugs to %s\n", prte_tool_basename,
                    PRTE_PROXY_PACKAGE_NAME, PRTE_PROXY_VERSION_STRING, PRTE_PROXY_BUGREPORT);
        } else {
            fprintf(stdout, "%s (%s) %s\n\nReport bugs to %s\n", prte_tool_basename,
                    "PMIx Reference RunTime Environment", PRTE_VERSION, PACKAGE_BUGREPORT);
        }
        exit(0);
    }

    /* Check for help request */
    if (prte_cmd_line_is_taken(prte_cmd_line, "help")) {
        char *str, *args = NULL;
        args = prte_cmd_line_get_usage_msg(prte_cmd_line, false);
        str = prte_show_help_string("help-prun.txt", "prun:usage", false, prte_tool_basename,
                                    "PRTE", PRTE_VERSION, prte_tool_basename, args,
                                    PACKAGE_BUGREPORT);
        if (NULL != str) {
            printf("%s", str);
            free(str);
        }
        free(args);

        /* If someone asks for help, that should be all we do */
        exit(0);
    }

    /* set debug flags */
    prte_debug_flag = prte_cmd_line_is_taken(prte_cmd_line, "debug");
    prte_debug_daemons_flag = prte_cmd_line_is_taken(prte_cmd_line, "debug-daemons");
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "debug-verbose", 0, 0))) {
        prte_debug_verbosity = pval->value.data.integer;
    }
    prte_debug_daemons_file_flag = prte_cmd_line_is_taken(prte_cmd_line, "debug-daemons-file");
    if (prte_debug_daemons_file_flag) {
        prte_debug_daemons_flag = true;
    }
    prte_leave_session_attached = prte_cmd_line_is_taken(prte_cmd_line, "leave-session-attached");
    /* if any debug level is set, ensure we output debug level dumps */
    if (prte_debug_flag || prte_debug_daemons_flag || prte_leave_session_attached) {
        prte_devel_level_output = true;
    }

    /* detach from controlling terminal
     * otherwise, remain attached so output can get to us
     */
    if (!prte_debug_flag && !prte_debug_daemons_flag
        && prte_cmd_line_is_taken(prte_cmd_line, "daemonize")) {
        pipe(wait_pipe);
        prte_state_base_parent_fd = wait_pipe[1];
        prte_daemon_init_callback(NULL, wait_dvm);
        close(wait_pipe[0]);
    } else {
#if defined(HAVE_SETSID)
        /* see if we were directed to separate from current session */
        if (prte_cmd_line_is_taken(prte_cmd_line, "set-sid")) {
            setsid();
        }
#endif
    }

    if (prte_cmd_line_is_taken(prte_cmd_line, "no-ready-msg")) {
        prte_state_base_ready_msg = false;
    }

    if (prte_cmd_line_is_taken(prte_cmd_line, "system-server")) {
        /* we should act as system-level PMIx server */
        prte_setenv("PRTE_MCA_pmix_system_server", "1", true, &environ);
    }
    /* always act as session-level PMIx server */
    prte_setenv("PRTE_MCA_pmix_session_server", "1", true, &environ);
    /* if we were asked to report a uri, set the MCA param to do so */
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "report-uri", 0, 0))) {
        prte_setenv("PMIX_MCA_ptl_base_report_uri", pval->value.data.string, true, &environ);
    }
    /* don't aggregate help messages as that will apply job-to-job */
    prte_setenv("PRTE_MCA_prte_base_help_aggregate", "0", true, &environ);

    /* if we are supporting a singleton, push its ID into the environ
     * so it can get picked up and registered by server init */
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "singleton", 0, 0))) {
        prte_setenv("PMIX_MCA_singleton", pval->value.data.string, true, &environ);
    }

    /* Setup MCA params */
    prte_register_params();

    /* save the environment for launch purposes. This MUST be
     * done so that we can pass it to any local procs we
     * spawn - otherwise, those local procs won't see any
     * non-MCA envars were set in the enviro prior to calling
     * prun
     */
    prte_launch_environ = prte_argv_copy(environ);

    /* setup PRTE infrastructure */
    if (PRTE_SUCCESS != (ret = prte_init(&pargc, &pargv, PRTE_PROC_MASTER))) {
        PRTE_ERROR_LOG(ret);
        return ret;
    }
    /* get my proc ID */
    ret = PMIx_Get(NULL, PMIX_PROCID, NULL, 0, &val);
    if (PMIX_SUCCESS != ret) {
        PMIX_ERROR_LOG(ret);
        PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
        goto DONE;
    }
    memcpy(&myproc, val->data.proc, sizeof(pmix_proc_t));
    PMIX_VALUE_RELEASE(val);

    /** setup callbacks for signals we should forward */
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "forward-signals", 0, 0))) {
        param = pval->value.data.string;
    } else {
        param = NULL;
    }
    if (PRTE_SUCCESS != (rc = prte_ess_base_setup_signals(param))) {
        PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
        goto DONE;
    }
    if (0 < (i = prte_list_get_size(&prte_ess_base_signals))) {
        forward_signals_events = (prte_event_t *) malloc(sizeof(prte_event_t) * i);
        if (NULL == forward_signals_events) {
            ret = PRTE_ERR_OUT_OF_RESOURCE;
            PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
            goto DONE;
        }
        i = 0;
        PRTE_LIST_FOREACH(sig, &prte_ess_base_signals, prte_ess_base_signal_t)
        {
            setup_sighandler(sig->signal, forward_signals_events + i, signal_forward_callback);
            ++i;
        }
    }
    signals_set = true;

    /* if we are supporting a singleton, add it to our jobs */
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "singleton", 0, 0))) {
        rc = prep_singleton(pval->value.data.string);
        if (PRTE_SUCCESS != ret) {
            PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
            goto DONE;
        }
    }

    /* check for launch directives in case we were launched by a
     * tool wanting to direct our operation - this needs to be
     * done prior to starting the DVM as it may include instructions
     * on the daemon executable, the fork/exec agent to be used by
     * the daemons, or other directives impacting the DVM itself. */
    PMIX_LOAD_PROCID(&pname, myproc.nspace, PMIX_RANK_WILDCARD);
    PMIX_INFO_LOAD(&info, PMIX_OPTIONAL, NULL, PMIX_BOOL);
    /*  Have to cycle over directives we support*/
    ret = PMIx_Get(&pname, PMIX_FORKEXEC_AGENT, &info, 1, &val);
    PMIX_INFO_DESTRUCT(&info);
    if (PMIX_SUCCESS == ret) {
        /* set our fork/exec agent */
        PMIX_VALUE_RELEASE(val);
    }

    /* start the DVM */

    /* get the daemon job object - was created by ess/hnp component */
    if (NULL == (jdata = prte_get_job_data_object(PRTE_PROC_MY_NAME->nspace))) {
        prte_show_help("help-prun.txt", "bad-job-object", true, prte_tool_basename);
        PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
        goto DONE;
    }
    /* ess/hnp also should have created a daemon "app" */
    if (NULL == (dapp = (prte_app_context_t *) prte_pointer_array_get_item(jdata->apps, 0))) {
        prte_show_help("help-prun.txt", "bad-app-object", true, prte_tool_basename);
        PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
        goto DONE;
    }

    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "map-by", 0, 0))) {
        if (NULL != strcasestr(pval->value.data.string, "DONOTLAUNCH")) {
            prte_set_attribute(&jdata->attributes, PRTE_JOB_DO_NOT_LAUNCH, PRTE_ATTR_GLOBAL, NULL,
                               PMIX_BOOL);
            donotlaunch = true;
        }
    }

    /* Did the user specify a prefix, or want prefix by default? */
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "prefix", 0, 0))
        || want_prefix_by_default) {
        if (NULL != pval) {
            param = strdup(pval->value.data.string);
        } else {
            /* --enable-prun-prefix-default was given to prun */
            param = strdup(prte_install_dirs.prefix);
        }
        /* "Parse" the param, aka remove superfluous path_sep. */
        param_len = strlen(param);
        while (0 == strcmp(PRTE_PATH_SEP, &(param[param_len - 1]))) {
            param[param_len - 1] = '\0';
            param_len--;
            if (0 == param_len) {
                prte_show_help("help-prun.txt", "prun:empty-prefix", true, prte_tool_basename,
                               prte_tool_basename);
                PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
                goto DONE;
            }
        }
        prte_set_attribute(&dapp->attributes, PRTE_APP_PREFIX_DIR, PRTE_ATTR_GLOBAL, param,
                           PMIX_STRING);
        free(param);
    } else {
        /* Check if called with fully-qualified path to prte.
           (Note: Put this second so can override with --prefix (above). */
        tpath = NULL;
        if ('/' == argv[0][0]) {
            char *tmp_basename = NULL;
            tpath = prte_dirname(argv[0]);

            if (NULL != tpath) {
                /* Quick sanity check to ensure we got
                   something/bin/<exec_name> and that the installation
                   tree is at least more or less what we expect it to
                   be */
                tmp_basename = prte_basename(tpath);
                if (0 == strcmp("bin", tmp_basename)) {
                    char *tmp = tpath;
                    tpath = prte_dirname(tmp);
                    free(tmp);
                } else {
                    free(tpath);
                    tpath = NULL;
                }
                free(tmp_basename);
            }
            prte_set_attribute(&dapp->attributes, PRTE_APP_PREFIX_DIR, PRTE_ATTR_GLOBAL, tpath,
                               PMIX_STRING);
        }
    }

    /* setup to listen for commands sent specifically to me, even though I would probably
     * be the one sending them! Unfortunately, since I am a participating daemon,
     * there are times I need to send a command to "all daemons", and that means *I* have
     * to receive it too
     */
    prte_rml.recv_buffer_nb(PRTE_NAME_WILDCARD, PRTE_RML_TAG_DAEMON, PRTE_RML_PERSISTENT,
                            prte_daemon_recv, NULL);

    /* setup to capture job-level info */
    PMIX_INFO_LIST_START(jinfo);

    /* see if we ourselves were spawned by someone */
    ret = PMIx_Get(&prte_process_info.myproc, PMIX_PARENT_ID, NULL, 0, &val);
    if (PMIX_SUCCESS == ret) {
        PMIX_LOAD_PROCID(&parent, val->data.proc->nspace, val->data.proc->rank);
        PMIX_VALUE_RELEASE(val);
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_REQUESTOR_IS_TOOL, NULL, PMIX_BOOL);
        /* record that this tool is connected to us */
        PRTE_PMIX_CONSTRUCT_LOCK(&mylock.lock);
        PMIX_INFO_CREATE(iptr, 2);
        PMIX_INFO_LOAD(&iptr[0], PMIX_NSPACE, parent.nspace, PMIX_STRING);
        PMIX_INFO_LOAD(&iptr[1], PMIX_RANK, &parent.rank, PMIX_PROC_RANK);
        pmix_tool_connected_fn(iptr, 2, toolcbfunc, &mylock);
        /* we have to cycle the event library here so we can process
         * this request */
        while (prte_event_base_active && mylock.lock.active) {
            prte_event_loop(prte_event_base, PRTE_EVLOOP_ONCE);
        }
        PRTE_ACQUIRE_OBJECT(&mylock.lock);
        PMIX_INFO_FREE(iptr, 2);
        PRTE_PMIX_DESTRUCT_LOCK(&mylock.lock);
    } else {
        PMIX_LOAD_PROCID(&parent, prte_process_info.myproc.nspace, prte_process_info.myproc.rank);
    }

    /* default to a persistent DVM */
    prte_persistent = true;

    /* if we are told to daemonize, then we cannot have apps */
    if (!prte_cmd_line_is_taken(prte_cmd_line, "daemonize")) {
        /* see if they want to run an application - let's parse
         * the cmd line to get it */
        rc = prte_parse_locals(prte_cmd_line, &apps, pargc, pargv, &hostfiles, &hosts);

        /* did they provide an app? */
        if (PMIX_SUCCESS != rc || 0 == prte_list_get_size(&apps)) {
            if (proxyrun) {
                prte_show_help("help-prun.txt", "prun:executable-not-specified", true,
                               prte_tool_basename, prte_tool_basename);
                PRTE_UPDATE_EXIT_STATUS(rc);
                goto DONE;
            }
            /* nope - just need to wait for instructions */
        } else {
            /* they did provide an app - this is only allowed
             * when running as a proxy! */
            if (!proxyrun) {
                prte_show_help("help-prun.txt", "prun:executable-incorrectly-given", true,
                               prte_tool_basename, prte_tool_basename);
                PRTE_UPDATE_EXIT_STATUS(rc);
                goto DONE;
            }
            /* mark that we are not a persistent DVM */
            prte_persistent = false;
        }
    }

    /* add any hostfile directives to the daemon job */
    if (prte_persistent) {
        if (0 < (j = prte_cmd_line_get_ninsts(prte_cmd_line, "hostfile"))) {
            if (1 < j) {
                prte_show_help("help-prun.txt", "prun:multiple-hostfiles", true, prte_tool_basename,
                               NULL);
                PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
                goto DONE;
            } else {
                pval = prte_cmd_line_get_param(prte_cmd_line, "hostfile", 0, 0);
                prte_set_attribute(&dapp->attributes, PRTE_APP_HOSTFILE, PRTE_ATTR_GLOBAL,
                                   pval->value.data.string, PMIX_STRING);
            }
        }
        if (0 < (j = prte_cmd_line_get_ninsts(prte_cmd_line, "machinefile"))) {
            if (1 < j
                || prte_get_attribute(&dapp->attributes, PRTE_APP_HOSTFILE, NULL, PMIX_STRING)) {
                prte_show_help("help-prun.txt", "prun:multiple-hostfiles", true, prte_tool_basename,
                               NULL);
                PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
                goto DONE;
            } else {
                pval = prte_cmd_line_get_param(prte_cmd_line, "machinefile", 0, 0);
                prte_set_attribute(&dapp->attributes, PRTE_APP_HOSTFILE, PRTE_ATTR_GLOBAL,
                                   pval->value.data.string, PMIX_STRING);
            }
        }

        /* Did the user specify any hosts? */
        if (0 < (j = prte_cmd_line_get_ninsts(prte_cmd_line, "host"))) {
            char **targ = NULL, *tval;
            for (i = 0; i < j; ++i) {
                pval = prte_cmd_line_get_param(prte_cmd_line, "host", i, 0);
                prte_argv_append_nosize(&targ, pval->value.data.string);
            }
            tval = prte_argv_join(targ, ',');
            prte_set_attribute(&dapp->attributes, PRTE_APP_DASH_HOST, PRTE_ATTR_GLOBAL, tval,
                               PMIX_STRING);
            prte_argv_free(targ);
            free(tval);
        }
    } else {
        /* the directives will be in the app(s) */
        if (NULL != hostfiles) {
            char *tval;
            tval = prte_argv_join(hostfiles, ',');
            prte_set_attribute(&dapp->attributes, PRTE_APP_HOSTFILE, PRTE_ATTR_GLOBAL, tval,
                               PMIX_STRING);
            free(tval);
            prte_argv_free(hostfiles);
        }
        if (NULL != hosts) {
            char *tval;
            tval = prte_argv_join(hosts, ',');
            prte_set_attribute(&dapp->attributes, PRTE_APP_DASH_HOST, PRTE_ATTR_GLOBAL, tval,
                               PMIX_STRING);
            free(tval);
            prte_argv_free(hosts);
        }
    }
    /* pickup any relevant envars that need to go on the DVM cmd line */
    rc = prte_schizo.parse_env(prte_cmd_line, environ, &pargv, true);
    if (PRTE_SUCCESS != rc) {
        PRTE_UPDATE_EXIT_STATUS(rc);
        goto DONE;
    }

    /* spawn the DVM - we skip the initial steps as this
     * isn't a user-level application */
    PRTE_ACTIVATE_JOB_STATE(jdata, PRTE_JOB_STATE_ALLOCATE);

    /* we need to loop the event library until the DVM is alive */
    while (prte_event_base_active && !prte_dvm_ready) {
        prte_event_loop(prte_event_base, PRTE_EVLOOP_ONCE);
    }
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "report-pid", 0, 0))) {
        /* if the string is a "-", then output to stdout */
        if (0 == strcmp(pval->value.data.string, "-")) {
            fprintf(stdout, "%lu\n", (unsigned long) getpid());
        } else if (0 == strcmp(pval->value.data.string, "+")) {
            /* output to stderr */
            fprintf(stderr, "%lu\n", (unsigned long) getpid());
        } else {
            char *leftover;
            int outpipe;
            /* see if it is an integer pipe */
            leftover = NULL;
            outpipe = strtol(pval->value.data.string, &leftover, 10);
            if (NULL == leftover || 0 == strlen(leftover)) {
                /* stitch together the var names and URI */
                prte_asprintf(&leftover, "%lu", (unsigned long) getpid());
                /* output to the pipe */
                rc = prte_fd_write(outpipe, strlen(leftover) + 1, leftover);
                free(leftover);
                close(outpipe);
            } else {
                /* must be a file */
                FILE *fp;
                fp = fopen(pval->value.data.string, "w");
                if (NULL == fp) {
                    prte_output(0, "Impossible to open the file %s in write mode\n",
                                pval->value.data.string);
                    PRTE_UPDATE_EXIT_STATUS(1);
                    goto DONE;
                }
                /* output my PID */
                fprintf(fp, "%lu\n", (unsigned long) getpid());
                fclose(fp);
                mypidfile = strdup(pval->value.data.string);
            }
        }
    }

    if (prte_persistent) {
        PMIX_INFO_LIST_RELEASE(jinfo);
        goto proceed;
    }

    /***** CHECK FOR LAUNCH DIRECTIVES - ADD THEM TO JOB INFO IF FOUND ****/
    PMIX_LOAD_PROCID(&pname, myproc.nspace, PMIX_RANK_WILDCARD);
    PMIX_INFO_LOAD(&info, PMIX_OPTIONAL, NULL, PMIX_BOOL);
    ret = PMIx_Get(&pname, PMIX_LAUNCH_DIRECTIVES, &info, 1, &val);
    PMIX_INFO_DESTRUCT(&info);
    if (PMIX_SUCCESS == ret) {
        iptr = (pmix_info_t *) val->data.darray->array;
        ninfo = val->data.darray->size;
        for (n = 0; n < ninfo; n++) {
            PMIX_INFO_LIST_XFER(ret, jinfo, &iptr[n]);
        }
        PMIX_VALUE_RELEASE(val);
    }

    /* pass the personality */
    PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_PERSONALITY, schizo->name, PMIX_STRING);

    /* get display options */
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "display", 0, 0))) {
        targv = prte_argv_split(pval->value.data.string, ',');

        for (int idx = 0; idx < prte_argv_count(targv); idx++) {
            if (0 == strcmp(targv[idx], "allocation")) {
                PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_MAPBY, ":DISPLAYALLOC", PMIX_STRING);
            }
            if (0 == strcmp(targv[idx], "map")) {
                PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_MAPBY, ":DISPLAY", PMIX_STRING);
            }
            if (0 == strcmp(targv[idx], "bind")) {
                PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_BINDTO, ":REPORT", PMIX_STRING);
            }
            if (0 == strcmp(targv[idx], "map-devel")) {
                PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_MAPBY, ":DISPLAYDEVEL", PMIX_STRING);
            }
            if (0 == strcmp(targv[idx], "topo")) {
                PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_MAPBY, ":DISPLAYTOPO", PMIX_STRING);
            }
        }
        prte_argv_free(targv);
    }

    /* cannot have both files and directory set for output */
    ptr = NULL;
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "output", 0, 0))) {
        targv = prte_argv_split(pval->value.data.string, ',');

        for (int idx = 0; idx < prte_argv_count(targv); idx++) {
            if (0 == strcmp(targv[idx], "tag")) {
                PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_TAG_OUTPUT, &flag, PMIX_BOOL);
            }
            if (0 == strcmp(targv[idx], "timestamp")) {
                PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_TIMESTAMP_OUTPUT, &flag, PMIX_BOOL);
            }
            if (0 == strcmp(targv[idx], "xml")) {
                PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_MAPBY, ":XMLOUTPUT", PMIX_STRING);
            }
            if (0 == strcmp(targv[idx], "merge-stderr-to-stdout")) {
                PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_MERGE_STDERR_STDOUT, &flag, PMIX_BOOL);
            }
            if (NULL != (ptr = strchr(targv[idx], ':'))) {
                ++ptr;
                ptr = strdup(ptr);
            }
        }
        prte_argv_free(targv);
    }

    param = NULL;
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "output-filename", 0, 0))) {
        param = pval->value.data.string;
    }
    if (NULL != param && NULL != ptr) {
        prte_show_help("help-prted.txt", "both-file-and-dir-set", true, param, ptr);
        return PRTE_ERR_FATAL;
    } else if (NULL != param) {
        /* if we were asked to output to files, pass it along. */
        /* if the given filename isn't an absolute path, then
         * convert it to one so the name will be relative to
         * the directory where prun was given as that is what
         * the user will have seen */
        if (!prte_path_is_absolute(param)) {
            char cwd[PRTE_PATH_MAX];
            if (NULL == getcwd(cwd, sizeof(cwd))) {
                PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
                goto DONE;
            }
            ptr = prte_os_path(false, cwd, param, NULL);
        } else {
            ptr = strdup(param);
        }
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_OUTPUT_TO_FILE, ptr, PMIX_STRING);
        free(ptr);
    } else if (NULL != ptr) {
        /* if we were asked to output to a directory, pass it along. */
        /* If the given filename isn't an absolute path, then
         * convert it to one so the name will be relative to
         * the directory where prun was given as that is what
         * the user will have seen */
        if (!prte_path_is_absolute(ptr)) {
            char cwd[PRTE_PATH_MAX];
            if (NULL == getcwd(cwd, sizeof(cwd))) {
                PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
                goto DONE;
            }
            param = prte_os_path(false, cwd, ptr, NULL);
        } else {
            param = strdup(ptr);
        }
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_OUTPUT_TO_DIRECTORY, param, PMIX_STRING);
        free(param);
    }

    /* check what user wants us to do with stdin */
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "stdin", 0, 0))) {
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_STDIN_TGT, pval->value.data.string, PMIX_STRING);
    }

    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "map-by", 0, 0))) {
        if (donotlaunch && NULL == strcasestr(pval->value.data.string, "donotlaunch")) {
            /* must add directive */
            char *tval;
            prte_asprintf(&tval, "%s:DONOTLAUNCH", pval->value.data.string);
            PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_MAPBY, tval, PMIX_STRING);
            free(tval);
        } else {
            PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_MAPBY, pval->value.data.string, PMIX_STRING);
        }
    }

    /* if the user specified a ranking policy, then set it */
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "rank-by", 0, 0))) {
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_RANKBY, pval->value.data.string, PMIX_STRING);
    }

    /* if the user specified a binding policy, then set it */
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "bind-to", 0, 0))) {
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_BINDTO, pval->value.data.string, PMIX_STRING);
    }
    /* mark if recovery was enabled on the cmd line */
    if (prte_cmd_line_is_taken(prte_cmd_line, "enable-recovery")) {
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_JOB_RECOVERABLE, NULL, PMIX_BOOL);
    }
    /* record the max restarts */
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "max-restarts", 0, 0))
        && 0 < pval->value.data.integer) {
        ui32 = pval->value.data.integer;
        PRTE_LIST_FOREACH(app, &apps, prte_pmix_app_t)
        {
            PMIX_INFO_LIST_ADD(ret, app->info, PMIX_MAX_RESTARTS, &ui32, PMIX_UINT32);
        }
    }
    /* if continuous operation was specified */
    if (prte_cmd_line_is_taken(prte_cmd_line, "continuous")) {
        /* mark this job as continuously operating */
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_JOB_CONTINUOUS, NULL, PMIX_BOOL);
    }

    /* if stop-on-exec was specified */
    if (prte_cmd_line_is_taken(prte_cmd_line, "stop-on-exec")) {
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_DEBUG_STOP_ON_EXEC, NULL, PMIX_BOOL);
    }

    /* check for a job timeout specification, to be provided in seconds
     * as that is what MPICH used
     */
    param = NULL;
    if (NULL != (pval = prte_cmd_line_get_param(prte_cmd_line, "timeout", 0, 0))
        || NULL != (param = getenv("MPIEXEC_TIMEOUT"))) {
        if (NULL != param) {
            i = strtol(param, NULL, 10);
            /* both cannot be present, or they must agree */
            if (NULL != pval && i != pval->value.data.integer) {
                prte_show_help("help-prun.txt", "prun:timeoutconflict", false, prte_tool_basename,
                               pval->value.data.integer, param);
                PRTE_UPDATE_EXIT_STATUS(1);
                goto DONE;
            }
        } else {
            i = pval->value.data.integer;
        }
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_TIMEOUT, &i, PMIX_INT);
    }
    if (prte_cmd_line_is_taken(prte_cmd_line, "get-stack-traces")) {
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_TIMEOUT_STACKTRACES, NULL, PMIX_BOOL);
    }
    if (prte_cmd_line_is_taken(prte_cmd_line, "report-state-on-timeout")) {
        PMIX_INFO_LIST_ADD(ret, jinfo, PMIX_TIMEOUT_REPORT_STATE, NULL, PMIX_BOOL);
    }

    /* give the schizo components a chance to add to the job info */
    prte_schizo.job_info(prte_cmd_line, jinfo);

    /* pickup any relevant envars */
    ninfo = 4;
    PMIX_INFO_CREATE(iptr, ninfo);
    flag = true;
    PMIX_INFO_LOAD(&iptr[0], PMIX_SETUP_APP_ENVARS, &flag, PMIX_BOOL);
    ui32 = geteuid();
    PMIX_INFO_LOAD(&iptr[1], PMIX_USERID, &ui32, PMIX_UINT32);
    ui32 = getegid();
    PMIX_INFO_LOAD(&iptr[2], PMIX_GRPID, &ui32, PMIX_UINT32);
    if (0 == strcasecmp(schizo->name, "prte")) {
        param = strdup("prte");
    } else {
        prte_asprintf(&param, "%s,prte", schizo->name);
    }
    PMIX_INFO_LOAD(&iptr[3], PMIX_PERSONALITY, param, PMIX_STRING);
    free(param);

    PRTE_PMIX_CONSTRUCT_LOCK(&mylock.lock);
    ret = PMIx_server_setup_application(prte_process_info.myproc.nspace, iptr, ninfo, setupcbfunc,
                                        &mylock);
    if (PMIX_SUCCESS != ret) {
        prte_output(0, "Error setting up application: %s", PMIx_Error_string(ret));
        PRTE_PMIX_DESTRUCT_LOCK(&mylock.lock);
        PRTE_UPDATE_EXIT_STATUS(ret);
        goto DONE;
    }
    PRTE_PMIX_WAIT_THREAD(&mylock.lock);
    PMIX_INFO_FREE(iptr, ninfo);
    if (PMIX_SUCCESS != mylock.status) {
        prte_output(0, "Error setting up application: %s", PMIx_Error_string(mylock.status));
        PRTE_UPDATE_EXIT_STATUS(mylock.status);
        PRTE_PMIX_DESTRUCT_LOCK(&mylock.lock);
        goto DONE;
    }
    PRTE_PMIX_DESTRUCT_LOCK(&mylock.lock);
    /* transfer any returned ENVARS to the job_info */
    if (NULL != mylock.info) {
        for (n = 0; n < mylock.ninfo; n++) {
            if (0 == strncmp(mylock.info[n].key, PMIX_SET_ENVAR, PMIX_MAX_KEYLEN)
                || 0 == strncmp(mylock.info[n].key, PMIX_ADD_ENVAR, PMIX_MAX_KEYLEN)
                || 0 == strncmp(mylock.info[n].key, PMIX_UNSET_ENVAR, PMIX_MAX_KEYLEN)
                || 0 == strncmp(mylock.info[n].key, PMIX_PREPEND_ENVAR, PMIX_MAX_KEYLEN)
                || 0 == strncmp(mylock.info[n].key, PMIX_APPEND_ENVAR, PMIX_MAX_KEYLEN)) {
                PMIX_INFO_LIST_XFER(ret, jinfo, &mylock.info[n]);
            }
        }
        PMIX_INFO_FREE(mylock.info, mylock.ninfo);
    }

    /* convert the job info into an array */
    PMIX_INFO_LIST_CONVERT(ret, jinfo, &darray);
    if (PMIX_ERR_EMPTY == ret) {
        iptr = NULL;
        ninfo = 0;
    } else if (PMIX_SUCCESS != ret) {
        PMIX_ERROR_LOG(ret);
        PRTE_UPDATE_EXIT_STATUS(rc);
        goto DONE;
    } else {
        iptr = (pmix_info_t *) darray.array;
        ninfo = darray.size;
    }
    PMIX_INFO_LIST_RELEASE(jinfo);

    /* convert the apps to an array */
    napps = prte_list_get_size(&apps);
    PMIX_APP_CREATE(papps, napps);
    n = 0;
    PRTE_LIST_FOREACH(app, &apps, prte_pmix_app_t)
    {
        papps[n].cmd = strdup(app->app.cmd);
        papps[n].argv = prte_argv_copy(app->app.argv);
        papps[n].env = prte_argv_copy(app->app.env);
        papps[n].cwd = strdup(app->app.cwd);
        papps[n].maxprocs = app->app.maxprocs;
        PMIX_INFO_LIST_CONVERT(ret, app->info, &darray);
        if (PMIX_SUCCESS != ret) {
            if (PMIX_ERR_EMPTY == ret) {
                papps[n].info = NULL;
                papps[n].ninfo = 0;
            } else {
                PMIX_ERROR_LOG(ret);
                PRTE_UPDATE_EXIT_STATUS(rc);
                goto DONE;
            }
        } else {
            papps[n].info = (pmix_info_t *) darray.array;
            papps[n].ninfo = darray.size;
        }
        /* pickup any relevant envars */
        rc = prte_schizo.parse_env(prte_cmd_line, environ, &papps[n].env, false);
        if (PRTE_SUCCESS != rc) {
            PRTE_UPDATE_EXIT_STATUS(rc);
            goto DONE;
        }
        ++n;
    }

    if (verbose) {
        prte_output(0, "Spawning job");
    }

    PRTE_PMIX_CONSTRUCT_LOCK(&lock);
    ret = pmix_server_spawn_fn(&parent, iptr, ninfo, papps, napps, spcbfunc, &lock);
    if (PRTE_SUCCESS != ret) {
        prte_output(0, "PMIx_Spawn failed (%d): %s", ret, PMIx_Error_string(ret));
        rc = ret;
        PRTE_UPDATE_EXIT_STATUS(rc);
        goto DONE;
    }
    /* we have to cycle the event library here so we can process
     * the spawn request */
    while (prte_event_base_active && lock.active) {
        prte_event_loop(prte_event_base, PRTE_EVLOOP_ONCE);
    }
    PRTE_ACQUIRE_OBJECT(&lock.lock);
    if (PMIX_SUCCESS != lock.status) {
        PRTE_UPDATE_EXIT_STATUS(lock.status);
        goto DONE;
    }
    PMIX_LOAD_NSPACE(spawnednspace, lock.msg);
    PRTE_PMIX_DESTRUCT_LOCK(&lock);
    if (verbose) {
        prte_output(0, "JOB %s EXECUTING", PRTE_JOBID_PRINT(spawnednspace));
    }
    /* need to "pull" the IOF from the spawned job since we didn't
     * go thru PMIx_Spawn to start it - and thus, PMIx didn't
     * "pull" it for us */
    PMIX_LOAD_PROCID(&pname, spawnednspace, PMIX_RANK_WILDCARD);
    ret = PMIx_IOF_pull(&pname, 1, NULL, 0,
                        PMIX_FWD_STDOUT_CHANNEL | PMIX_FWD_STDERR_CHANNEL
                            | PMIX_FWD_STDDIAG_CHANNEL,
                        NULL, NULL, NULL);
    if (PMIX_SUCCESS != ret && PMIX_OPERATION_SUCCEEDED != ret) {
        prte_output(0, "IOF pull failed: %s", PMIx_Error_string(ret));
    }

    /* push our stdin to the apps */
    PMIX_LOAD_PROCID(&pname, spawnednspace, 0); // forward stdin to rank=0
    PMIX_INFO_CREATE(iptr, 1);
    PMIX_INFO_LOAD(&iptr[0], PMIX_IOF_PUSH_STDIN, NULL, PMIX_BOOL);
    PRTE_PMIX_CONSTRUCT_LOCK(&lock);
    ret = PMIx_IOF_push(&pname, 1, NULL, iptr, 1, opcbfunc, &lock);
    if (PMIX_SUCCESS != ret && PMIX_OPERATION_SUCCEEDED != ret) {
        prte_output(0, "IOF push of stdin failed: %s", PMIx_Error_string(ret));
    } else if (PMIX_SUCCESS == ret) {
        PRTE_PMIX_WAIT_THREAD(&lock);
    }
    PRTE_PMIX_DESTRUCT_LOCK(&lock);
    PMIX_INFO_FREE(iptr, 1);

    pmix_status_t rc_define=PMIX_RC_DEFINE;
    /* register the resource change cmd handler */
    PMIx_Register_event_handler(&rc_define, 1, NULL, 0, rchandler, NULL, NULL);

    /* create a pset for the job */

    pmix_data_buffer_t *buf;
    

    prte_daemon_cmd_flag_t cmd = PRTE_DYNRES_DEFINE_PSET;

    char *pset_name = strdup("test1");
    size_t nprocs = prte_get_job_data_object(spawnednspace)->num_procs;
    prte_pointer_array_t *pset_procs_parray = prte_get_job_data_object(spawnednspace)->procs;
     


    int ndaemons = prte_process_info.num_daemons;
    pmix_proc_t daemon_procid;
    PMIX_LOAD_PROCID(&daemon_procid, PRTE_PROC_MY_HNP->nspace, 0);
    for(i = 0; i < ndaemons; i++){
        PMIX_DATA_BUFFER_CREATE(buf);
        ret = PMIx_Data_pack(NULL, buf, &cmd, 1, PMIX_UINT8);
        ret = PMIx_Data_pack(NULL, buf, &nprocs, 1, PMIX_SIZE);
    
        ret = PMIx_Data_pack(NULL, buf, &pset_name, 1, PMIX_STRING);
    
        for(n = 0; n < nprocs; n++){
            pmix_proc_t pset_proc;
            prte_proc_t *prte_proc = (prte_proc_t *) pset_procs_parray->addr[n];
            PMIX_PROC_LOAD(&pset_proc, prte_proc->name.nspace, prte_proc->name.rank);
            ret = PMIx_Data_pack(NULL, buf, &pset_proc, 1, PMIX_PROC);
        }

        daemon_procid.rank = i;
        prte_rml.send_buffer_nb(&daemon_procid, buf, PRTE_RML_TAG_MALLEABILITY, prte_rml_send_callback, NULL);
    }
    free(pset_name);
    //printf("\nPRRTE HNP Server pid:\n %lu\n\n", (unsigned long) getpid());
    master_timing_list = (node_t *)calloc(1, sizeof(node_t));
    
    timings_my_rank = 0;
proceed:
    /* loop the event lib until an exit event is detected */
    while (prte_event_base_active) {
        prte_event_loop(prte_event_base, PRTE_EVLOOP_ONCE);
    }
    PRTE_ACQUIRE_OBJECT(prte_event_base_active);

    /* close the push of our stdin */
    PMIX_INFO_LOAD(&info, PMIX_IOF_COMPLETE, NULL, PMIX_BOOL);
    PRTE_PMIX_CONSTRUCT_LOCK(&lock);
    ret = PMIx_IOF_push(NULL, 0, NULL, &info, 1, opcbfunc, &lock);
    if (PMIX_SUCCESS != ret && PMIX_OPERATION_SUCCEEDED != ret) {
        prte_output(0, "IOF close of stdin failed: %s", PMIx_Error_string(ret));
    } else if (PMIX_SUCCESS == ret) {
        PRTE_PMIX_WAIT_THREAD(&lock);
    }
    PRTE_PMIX_DESTRUCT_LOCK(&lock);
    PMIX_INFO_DESTRUCT(&info);

DONE:
    /* cleanup and leave */
    prte_finalize();

    if (NULL != mypidfile) {
        unlink(mypidfile);
    }

    if (prte_debug_flag) {
        fprintf(stderr, "exiting with status %d\n", prte_exit_status);
    }
    exit(prte_exit_status);
}

static void clean_abort(int fd, short flags, void *arg)
{
    /* if we have already ordered this once, don't keep
     * doing it to avoid race conditions
     */
    if (prte_mutex_trylock(&prun_abort_inprogress_lock)) { /* returns 1 if already locked */
        if (forcibly_die) {
            /* exit with a non-zero status */
            exit(1);
        }
        fprintf(stderr,
                "%s: abort is already in progress...hit ctrl-c again to forcibly terminate\n\n",
                prte_tool_basename);
        forcibly_die = true;
        /* reset the event */
        prte_event_add(&term_handler, NULL);
        return;
    }

    fflush(stderr);
    /* ensure we exit with a non-zero status */
    PRTE_UPDATE_EXIT_STATUS(PRTE_ERROR_DEFAULT_EXIT_CODE);
    /* ensure that the forwarding of stdin stops */
    prte_job_term_ordered = true;
    /* tell us to be quiet - hey, the user killed us with a ctrl-c,
     * so need to tell them that!
     */
    prte_execute_quiet = true;
    /* We are in an event handler; the job completed procedure
     will delete the signal handler that is currently running
     (which is a Bad Thing), so we can't call it directly.
     Instead, we have to exit this handler and setup to call
     job_completed() after this. */
    prte_plm.terminate_orteds();
}

static struct timeval current, last = {0, 0};
static bool first = true;

/*
 * Attempt to terminate the job and wait for callback indicating
 * the job has been aborted.
 */
static void abort_signal_callback(int fd)
{
    uint8_t foo = 1;
    char *msg
        = "Abort is in progress...hit ctrl-c again within 5 seconds to forcibly terminate\n\n";

    /* if this is the first time thru, just get
     * the current time
     */
    if (first) {
        first = false;
        gettimeofday(&current, NULL);
    } else {
        /* get the current time */
        gettimeofday(&current, NULL);
        /* if this is within 5 seconds of the
         * last time we were called, then just
         * exit - we are probably stuck
         */
        if ((current.tv_sec - last.tv_sec) < 5) {
            exit(1);
        }
        if (-1 == write(1, (void *) msg, strlen(msg))) {
            exit(1);
        }
    }
    /* save the time */
    last.tv_sec = current.tv_sec;
    /* tell the event lib to attempt to abnormally terminate */
    if (-1 == write(term_pipe[1], &foo, 1)) {
        exit(1);
    }
}

static int prep_singleton(const char *name)
{
    char *ptr, *p1;
    prte_job_t *jdata;
    prte_node_t *node;
    prte_proc_t *proc;
    int rc;
    pmix_rank_t rank;
    prte_app_context_t *app;
    char cwd[PRTE_PATH_MAX];

    ptr = strdup(name);
    p1 = strrchr(ptr, '.');
    *p1 = '\0';
    ++p1;
    rank = strtoul(p1, NULL, 10);
    jdata = PRTE_NEW(prte_job_t);
    PMIX_LOAD_NSPACE(jdata->nspace, ptr);
    free(ptr);
    rc = prte_set_job_data_object(jdata);
    if (PRTE_SUCCESS != rc) {
        PRTE_UPDATE_EXIT_STATUS(PRTE_ERR_FATAL);
        PRTE_RELEASE(jdata);
        return PRTE_ERR_FATAL;
    }
    /* must have an app */
    app = PRTE_NEW(prte_app_context_t);
    app->app = strdup(jdata->nspace);
    app->num_procs = 1;
    prte_argv_append_nosize(&app->argv, app->app);
    getcwd(cwd, sizeof(cwd));
    app->cwd = strdup(cwd);
    prte_pointer_array_set_item(jdata->apps, 0, app);
    jdata->num_apps = 1;

    /* add a map */
    jdata->map = PRTE_NEW(prte_job_map_t);
    /* add our node to the map since the singleton must
     * be here */
    node = (prte_node_t *) prte_pointer_array_get_item(prte_node_pool, PRTE_PROC_MY_NAME->rank);
    PRTE_RETAIN(node);
    prte_pointer_array_add(jdata->map->nodes, node);
    ++(jdata->map->num_nodes);

    /* create a proc for the singleton */
    proc = PRTE_NEW(prte_proc_t);
    PMIX_LOAD_PROCID(&proc->name, jdata->nspace, rank);
    proc->rank = proc->name.rank;
    proc->parent = PRTE_PROC_MY_NAME->rank;
    proc->app_idx = 0;
    proc->app_rank = rank;
    proc->local_rank = 0;
    proc->node_rank = 0;
    proc->state = PRTE_PROC_STATE_RUNNING;
    /* link it to the job */
    PRTE_RETAIN(jdata);
    proc->job = jdata;
    /* link it to the app */
    PRTE_RETAIN(proc);
    prte_pointer_array_set_item(&app->procs, rank, proc);
    app->first_rank = rank;
    /* link it to the node */
    PRTE_RETAIN(node);
    proc->node = node;
    /* add it to the job */
    prte_pointer_array_set_item(jdata->procs, rank, proc);
    jdata->num_procs = 1;
    jdata->num_local_procs = 1;
    /* add it to the node */
    PRTE_RETAIN(proc);
    prte_pointer_array_add(node->procs, proc);
    node->num_procs = 1;
    node->slots_inuse = 1;

    return PRTE_SUCCESS;
}

static void signal_forward_callback(int signum, short args, void *cbdata)
{
    pmix_status_t rc;
    pmix_proc_t proc;
    pmix_info_t info;

    if (verbose) {
        fprintf(stderr, "%s: Forwarding signal %d to job\n", prte_tool_basename, signum);
    }

    /* send the signal out to the processes */
    PMIX_LOAD_PROCID(&proc, spawnednspace, PMIX_RANK_WILDCARD);
    PMIX_INFO_LOAD(&info, PMIX_JOB_CTRL_SIGNAL, &signum, PMIX_INT);
    rc = PMIx_Job_control(&proc, 1, &info, 1, NULL, NULL);
    if (PMIX_SUCCESS != rc && PMIX_OPERATION_SUCCEEDED != rc) {
        fprintf(stderr, "Signal %d could not be sent to job %s (returned %s)", signum,
                spawnednspace, PMIx_Error_string(rc));
    }
}

/**
 * Deal with sigpipe errors
 */
static int sigpipe_error_count = 0;
static void epipe_signal_callback(int fd, short args, void *cbdata)
{
    sigpipe_error_count++;

    if (10 < sigpipe_error_count) {
        /* time to abort */
        prte_output(0, "%s: SIGPIPE detected - aborting", prte_tool_basename);
        clean_abort(0, 0, NULL);
    }

    return;
}



static pmix_status_t parse_rc_cmd(char *cmd, pmix_res_change_type_t *_rc_type, size_t *nprocs){

    char * token= strtok(cmd, " ");
    if(token==NULL || 0!=strncmp(token, "pmix_session", 12)){
        return PMIX_ERR_BAD_PARAM;
    }

    token= strtok(NULL, " ");
    if(token == NULL)return PMIX_ERR_BAD_PARAM;
    if(0 == strncmp(token, "add", 3))*_rc_type = PMIX_RES_CHANGE_ADD;
    else if(0 == strncmp(token, "sub", 3))*_rc_type = PMIX_RES_CHANGE_SUB;
    else return PMIX_ERR_BAD_PARAM;

    token = strtok(NULL, " ");
    return ((*nprocs = atoi(token)) <= 0) ? PMIX_ERR_BAD_PARAM : PMIX_SUCCESS;



}

static void setup_resource_add(prte_job_t *job_data, size_t rc_nprocs, prte_proc_t ***_delta_procs){
    
    size_t n;
    prte_proc_t *proc;
    
    /* determine highest rank in this job. 
     * For simplicity we use a stack approach for resources to avoid fragmentation of the rank space.
     * Until a consens dealing with the consistency of job updates on the client side,
     * this also avoids dealing with recycled process ids on the client side 
     * Moreover the prte proc state machine requires unique ids*/ 
        
    pmix_rank_t num_procs = job_data->num_procs;
    pmix_rank_t highest_rank = 0;
    for(n = 0; n < job_data->procs->size; n++){
        if(NULL == (proc = prte_pointer_array_get_item(job_data->procs, n))){
            continue;
        }
        if(proc->rank > highest_rank){
            highest_rank = proc->rank;
        }
    }
    if(highest_rank_global > highest_rank){
        highest_rank = highest_rank_global;
    }
    //printf("PRRTE: highest rank: %d\n", highest_rank);


    /* create the delta pset starting at highest rank */
    prte_proc_t **delta_procs = malloc(rc_nprocs * sizeof(prte_proc_t*));
    *_delta_procs = delta_procs;
    for(n = 0; n < rc_nprocs; n++){
        delta_procs[n] = PRTE_NEW(prte_proc_t);
        PMIX_LOAD_NSPACE(delta_procs[n]->name.nspace, spawnednspace);
        delta_procs[n]->name.rank = delta_procs[n]->rank = delta_procs[n]->app_rank =  highest_rank + 1 + n;
        delta_procs[n]->app_idx = 0;
        delta_procs[n]->state = PRTE_PROC_STATE_INIT;
    }

    /* set the deamon vpid, node, ranks etc. for every proc */
    prte_node_t **job_nodes = job_data->map->nodes->addr; 
    size_t proc_index = 0, node_index, daemon_index;

    while(proc_index < rc_nprocs){
        prte_node_t *node;

        /* first fill up nodes allocated to the job */
        for(node_index = 0; node_index < job_data->map->nodes->size; node_index++){
            if(NULL == (node = job_nodes[node_index])){
                continue;
            }

            /* Get the daemon of this node to set the parent of the proc*/
            prte_proc_t *daemon_proc = NULL;
            pmix_rank_t parent_vpid = PMIX_RANK_INVALID;
            for(daemon_index = 0; daemon_index < prte_get_job_data_object(PRTE_PROC_MY_PROCID->nspace)->procs->size; daemon_index++){
                daemon_proc = prte_pointer_array_get_item(prte_get_job_data_object(PRTE_PROC_MY_PROCID->nspace)->procs, daemon_index);
                if(NULL != daemon_proc && (0 == strcmp(daemon_proc->node->name,node->name))){
                    parent_vpid = daemon_proc->name.rank;
                }
            }

            int32_t cur_slot = node->slots_inuse;
            for(; cur_slot < node->slots && proc_index < rc_nprocs; ){
                /* set parent */
                delta_procs[proc_index]->parent = parent_vpid;
               /* set node */
                PRTE_RETAIN(node);
                delta_procs[proc_index]->node = node;
                delta_procs[proc_index]->node_rank = cur_slot;
                delta_procs[proc_index]->local_rank = cur_slot;
                printf("MASTER: ADD proc on node %d, with local rank %d, node slots = %d\n", node_index, cur_slot, node->slots);
                /* add proc to node */
                PRTE_RETAIN(delta_procs[proc_index]);
                prte_pointer_array_add(node->procs, delta_procs[proc_index]);
                node->num_procs++;
                node->slots_inuse++;
                node->slots_available--;
                /* and connect it back to its job object, if not already done */
                if (NULL == delta_procs[proc_index]->job) {
                    PRTE_RETAIN(job_data);
                    delta_procs[proc_index]->job = job_data;
                }
                cur_slot++;
                proc_index++;
            }
            if(proc_index == rc_nprocs){
                break;
            }
        }

        
        /* If we reach this, there were not enough nodes to add all processes. 
         * So we try to add another node from the daemon job 
         */
        if(proc_index < rc_nprocs){
            bool node_added = false, already_allocated;
            prte_job_t *djob = prte_get_job_data_object(PRTE_PROC_MY_NAME->nspace);
            prte_node_t *dnode;

            /* Find a node from the DVM that is not yet assigned to the job */
            for(n = 0; n < djob->map->nodes->size; n++){
                if(NULL == (dnode = prte_pointer_array_get_item(djob->map->nodes, n))){
                    continue;
                }
                already_allocated = false;
                for(node_index = 0; node_index < job_data->map->nodes->size; node_index++){
                    if(NULL == (node = job_nodes[node_index])){
                        continue;
                    }
                    if(0 == strcmp(dnode->name, node->name)){
                        already_allocated = true;
                        break;
                    }
                }
                /* add the node to the job */
                if(!already_allocated){
                    //printf("NOT ENOUGH NODES: Adding new node %s to job %s to fullfill the request\n", dnode->name, job_data->nspace);
                    PRTE_RETAIN(dnode);
                    prte_pointer_array_add(job_data->map->nodes, dnode);
                    job_data->total_slots_alloc += dnode->slots_available;
                    job_data->map->num_nodes++;
                    job_data->num_daemons_reported++;
                    node_added = true;
                    break;
                }

            }
            /* We couldn't add another node to satisfy the request so leave now */
            if(!node_added){
                break;
            }
        }


    }

    /* TODO: cleanly exit res_change handler */
    if(proc_index < rc_nprocs){
        printf("Not enough nodes/slots available for this request.\n");
    }
    /* add the procs to the job and app context */
    prte_app_context_t *app = job_data->apps->addr[0];
    for(n = 0; n < rc_nprocs; n++){
        prte_pointer_array_add(job_data->procs, delta_procs[n]);
        PRTE_RETAIN(delta_procs[n]);
        prte_pointer_array_add(&app->procs, delta_procs[n]);
        job_data->num_procs++;
        app->num_procs++;
    }
    bool fully_described = true;
    prte_set_attribute(&job_data->attributes, PRTE_JOB_FULLY_DESCRIBED, PRTE_ATTR_GLOBAL, &fully_described, PMIX_BOOL);
    prte_set_attribute(&job_data->attributes, PRTE_JOB_FIXED_DVM, PRTE_ATTR_GLOBAL, &fully_described, PMIX_BOOL);

    highest_rank_global = highest_rank + rc_nprocs;
}

static void setup_resource_sub(pmix_nspace_t job, prte_job_t *job_data, size_t rc_nprocs, prte_proc_t ***_delta_procs){
    
    size_t n, i, rm_nodes = 0;
    prte_job_t *job_data_cpy = prte_get_job_data_object(job);
    prte_attribute_t *attr;

    PRTE_RELEASE(job_data->apps);
    PRTE_RELEASE(job_data->procs);
    //PRTE_DESTRUCT(&job_data->attributes);
    /* create a copy of our job data object 
     * We need to use a copy as we do not want to change our stored job data,
     * instead we want to send an adjusted job data object only to update the pmix sever's job data
     */
    memcpy(job_data, job_data_cpy, sizeof(prte_job_t));

    /*TODO list foreach */
    /* copy the attributes list */
    memset(&job_data->attributes, 0, sizeof(prte_list_t));
    PRTE_CONSTRUCT(&job_data->attributes, prte_list_t);

    /* set these attributes, just in case they aren't set yet */
    bool fully_described = true;
    prte_set_attribute(&job_data->attributes, PRTE_JOB_FULLY_DESCRIBED, PRTE_ATTR_GLOBAL, &fully_described, PMIX_BOOL);
    prte_set_attribute(&job_data->attributes, PRTE_JOB_FIXED_DVM, PRTE_ATTR_GLOBAL, &fully_described, PMIX_BOOL);
    prte_set_attribute(&job_data->attributes, PRTE_JOB_LAUNCH_PROXY, PRTE_ATTR_GLOBAL, &prte_process_info.myproc, PMIX_PROC);

    //size_t list_length = job_data_cpy->attributes.prte_list_length;
    //int ctr=0;
    //PRTE_LIST_FOREACH(attr, &job_data_cpy->attributes, prte_attribute_t){
    //        prte_add_attribute(&new_attributes_list, attr->key, attr->local, &attr->data.data, attr->data.type);
    //        printf("prte attribute: %d\n", attr->key);
    //}


    /* copy app contexts (as we adjust the num procs) */
    job_data->apps = PRTE_NEW(prte_pointer_array_t);
    prte_pointer_array_init(job_data->apps, job_data_cpy->apps->size, PRTE_GLOBAL_ARRAY_MAX_SIZE, 2);
    for(n = 0; n < job_data_cpy->apps->size; n++){
        prte_app_context_t *app_ptr = prte_pointer_array_get_item(job_data_cpy->apps, n);
        if(NULL != app_ptr){
            prte_app_context_t *app_cpy = PRTE_NEW(prte_app_context_t);
            memcpy(app_cpy, app_ptr, sizeof(prte_app_context_t));
            for(i = 0; i< app_cpy->num_procs; i++){
                prte_proc_t *app_proc;
                // Protect the procs so we can release the job data later
                if(NULL != (app_proc = prte_pointer_array_get_item(&app_cpy->procs, i))){
                    PRTE_RETAIN(app_proc);
                }
            }
            prte_pointer_array_add(job_data->apps, app_cpy);
        }
    }
    /* copy the procs array so we can adjust it without changing the actual job data */
    job_data->procs = PRTE_NEW(prte_pointer_array_t);
    prte_pointer_array_init(job_data->procs, 4, PRTE_GLOBAL_ARRAY_MAX_SIZE,
                            PRTE_GLOBAL_ARRAY_BLOCK_SIZE);
    for(n = 0; n < job_data_cpy->procs->size; n++){
        prte_proc_t *proc_ptr = prte_pointer_array_get_item(job_data_cpy->procs, n);
        if(NULL != proc_ptr){
            int ret = prte_pointer_array_add(job_data->procs, proc_ptr);
            /* Protect the procs when we release the array later when releasing the job data */
            PRTE_RETAIN(proc_ptr); 
        }
    }
    prte_node_t *node;
    /* copy the map and create copies of the nodes in the list*/
    job_data->map = PRTE_NEW(prte_job_map_t);
    memcpy(job_data->map, job_data_cpy->map, sizeof(prte_job_map_t));
    job_data->map->nodes = PRTE_NEW(prte_pointer_array_t);
    for(i = 0; i < job_data_cpy->map->nodes->size; i++){
        if(NULL == (node = prte_pointer_array_get_item(job_data_cpy->map->nodes, i))){
            continue;
        }
        prte_node_t *node_cpy = PRTE_NEW(prte_node_t);
        memcpy(node_cpy, node, sizeof(prte_node_t));
        node_cpy->procs = PRTE_NEW(prte_pointer_array_t);
        for(int k = 0; k < node->procs->size; k++){
            prte_proc_t *n_proc;
            if(NULL == (n_proc = prte_pointer_array_get_item(node->procs, k))){
                continue;
            }
            prte_pointer_array_add(node_cpy->procs, n_proc);
        }

        prte_pointer_array_add(job_data->map->nodes, node_cpy);
    }

    PRTE_RETAIN(job_data->map->nodes);

    /* create the delta pset, i.e. determine the procs to be finalized */
    prte_proc_t **delta_procs = malloc(rc_nprocs * sizeof(prte_proc_t*));
    *_delta_procs = delta_procs;

    /* We traverse the list of nodes and their procs in reverse order and choose the ones of our job */
    prte_app_context_t *app = job_data->apps->addr[0];
    prte_node_t **job_nodes = job_data->map->nodes->addr; // or should we use dvm nodes?
    
    size_t proc_index = 0, node_index;
    for(node_index = job_data->map->num_nodes -1; node_index >= 0; node_index--){


        if(NULL == (node = job_nodes[node_index])){
            continue;
        }

        int32_t cur_slot = node->slots_inuse - 1;
        printf("SUB: node index: %d, slots_in_use: %d\n", node_index, node->slots_inuse);
        for(; cur_slot >= 0 && proc_index < rc_nprocs; ){
            printf("cur_slot: %d\n", cur_slot);
            for(i = node->procs->size; i >= 0; i--){
                if(NULL == (delta_procs[proc_index] = prte_pointer_array_get_item(node->procs, i))){
                    continue;
                }

                if(delta_procs[proc_index]->name.rank > highest_rank_global){
                    highest_rank_global = delta_procs[proc_index]->name.rank;
                }

                /* found a proc to remove */
                if( PMIX_CHECK_NSPACE(delta_procs[proc_index]->name.nspace, job) && 
                    delta_procs[proc_index]->node_rank == cur_slot){
                        //PRTE_FLAG_TEST
                    break;
                }

                if(i == 0){
                    delta_procs[proc_index] = NULL;
                }

            }

            --cur_slot;

            if(NULL == delta_procs[proc_index]){
                continue;
            }

            /* we added this proc to the delta pset. Now remove it from the job & app proc lists */

            for(n = 0; n < job_data->procs->size; n++){
                prte_proc_t * proct;
                if(NULL == (proct = prte_pointer_array_get_item(job_data->procs, n))){
                    continue;
                }
                if(proct->name.rank == delta_procs[proc_index]->name.rank){
                    int ret = prte_pointer_array_set_item(job_data->procs, n, NULL);
                    if(NULL != prte_pointer_array_get_item(job_data->procs, n));
                    PRTE_RELEASE(proct); // for referenece counting
                }
            }

            for(n = 0; n < app->procs.size; n++){
                prte_proc_t * proct; 
                if(NULL == (proct = prte_pointer_array_get_item(job_data->procs, n))){
                    continue;
                }
                if(proct->name.rank == delta_procs[proc_index]->name.rank){ 

                    prte_pointer_array_set_item(&app->procs, n, NULL);
                    PRTE_RELEASE(proct); // for referenece counting
                }
            }
            /*if(--cur_slot == 0){
                ++rm_nodes;
            }*/
            proc_index++;
        }
        if(proc_index == rc_nprocs){
            break;
        }
    }
    if(proc_index < rc_nprocs){
        printf("Not enough nodes/slots available for this request\n");
    }

    /* set the number of nodes and procs accordingly */
    job_data->num_procs -= rc_nprocs;
    app->num_procs -= rc_nprocs;

    /* TODO: Need to respect procs from other jobs */
    //job_data->map->num_nodes -= rm_nodes;

    
    

    memset(&job_data->children, 0, sizeof(prte_list_t));
    PRTE_CONSTRUCT(&job_data->children, prte_list_t);
    /* prepend the launch message with the sub command, so it is handle correctly */ 
    prte_daemon_cmd_flag_t command = PRTE_DAEMON_DVM_SUB_PROCS;
    PMIx_Data_pack(NULL, &job_data->launch_msg, &command, 1, PMIX_UINT8);
}

static void _rchandler(int sd, short args, void *cbdata)
{
    prte_pmix_server_op_caddy_t *scd = (prte_pmix_server_op_caddy_t *) cbdata;
    pmix_info_t *info = scd->info;
    size_t ninfo = scd->ninfo;
    size_t n, i, ret, sz, rc_nprocs;
    pmix_status_t rc=PMIX_SUCCESS;
    char *recv_cmd = NULL;
    char *delta_pset_name;
    pmix_res_change_type_t rc_type;
    pmix_data_buffer_t *buf;
    prte_grpcomm_signature_t *sig;

    if(0 < prte_list_get_size(&prte_pmix_server_globals.res_changes)){
        scd->evncbfunc(PMIX_ERR_BAD_PARAM, NULL, 0, NULL, NULL, cbdata);
        return;
    }

    /* parse the input */
    for(n=0; n < ninfo; n++){
        if(0 == strcmp(info[n].key, "PMIX_RC_CMD")){
            PMIX_VALUE_UNLOAD(rc, &info[n].value, (void**)&recv_cmd, &sz);
            //prte_output(2, "SERVER: RECEIVED RC_CMD %s\n", recv_cmd);
            
            rc=parse_rc_cmd(recv_cmd, &rc_type, &rc_nprocs);
            if(rc != PMIX_SUCCESS){
                printf("Error parsing rc command\n");
                scd->evncbfunc(PMIX_ERR_BAD_PARAM, NULL, 0, NULL, NULL, cbdata);
                return;
            }
        }
    }
    free(recv_cmd);

    init_add_timing(master_timing_list, &cur_master_timing_frame, sizeof(timing_frame_master_t));
    make_timestamp_base(&cur_master_timing_frame->rc_start);


    /* get the job data object & do a sanity check */
    prte_job_t *cur_job_data = prte_get_job_data_object(spawnednspace);
    if(rc_type == PMIX_RES_CHANGE_SUB && rc_nprocs >= cur_job_data->num_procs){
        scd->evncbfunc(PMIX_ERR_BAD_PARAM, NULL, 0, NULL, NULL, cbdata);
        return;
    }

    make_timestamp_base(&cur_master_timing_frame->jdata_start);
    /* create the delta pset and adjust the job_data */
    prte_job_t *job_data;
    prte_proc_t **delta_procs;
    if(rc_type == PMIX_RES_CHANGE_ADD){
        job_data = prte_get_job_data_object(spawnednspace);
        setup_resource_add(job_data, rc_nprocs, &delta_procs);
    }else{
        job_data = PRTE_NEW(prte_job_t);
        setup_resource_sub(spawnednspace, job_data, rc_nprocs, &delta_procs);
    }
    make_timestamp_base(&cur_master_timing_frame->jdata_end);


    /* we are the master so update the job data now */
    //prte_set_or_replace_job_data_object(job_data);

    //prte_list_item_t *item = prte_list_get_last(&job_data->attributes);
    //item->prte_list_next = &(job_data->attributes.prte_list_sentinel);


    delta_pset_name = (char*) malloc(256);
    sprintf(delta_pset_name, "rc%d", res_change_cnt++);

    set_res_change_id(&cur_master_timing_frame->res_change_id, delta_pset_name);
    cur_master_timing_frame->res_change_type = rc_type;
    cur_master_timing_frame->res_change_size = rc_nprocs;
        
    //printf("Received resource change '%s' of type %s\n Delta Pset:\n", delta_pset_name, rc_type == PMIX_RES_CHANGE_ADD ? "ADD" : "SUB");
    //for(n = 0; n < rc_nprocs; n++){
    //    printf("    [%d: %s]\n", n, PRTE_NAME_PRINT(&delta_procs[n]->name));
    //}
    //printf("--> job size after resource change will be: %d \n\n",job_data->num_procs);

    make_timestamp_base(&cur_master_timing_frame->pset_start);
    /* Send the 'Define Pset' command */
    prte_daemon_cmd_flag_t cmd = PRTE_DYNRES_DEFINE_PSET;
    int ndaemons = prte_process_info.num_daemons;
    pmix_proc_t daemon_procid;
    PMIX_LOAD_PROCID(&daemon_procid, PRTE_PROC_MY_HNP->nspace, 0);
    for(n=0; n < ndaemons; n++){
        PMIX_DATA_BUFFER_CREATE(buf);
        daemon_procid.rank = n;
        ret = PMIx_Data_pack(NULL, buf, &cmd, 1, PMIX_UINT8);
        ret = PMIx_Data_pack(NULL, buf, &rc_nprocs, 1, PMIX_SIZE);
        ret = PMIx_Data_pack(NULL, buf, (void*)&delta_pset_name, 1, PMIX_STRING);
        for(i = 0; i < rc_nprocs; i++){
            pmix_proc_t pset_proc;
            prte_proc_t *prte_proc = (prte_proc_t *) delta_procs[i];
            PMIX_PROC_LOAD(&pset_proc, prte_proc->name.nspace, prte_proc->name.rank);
            ret = PMIx_Data_pack(NULL, buf, &pset_proc, 1, PMIX_PROC);
        }
        prte_rml.send_buffer_nb(&daemon_procid, buf, PRTE_RML_TAG_MALLEABILITY, prte_rml_send_callback, NULL);
    }

    if(rc_type == PMIX_RES_CHANGE_SUB){
        prte_state_caddy_t *cd = PRTE_NEW(prte_state_caddy_t);
        cd->jdata = job_data;
        cd->job_state = (rc_type == PMIX_RES_CHANGE_ADD) ? PRTE_JOB_STATE_LAUNCH_APPS : PRTE_JOB_STATE_SUB;
        prte_plm_base_launch_apps(0,0, cd);
        /* message goes to all daemons */
        sig = PRTE_NEW(prte_grpcomm_signature_t);
        sig->signature = (pmix_proc_t *) malloc(sizeof(pmix_proc_t));
        PMIX_LOAD_PROCID(&sig->signature[0], PRTE_PROC_MY_NAME->nspace, PMIX_RANK_WILDCARD);
        sig->sz = 1;
        if (PRTE_SUCCESS != (rc = prte_grpcomm.xcast(sig, PRTE_RML_TAG_DAEMON, &job_data->launch_msg))) {
            PRTE_ERROR_LOG(rc);
            PRTE_RELEASE(sig);
            return;
        }
        PMIX_DATA_BUFFER_DESTRUCT(&job_data->launch_msg);
        PMIX_DATA_BUFFER_CONSTRUCT(&job_data->launch_msg);
        /* maintain accounting */
        PRTE_RELEASE(sig);
    }
       
    make_timestamp_base(&cur_master_timing_frame->rc_publish_start);
    /* Inform daemons about res change so they can answer queries */
    cmd = PRTE_DYNRES_DEFINE_RES_CHANGE;
    pmix_data_buffer_t *buf2;
    for(n=0; n < ndaemons; n++){
        PMIX_DATA_BUFFER_CREATE(buf2);
        daemon_procid.rank = n;

        ret = PMIx_Data_pack(NULL, buf2, &cmd, 1, PMIX_UINT8);

        ret = PMIx_Data_pack(NULL, buf2, &rc_type, 1, PMIX_UINT8);
        
        ret = PMIx_Data_pack(NULL, buf2, (void*)&delta_pset_name, 1, PMIX_STRING);
        prte_rml.send_buffer_nb(&daemon_procid, buf2, PRTE_RML_TAG_MALLEABILITY, prte_rml_send_callback, NULL);
        //printf("send rc define to daemon: %d\n", n);
    }

    //prte_state_caddy_t *cd = PRTE_NEW(prte_state_caddy_t);
    //cd->jdata = job_data;
    //cd->job_state = (rc_type == PMIX_RES_CHANGE_ADD) ? PRTE_JOB_STATE_LAUNCH_APPS : PRTE_JOB_STATE_SUB;


    make_timestamp_base(&cur_master_timing_frame->apply_start);    
    /* retrieve the data needed by the launcher */
    if(rc_type == PMIX_RES_CHANGE_ADD){
        prte_state_caddy_t *cd = PRTE_NEW(prte_state_caddy_t);
        cd->jdata = job_data;
        cd->job_state = (rc_type == PMIX_RES_CHANGE_ADD) ? PRTE_JOB_STATE_LAUNCH_APPS : PRTE_JOB_STATE_SUB;
        prte_plm_base_launch_apps(0,0, cd);
    }
    //if(rc_type == PMIX_RES_CHANGE_SUB){
    //    /* message goes to all daemons */
    //    sig = PRTE_NEW(prte_grpcomm_signature_t);
    //    sig->signature = (pmix_proc_t *) malloc(sizeof(pmix_proc_t));
    //    PMIX_LOAD_PROCID(&sig->signature[0], PRTE_PROC_MY_NAME->nspace, PMIX_RANK_WILDCARD);
    //    sig->sz = 1;
    //    if (PRTE_SUCCESS != (rc = prte_grpcomm.xcast(sig, PRTE_RML_TAG_DAEMON, &job_data->launch_msg))) {
    //        PRTE_ERROR_LOG(rc);
    //        PRTE_RELEASE(sig);
    //        return;
    //    }
    //    PMIX_DATA_BUFFER_DESTRUCT(&job_data->launch_msg);
    //    PMIX_DATA_BUFFER_CONSTRUCT(&job_data->launch_msg);
    //    /* maintain accounting */
    //    PRTE_RELEASE(sig);
    //}

    /*
    pmix_data_buffer_t *buf3;
    for(n=0; n < ndaemons; n++){
        PMIX_DATA_BUFFER_CREATE(buf3);
        PMIX_DATA
        *buf3 = job_data->launch_msg;
        daemon_procid.rank = n;
       
        prte_rml.send_buffer_nb(&daemon_procid, buf3, PRTE_RML_TAG_DAEMON, prte_rml_send_callback, NULL);
    }
    */
    
    /* maintain accounting */
    //PRTE_RELEASE(sig);    
    free(delta_procs);
    free(delta_pset_name);
    //memset(&job_data->launch_msg, 0, sizeof(pmix_data_buffer_t));
    //PMIX_DATA_BUFFER_CONSTRUCT(&job_data->launch_msg);
    ///* We created a copy of the job data so release it here */
    //if(rc_type == PMIX_RES_CHANGE_SUB){
    //    PRTE_RELEASE(job_data);
    //}

    scd->evncbfunc(PMIX_SUCCESS, NULL, 0, NULL, NULL, scd->cbdata);
    make_timestamp_base(&cur_master_timing_frame->rc_end1);
}

static void rchandler(size_t evhdlr_registration_id, pmix_status_t status,
                       const pmix_proc_t *source, pmix_info_t info[], size_t ninfo,
                       pmix_info_t results[], size_t nresults,
                       pmix_event_notification_cbfunc_fn_t cbfunc, void *cbdata){
    prte_pmix_server_op_caddy_t *cd;
    cd = PRTE_NEW(prte_pmix_server_op_caddy_t);
    cd->proc = *source;
    cd->ev.ev_base = prte_event_base;
    cd->codes = &status;
    cd->ncodes = 1;
    cd->info = (pmix_info_t *) info;
    cd->ninfo = ninfo;
    cd->evncbfunc = cbfunc;
    cd->cbdata = cbdata;
    prte_event_set(prte_event_base, &(cd->ev), -1, PRTE_EV_WRITE, _rchandler, cd);
    prte_event_set_priority(&(cd->ev), PRTE_MSG_PRI);
    PRTE_POST_OBJECT(cd);
    prte_event_active(&(cd->ev), PRTE_EV_WRITE, 1);

}
