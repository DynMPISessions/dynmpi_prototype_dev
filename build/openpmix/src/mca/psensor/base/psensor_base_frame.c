/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil -*- */
/*
 * Copyright (c) 2010      Cisco Systems, Inc.  All rights reserved.
 * Copyright (c) 2012-2013 Los Alamos National Security, Inc. All rights reserved.
 * Copyright (c) 2017-2020 Intel, Inc.  All rights reserved.
 * Copyright (c) 2020      Research Organization for Information Science
 *                         and Technology (RIST).  All rights reserved.
 * Copyright (c) 2021-2022 Nanook Consulting  All rights reserved.
 * $COPYRIGHT$
 *
 * Additional copyrights may follow
 *
 * $HEADER$
 */

#include "src/include/pmix_config.h"

#include "pmix_common.h"

#include <pthread.h>
#include <event.h>

#include "src/class/pmix_list.h"
#include "src/include/pmix_types.h"
#include "src/mca/base/pmix_base.h"
#include "src/mca/mca.h"
#include "src/runtime/pmix_progress_threads.h"

#include "src/mca/psensor/base/base.h"

/*
 * The following file was created by configure.  It contains extern
 * statements and the definition of an array of pointers to each
 * component's public mca_base_component_t struct.
 */

#include "src/mca/psensor/base/static-components.h"

/*
 * Global variables
 */
pmix_psensor_base_module_t pmix_psensor = {
    .start = pmix_psensor_base_start,
    .stop = pmix_psensor_base_stop
};

pmix_psensor_base_t pmix_psensor_base = {
    .actives = PMIX_LIST_STATIC_INIT,
    .evbase = NULL,
    .selected = false
};

static bool use_separate_thread = false;

static int pmix_psensor_register(pmix_mca_base_register_flag_t flags)
{
    (void) flags;
    (void) pmix_mca_base_var_register("pmix", "psensor", "base", "use_separate_thread",
                                      "Use a separate thread for monitoring local procs",
                                      PMIX_MCA_BASE_VAR_TYPE_BOOL,
                                      &use_separate_thread);
    return PMIX_SUCCESS;
}

static int pmix_psensor_base_close(void)
{
    pmix_psensor_base.selected = false;
    PMIX_LIST_DESTRUCT(&pmix_psensor_base.actives);

    if (use_separate_thread && NULL != pmix_psensor_base.evbase) {
        (void) pmix_progress_thread_stop("PSENSOR");
    }

    /* Close all remaining available components */
    return pmix_mca_base_framework_components_close(&pmix_psensor_base_framework, NULL);
}

/**
 * Function for finding and opening either all MCA components, or the one
 * that was specifically requested via a MCA parameter.
 */
static int pmix_psensor_base_open(pmix_mca_base_open_flag_t flags)
{
    /* construct the list of modules */
    PMIX_CONSTRUCT(&pmix_psensor_base.actives, pmix_list_t);

    if (use_separate_thread) {
        /* create an event base and progress thread for us */
        pmix_psensor_base.evbase = pmix_progress_thread_init("PSENSOR");
        if (NULL == pmix_psensor_base.evbase) {
            return PMIX_ERROR;
        }

    } else {
        pmix_psensor_base.evbase = pmix_globals.evbase;
    }

    /* Open up all available components */
    return pmix_mca_base_framework_components_open(&pmix_psensor_base_framework, flags);
}

PMIX_MCA_BASE_FRAMEWORK_DECLARE(pmix, psensor, "PMIx Monitoring Sensors", pmix_psensor_register,
                                pmix_psensor_base_open, pmix_psensor_base_close,
                                pmix_mca_psensor_base_static_components,
                                PMIX_MCA_BASE_FRAMEWORK_FLAG_DEFAULT);

PMIX_CLASS_INSTANCE(pmix_psensor_active_module_t, pmix_list_item_t, NULL, NULL);
