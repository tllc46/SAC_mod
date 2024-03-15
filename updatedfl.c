/** 
 * @file   updatedfl.c
 * 
 * @brief  Replace or append to the filelist
 * 
 */

#include <stdio.h>
#include <string.h>

#include "dfm.h"
#include "amf.h"
#include "hdr.h"

#include "extfunc.h"

#include "errors.h"

#include "msg.h"
#include "clf.h"
#include "ucf.h"
#include "dff.h"

#include <fstr.h>

DFM_EXTERN

/** 
 * Replace or append to the filelist in memory
 * 
 * @param call_data 
 *    Structure containing the descriptions of the sac files
 *    @see extfunc.h
 * @param update 
 *    - REPLACE to replace the current files in memory
 *    - FALSE to append to the current file list
 * @param nerr 
 *    Error Return Flag
 *    - 0 on Success
 *
 * @bug This routine assumes the size of the header does not change.
 *
 * @date   960229:  Original version.
 *
 */
void
updatedfl(sac_files call_data, int update, int *nerr, char (*filen)[MCPFN]) { //220306 Donglin Choi: LOAD store filenames

    char kfile[MCPFN + 1];
    char temp[call_data.nfiles][MCPFN]; //220306 Donglin Choi: LOAD store filenames
    int i;
    sac_header *this_header;
    float *ydata, *xdata;
    sac *s;
    *nerr = 0;
    for(i=0;i<call_data.nfiles;i++){ //220306 Donglin Choi: LOAD store filenames
	strcpy(temp[i],filen[i]); //220306 Donglin Choi
    } //220306 Donglin Choi

    if (update == REPLACE) {
        sacclear();
    }

    cmdfm.ndfl += call_data.nfiles;

    for (i = 0; i < call_data.nfiles; i++) {
        this_header = call_data.ext_hdrs[i];
        ydata = call_data.ext_yvalues[i];
        xdata = call_data.ext_xvalues[i];

        s = sac_new();
        sacput(s);
        memcpy(s->h, this_header, sizeof(sac_hdr));
        sac_alloc(s);
	sprintf(kfile, "%s", temp[i]); //220306 Donglin Choi: LOAD store filenames

        /* filename to storage */
        s->m->filename = fstrdup(kfile, MCPFN + 1);

        /* store the data */
        memcpy(s->y, ydata, s->h->npts * sizeof(float));
        if (s->x) {
            memcpy(s->x, xdata, s->h->npts * sizeof(float));
        }

        extrma(s->y, 1, s->h->npts, &s->h->depmin, &s->h->depmax,
               &s->h->depmen);

    }

    return;

}
