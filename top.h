/** 
 * @file   top.h
 * 
 * @brief  Top Level Commands
 *
 * 
 */

#ifndef _TOP_H_
#define _TOP_H_

void executecommand(int module, int index, int *nerr, char *kinmsg); //220306 Donglin Choi: LOAD EXTERNAL_COMMAND
void findcommand(char *kcommand, int *lfind, int *module, int *index);
void initcommon(void);
void initsac(void);
void initblkdata(void);
void setup_pager(void);
void sac_report_files_in_memory(int *nerr);
void saccommands(char *kinmsg, int kinmsg_s, int *nerr);

#endif /* _TOP_H_ */
