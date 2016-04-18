#include <time.h>

void udate_(int pnow[]);

/*****************************************************************************/
/*                                                                           */
/* Routine:        udate_                                                    */
/*                                                                           */
/* Purpose:        Return date information in integer format                 */
/*                                                                           */
/* Revision:       0 : Initial Documented Version, 8/89, DLM                 */
/*                                                                           */
/* Usage:          udate_( pnow )                                            */
/*                 int pnow[6]; (O: Integer date information)                */
/*                                                                           */
/* Notes:          udate_() returns integer date information obtained from a */
/*                 standard UNIX time()/localtime() calls in the order year, */
/*                 month, day of month, hour, minute, second.                */
/*                                                                           */
/*****************************************************************************/


void udate_(int pnow[])
{
  long sec_since_1_1_70;
  struct tm *local_time;

  /* Initialize clock */

  time(&sec_since_1_1_70);
  local_time = localtime(&sec_since_1_1_70);
  
  /*  Year, month, day  */
  pnow[0] = 1900 + local_time->tm_year;
  pnow[1] = 1    + local_time->tm_mon;
  pnow[2] = local_time->tm_mday;

  /*  Hour, minute, second  */
  pnow[3] = local_time->tm_hour;
  pnow[4] = local_time->tm_min;
  pnow[5] = local_time->tm_sec;
}
