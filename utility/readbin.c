/* Read compressed HISTORY file with */
/* traj 1 nstep 1 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/param.h>

main()
{
 char filename[80],title[80], dummy[15],dummy2[80];
 char buffer[80];
 int timestep, total, dummyi, stat, n, conf=0, save_ok;
 double step=0.0, dummyf;
 double x, y, z;
 FILE *fp, *out;

 out=fopen("output","w");

/* Use system command zcat to read compressed file */
 sprintf(buffer,"/usr/bin/zcat %s\0", "HISTORY.Z");

/* Retrieve data from buffer via pipe */
 fp=popen(buffer,"r");

/* Get titles time step etc */
 fgets(dummy2,80,fp);
 fgets(dummy2,80,fp);

 for(;;) {

    if((stat=fscanf(fp,"%s   %d   %d   %d   %d  %lf\n",
                   &dummy,&timestep,&total,&dummyi,&dummyi,&step)) != EOF) {
       printf("Timestep=%d (%f ps)   Total=%d\n",
               timestep, (double) timestep*step, total);

       conf++;

/* Save every 10 configuration (to make xmol movie) */
       if((conf % 10) == 0)
          save_ok=1;
       else
          save_ok=0;

       if(save_ok)
          fprintf(out," %d  time=%d\n\n", total,timestep);
       for(n=1;n<=total;n++) {
          fgets(dummy, 15, fp);                      /* read atom label */
          fscanf(fp,"%d  %lf  %lf\n", &dummyi, &dummyf, &dummyf); /* read atom number, mass, charge */
          fscanf(fp,"%lf  %lf  %lf\n", &x, &y, &z);  /* Read coor */
          fscanf(fp,"%lf  %lf  %lf\n", &dummyf, &dummyf, &dummyf); /* Read force */
          if(!save_ok) continue;
          fprintf(out,"%c  %lf  %lf  %lf\n", dummy[0], x, y, z);

          /* some other analysis on x, y, z coordinates of a particular */
          /* system configurations */

       }  /* end for */
    }     /* end if */
    else {
       pclose(fp);
       fclose(out);
       break;
    }
 }  /* end for */
}
