
#ifndef CORRESP_H
#define CORRESP_H

#include "epi.h"

typedef struct
{
  int     *points; /* one per camera, indexing the targets for that camera */
  double  corr;    /* normalized feature based correlation coefficient */
}
n_tupel;

typedef struct
{
  int n;          /* # of candidates */
  int p2[maxcand];    /* point numbers of candidates */
  double corr[maxcand];  /* feature based correlation coefficient */
  double dist[maxcand];  /* distance perpendicular to epipolar line */
} 
correspond;         /* correspondence candidates */

#endif
