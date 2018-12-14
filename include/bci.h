/*
libbci is an interpolation library based on the bicubic interpolation method
Copyright (C) 2004 Andrew E. Shevtsov

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
*/
#ifndef __BCI_H__
#define __BCI_H__

#if defined(__cplusplus)
#define __BEGIN_DECLS extern "C" {
#define __END_DECLS   }
#else
#define __BEGIN_DECLS
#define __END_DECLS
#endif

#define pow_2(x)	((x)*(x))
#define pow_3(x)	((x)*(x)*(x))

struct xyz {
    double x, y, z;
};

struct xyzwder {
    double x, y, z, dx, dy, dxdy;
};

struct der {
    double dx, dy, dxdy;
};

typedef struct der doubleder;
typedef struct xyz doublexyz;
typedef struct xyzwder doublexyzwder;

__BEGIN_DECLS
/*
 auxiliary functions (tools.c)
 */
void printm(char *title, double *data, int cols, int rows);

void triangm(double *mtrx, int cols, double *g);

double pow_int(double x, int n);

/*
 low level functions (tdspl.c)
 */
void td_solve4by4(double c[16], doublexyz x[16]);

void td_solveonespline4(double c[16], const doublexyzwder xyzwder[4]);

/*
 high level functions (bci.c)
 */
void td_fillgrid(const doublexyz *xyz, int cols, int rows,
                 doublexyz *xyzf, int fcols, int frows);

void td_fillgrid_wo_nulls(const doublexyz *xyz, int cols, int rows,
                          doublexyz *xyzf, int fcols, int frows);

__END_DECLS
#endif /* __BCI_H__ */
