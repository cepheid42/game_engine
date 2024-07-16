//
// Created by cepheid on 7/13/24.
//

#ifndef SCHNIEDER_H
#define SCHNIEDER_H

struct Grid;

#define SizeX 10
#define SizeY 10
#define SizeZ 10

#define Hy1
#define Chyh1
#define Chye1
#define Ez1
#define Hx2
#define Chxh2
#define Chxe2
#define Ez2
#define Hy2
#define Chyh2
#define Chye2
#define Ez2
#define Hz2
#define Chzh2
#define Chze2
#define Ex2
#define Ey2
#define Hx
#define Chxh
#define Chxe
#define Hy
#define Chyh
#define Chye
#define Hz
#define Chzh
#define Chze
#define Ex
#define Ey
#define Ez
#define Cexh1
#define Cexe1
#define Ceyh1
#define Ceye1
#define Cezh1
#define Ceze1
#define Cexh2
#define Cexe2
#define Ceyh2
#define Ceye2
#define Cezh2
#define Ceze2
#define Cexh
#define Cexe
#define Ceyh
#define Ceye
#define Cezh
#define Ceze
enum GridType { oneDGrid, tmZGrid, teZGrid, threeDGrid };
inline GridType Type = tmZGrid;

template<size_t DIM>
struct


/* update magnetic field */
void updateH(Grid *g) {
  int mm, nn, pp;

  if (Type == oneDGrid) {
    for (mm = 0; mm < SizeX - 1; mm++)
      Hy1(mm) = Chyh1(mm) * Hy1(mm) + Chye1(mm) * (Ez1(mm + 1) - Ez1(mm));
  } else if (Type == tmZGrid) {
    for (mm = 0; mm < SizeX; mm++)
      for (nn = 0; nn < SizeY - 1; nn++)
        Hx2(mm, nn) = Chxh2(mm, nn) * Hx2(mm, nn) - Chxe2(mm, nn) * (Ez2(mm, nn + 1) - Ez2(mm, nn));

    for (mm = 0; mm < SizeX - 1; mm++)
      for (nn = 0; nn < SizeY; nn++)
        Hy2(mm, nn) = Chyh2(mm, nn) * Hy2(mm, nn) + Chye2(mm, nn) * (Ez2(mm + 1, nn) - Ez2(mm, nn));

  } else if (Type == teZGrid) {
    for (mm = 0; mm < SizeX - 1; mm++)
      for (nn = 0; nn < SizeY - 1; nn++)
        Hz2(mm, nn) = Chzh2(mm, nn) * Hz2(mm, nn) -
                      Chze2(mm, nn) * ((Ey2(mm + 1, nn) - Ey2(mm, nn)) -
                                       (Ex2(mm, nn + 1) - Ex2(mm, nn)));
  } else {
    for (mm = 0; mm < SizeX; mm++)
      for (nn = 0; nn < SizeY - 1; nn++)
        for (pp = 0; pp < SizeZ - 1; pp++)
          Hx(mm, nn, pp) = Chxh(mm, nn, pp) * Hx(mm, nn, pp) +
                           Chxe(mm, nn, pp) * ((Ey(mm, nn, pp + 1) - Ey(mm, nn, pp)) -
                                               (Ez(mm, nn + 1, pp) - Ez(mm, nn, pp)));

    for (mm = 0; mm < SizeX - 1; mm++)
      for (nn = 0; nn < SizeY; nn++)
        for (pp = 0; pp < SizeZ - 1; pp++)
          Hy(mm, nn, pp) = Chyh(mm, nn, pp) * Hy(mm, nn, pp) +
                           Chye(mm, nn, pp) * ((Ez(mm + 1, nn, pp) - Ez(mm, nn, pp)) -
                                               (Ex(mm, nn, pp + 1) - Ex(mm, nn, pp)));

    for (mm = 0; mm < SizeX - 1; mm++)
      for (nn = 0; nn < SizeY - 1; nn++)
        for (pp = 0; pp < SizeZ; pp++)
          Hz(mm, nn, pp) = Chzh(mm, nn, pp) * Hz(mm, nn, pp) +
                           Chze(mm, nn, pp) * ((Ex(mm, nn + 1, pp) - Ex(mm, nn, pp)) - (Ey
                                                 (mm + 1, nn, pp) - Ey(mm, nn, pp)));
  }
} /* end updateH() */

/* update electric field */
void updateE(Grid *g) {
  int mm, nn, pp;

  if (Type == oneDGrid) {
    for (mm = 1; mm < SizeX - 1; mm++)
      Ez1(mm) = Ceze1(mm) * Ez1(mm)
                + Cezh1(mm) * (Hy1(mm) - Hy1(mm - 1));
  } else if (Type == tmZGrid) {
    for (mm = 1; mm < SizeX - 1; mm++)
      for (nn = 1; nn < SizeY - 1; nn++)
        Ez2(mm, nn) = Ceze2(mm, nn) * Ez2(mm, nn) +
                      Cezh2(mm, nn) * ((Hy2(mm, nn) - Hy2(mm - 1, nn)) -
                                       (Hx2(mm, nn) - Hx2(mm, nn - 1)));
  } else if (Type == teZGrid) {
    for (mm = 1; mm < SizeX - 1; mm++)
      for (nn = 1; nn < SizeY - 1; nn++)
        Ex2(mm, nn) = Cexe2(mm, nn) * Ex2(mm, nn) +
                      Cexh2(mm, nn) * (Hz2(mm, nn) - Hz2(mm, nn - 1));

    for (mm = 1; mm < SizeX - 1; mm++)
      for (nn = 1; nn < SizeY - 1; nn++)
        Ey2(mm, nn) = Ceye2(mm, nn) * Ey2(mm, nn) -
                      Ceyh2(mm, nn) * (Hz2(mm, nn) - Hz2(mm - 1, nn));
  } else {
    for (mm = 0; mm < SizeX - 1; mm++)
      for (nn = 1; nn < SizeY - 1; nn++)
        for (pp = 1; pp < SizeZ - 1; pp++)
          Ex(mm, nn, pp) = Cexe(mm, nn, pp) * Ex(mm, nn, pp) +
                           Cexh(mm, nn, pp) * ((Hz(mm, nn, pp) - Hz(mm, nn - 1, pp)) - (Hy
                                                 (mm, nn, pp) - Hy(mm, nn, pp - 1)));

    for (mm = 1; mm < SizeX - 1; mm++)
      for (nn = 0; nn < SizeY - 1; nn++)
        for (pp = 1; pp < SizeZ - 1; pp++)
          Ey(mm, nn, pp) = Ceye(mm, nn, pp) * Ey(mm, nn, pp) +
                           Ceyh(mm, nn, pp) * ((Hx(mm, nn, pp) - Hx(mm, nn, pp - 1)) -
                                               (Hz(mm, nn, pp) - Hz(mm - 1, nn, pp)));

    for (mm = 1; mm < SizeX - 1; mm++)
      for (nn = 1; nn < SizeY - 1; nn++)
        for (pp = 0; pp < SizeZ - 1; pp++)
          Ez(mm, nn, pp) = Ceze(mm, nn, pp) * Ez(mm, nn, pp) +
                           Cezh(mm, nn, pp) * ((Hy(mm, nn, pp) - Hy(mm - 1, nn, pp)) -
                                               (Hx(mm, nn, pp) - Hx(mm, nn - 1, pp)));
  }
} /* end updateE() */

#endif //SCHNIEDER_H
