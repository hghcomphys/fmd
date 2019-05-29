/*
  forces.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani, Hossein Ghorbanfekr

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include "forces.h"
#include "base.h"
#include "cspline.h"
#include "md_ghost.h"
#include "list.h"

#ifdef USE_CSPLINE
#define EAM_PAIR_UPDATE_rho_host                                       \
    {                                                                  \
        COMPUTE_rv_AND_r2;                                             \
                                                                       \
        eam = (eam_t *)pottable[atomkind1][atomkind2].data;            \
                                                                       \
        if (r2 < eam->cutoff_sqr)                                      \
        {                                                              \
            h = eam->dr2;                                              \
            ir2 = (int)(r2 / h);                                       \
            ir2_h = ir2 + 1;                                           \
            a = ir2_h - r2/h;                                          \
            b=1-a;                                                     \
                                                                       \
            unsigned jloc = pottable[atomkind1][atomkind2].jloc;       \
            rho = eam->elements[jloc].rho;                             \
            rhoDD = eam->elements[jloc].rhoDD;                         \
            rho_host += SPLINE_VAL(a,b,rho,ir2,ir2_h,rhoDD,h);         \
        }                                                              \
    }
#else
#define EAM_PAIR_UPDATE_rho_host                                       \
    {                                                                  \
        COMPUTE_rv_AND_r2;                                             \
                                                                       \
        eam = (eam_t *)pottable[atomkind1][atomkind2].data;            \
                                                                       \
        if (r2 < eam->cutoff_sqr)                                      \
        {                                                              \
            h = eam->dr2;                                              \
            ir2 = (int)(r2 / h);                                       \
            ir2_h = ir2 + 1;                                           \
            a = ir2_h - r2/h;                                          \
            b=1-a;                                                     \
                                                                       \
            unsigned jloc = pottable[atomkind1][atomkind2].jloc;       \
            rho = eam->elements[jloc].rho;                             \
            rho_host += rho[ir2]*a + rho[ir2_h]*b;                     \
        }                                                              \
    }
#endif

#ifdef USE_CSPLINE
#define EAM_COMPUTE_FembPrime_AND_UPDATE_Femb_sum                                \
    {                                                                            \
        if (rho_host == 0.0)                                                     \
        { /* the first atom didn't see any other atom within its cutoff sphere */\
            Femb_sum += sysp->potsys.atomkinds[atomkind1].aux->eam_F0;           \
        }                                                                        \
        else                                                                     \
        {                                                                        \
            h = eam->drho;                                                       \
            irho = (int)(rho_host / h);                                          \
            assert(irho < eam->Nrho - 1);                                        \
            irho_h = irho + 1;                                                   \
                                                                                 \
            unsigned iloc = pottable[atomkind1][atomkind2].iloc;                 \
            F = eam->elements[iloc].F;                                           \
            F_DD = eam->elements[iloc].F_DD;                                     \
            a = irho_h - rho_host/h;                                             \
            b = 1-a;                                                             \
            item1_p->FembPrime = SPLINE_DERIV(a,b,F,irho,irho_h,F_DD,h);         \
            Femb_sum += SPLINE_VAL(a,b,F,irho,irho_h,F_DD,h);                    \
        }                                                                        \
    }
#else
#define EAM_COMPUTE_FembPrime_AND_UPDATE_Femb_sum                                \
    {                                                                            \
        if (rho_host == 0.0)                                                     \
        { /* the first atom didn't see any other atom within its cutoff sphere */\
            Femb_sum += sysp->potsys.atomkinds[atomkind1].aux->eam_F0;           \
        }                                                                        \
        else                                                                     \
        {                                                                        \
            h = eam->drho;                                                       \
            irho = (int)(rho_host / h);                                          \
            assert(irho < eam->Nrho - 1);                                        \
            irho_h = irho + 1;                                                   \
                                                                                 \
            unsigned iloc = pottable[atomkind1][atomkind2].iloc;                 \
            F = eam->elements[iloc].F;                                           \
            item1_p->FembPrime = (F[irho_h] - F[irho]) / h;                      \
            Femb_sum += F[irho] + (rho_host - irho * h) * item1_p->FembPrime;    \
        }                                                                        \
    }
#endif

static void compute_hybrid_pass1(fmd_sys_t *sysp, double *FembSum_p)
{
    int jc[3], kc[3];
    int d, ir2, irho, ir2_h, irho_h;
    TParticleListItem *item1_p, *item2_p;
    double r2, rv[3];
    double *rho, *rhoDD, *F, *F_DD;
    double a, b, h;
    int ic0, ic1, ic2;
    double Femb_sum=0.0;
    potpair_t **pottable = sysp->potsys.pottable;

    // iterate over all cells(lists)
    #pragma omp parallel for private(ic0,ic1,ic2,item1_p,kc,jc,item2_p,d,rv,r2,h,ir2,ir2_h,a,b,rho, \
      rhoDD,F,F_DD,irho,irho_h) shared(sysp,pottable) default(none) collapse(3) reduction(+:Femb_sum) \
      schedule(static,1)
    for (ic0 = sysp->subDomain.ic_start[0]; ic0 < sysp->subDomain.ic_stop[0]; ic0++)
    for (ic1 = sysp->subDomain.ic_start[1]; ic1 < sysp->subDomain.ic_stop[1]; ic1++)
    for (ic2 = sysp->subDomain.ic_start[2]; ic2 < sysp->subDomain.ic_stop[2]; ic2++)
    {
        // iterate over all items in cell ic
        for (item1_p = sysp->subDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            if (!(sysp->activeGroup == -1 || item1_p->P.groupID == sysp->activeGroup))
                continue;

            eam_t *eam;
            unsigned atomkind1, atomkind2;
            atomkind1 = item1_p->P.elementID;

            if (!sysp->potsys.atomkinds[atomkind1].usesEAM)
                continue;

            double rho_host = 0.0;
            // iterate over neighbor cells of cell ic
            for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0)
                for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1)
                    for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2)
                        // iterate over all items in cell jc
                        for (item2_p = sysp->subDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
                        {
                            if (!(sysp->activeGroup == -1 || item2_p->P.groupID == sysp->activeGroup))
                                continue;

                            if (item1_p != item2_p)
                            {
                                atomkind2 = item2_p->P.elementID;
                                if (pottable[atomkind1][atomkind2].kind == POTKIND_EAM_ALLOY)
                                    EAM_PAIR_UPDATE_rho_host;
                            }
                        }
                    }
                }
            }

            EAM_COMPUTE_FembPrime_AND_UPDATE_Femb_sum;
        }
    }

    *FembSum_p=Femb_sum;
}

static void computeEAM_pass1(fmd_sys_t *sysp, double *FembSum_p)
{
    int jc[3], kc[3];
    int d, ir2, irho, ir2_h, irho_h;
    TParticleListItem *item1_p, *item2_p;
    double r2, rv[3];
    double *rho, *rhoDD, *F, *F_DD;
    double a, b, h;
    int ic0, ic1, ic2;
    double Femb_sum=0;
    potpair_t **pottable = sysp->potsys.pottable;

    // iterate over all cells(lists)
    #pragma omp parallel for private(ic0,ic1,ic2,item1_p,kc,jc,item2_p,d,rv,r2,h,ir2,ir2_h,a,b,rho, \
      rhoDD,F,F_DD,irho,irho_h) shared(sysp,pottable) default(none) collapse(3) reduction(+:Femb_sum) \
      schedule(static,1)
    for (ic0 = sysp->subDomain.ic_start[0]; ic0 < sysp->subDomain.ic_stop[0]; ic0++)
    for (ic1 = sysp->subDomain.ic_start[1]; ic1 < sysp->subDomain.ic_stop[1]; ic1++)
    for (ic2 = sysp->subDomain.ic_start[2]; ic2 < sysp->subDomain.ic_stop[2]; ic2++)
    {
        // iterate over all items in cell ic
        for (item1_p = sysp->subDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            if (!(sysp->activeGroup == -1 || item1_p->P.groupID == sysp->activeGroup))
                continue;

            eam_t *eam;
            unsigned atomkind1, atomkind2;
            atomkind1 = item1_p->P.elementID;

            double rho_host = 0.0;
            // iterate over neighbor cells of cell ic
            for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0)
                for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1)
                    for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2)
                        // iterate over all items in cell jc
                        for (item2_p = sysp->subDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
                        {
                            if (!(sysp->activeGroup == -1 || item2_p->P.groupID == sysp->activeGroup))
                                continue;

                            if (item1_p != item2_p)
                            {
                                atomkind2 = item2_p->P.elementID;
                                EAM_PAIR_UPDATE_rho_host;
                            }
                        }
                    }
                }
            }

            EAM_COMPUTE_FembPrime_AND_UPDATE_Femb_sum;
        }
    }
    *FembSum_p=Femb_sum;
}

static void computeEAM_pass0(fmd_sys_t *sysp, double FembSum)
{
    int jc[3], kc[3];
    int d, ir2, ir2_h;
    TParticleListItem *item1_p, *item2_p;
    double r2, rv[3];
    double *rho_i, *rho_j, *phi;
    double *rho_iDD, *rho_jDD, *phiDD;
    double rho_ip, rho_jp;
    double mag;
    double sum[3];
    double phi_deriv;
    double a, b, h;
    int ic0, ic1, ic2;
    potpair_t **pottable = sysp->potsys.pottable;
#ifdef USE_TTM
    double mass;
    int ttm_index;
    double dx;
    double pxx = 0.0;
#endif
    double potEnergy = 0.0;

    // iterate over all cells(lists)
#ifdef USE_TTM
    #pragma omp parallel for private(ic0,ic1,ic2,ttm_index,item1_p,d,element_i,rho_i,rho_iDD,kc,jc,item2_p,rv,r2,h,ir2, \
      ir2_h,element_j,phi,phiDD,a,b,phi_deriv,rho_ip,rho_jp,rho_jDD,rho_j,mag,mass,dx,sum) \
      shared(sysp,ttm_lattice_aux,ttm_useSuction,ttm_suctionWidth,ttm_suctionIntensity,ttm_pxx_compute, \
      ttm_pxx_pos) default(none) collapse(3) reduction(+:potEnergy,pxx) schedule(static,1)
#else
    #pragma omp parallel for private(ic0,ic1,ic2,item1_p,d,rho_i,rho_iDD,kc,jc,item2_p,rv,r2,h,ir2, \
      ir2_h,phi,phiDD,a,b,phi_deriv,rho_ip,rho_jp,rho_jDD,rho_j,mag,sum) \
      shared(sysp,pottable) default(none) collapse(3) reduction(+:potEnergy) schedule(static,1)
#endif
    for (ic0 = sysp->subDomain.ic_start[0]; ic0 < sysp->subDomain.ic_stop[0]; ic0++)
    for (ic1 = sysp->subDomain.ic_start[1]; ic1 < sysp->subDomain.ic_stop[1]; ic1++)
    for (ic2 = sysp->subDomain.ic_start[2]; ic2 < sysp->subDomain.ic_stop[2]; ic2++)
    {
#ifdef USE_TTM
        ttm_index = ic0 - sysp->subDomain.ic_start[0] + 1;
#endif
        // iterate over all items in cell ic
        for (item1_p = sysp->subDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            if (!(sysp->activeGroup == -1 || item1_p->P.groupID == sysp->activeGroup))
                continue;

            for (d=0; d<3; d++)
                sum[d] = 0.0;

            eam_t *eam;
            unsigned atomkind1, atomkind2;
            atomkind1 = item1_p->P.elementID;

            // iterate over neighbor cells of cell ic
            for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0)
                for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1)
                    for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2)
                        // iterate over all items in cell jc
                        for (item2_p = sysp->subDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
                        {
                            if (!(sysp->activeGroup == -1 || item2_p->P.groupID == sysp->activeGroup))
                                continue;
                            if (item1_p != item2_p)
                            {
                                COMPUTE_rv_AND_r2;

                                atomkind2 = item2_p->P.elementID;
                                eam = (eam_t *)pottable[atomkind1][atomkind2].data;

                                if (r2 < eam->cutoff_sqr)
                                {
                                    h = eam->dr2;
                                    ir2 = (int)(r2 / h);
                                    ir2_h = ir2 + 1;
                                    unsigned iloc = pottable[atomkind1][atomkind2].iloc;
                                    unsigned jloc = pottable[atomkind1][atomkind2].jloc;

                                    rho_i = eam->elements[iloc].rho;
#ifdef USE_CSPLINE
                                    rho_iDD = eam->elements[iloc].rhoDD;
#endif
                                    phi = eam->elements[iloc].phi[jloc];
#ifdef USE_CSPLINE
                                    phiDD = eam->elements[iloc].phiDD[jloc];
                                    a = ir2_h - r2/h;
                                    b = 1-a;
                                    phi_deriv = SPLINE_DERIV(a,b,phi,ir2,ir2_h,phiDD,h);
                                    rho_ip = SPLINE_DERIV(a,b,rho_i,ir2,ir2_h,rho_iDD,h);
                                    if (jloc == iloc)
                                        rho_jp = rho_ip;
                                    else
                                    {
                                        rho_j = eam->elements[jloc].rho;
                                        rho_jDD = eam->elements[jloc].rhoDD;
                                        rho_jp = SPLINE_DERIV(a,b,rho_j,ir2,ir2_h,rho_jDD,h);
                                    }

                                    mag = 2 * (item1_p->FembPrime * rho_jp +
                                          item2_p->FembPrime * rho_ip + phi_deriv);
                                    potEnergy += SPLINE_VAL(a,b,phi,ir2,ir2_h,phiDD,h);
#else
                                    rho_j = eam->elements[jloc].rho;
                                    mag = 2 * (item1_p->FembPrime * (rho_j[ir2_h] - rho_j[ir2]) +
                                               item2_p->FembPrime * (rho_i[ir2_h] - rho_i[ir2]) +
                                                                        (phi[ir2_h] - phi[ir2])) / h;
                                    potEnergy += phi[ir2] + (r2/h - ir2) * (phi[ir2_h] - phi[ir2]);
#endif
                                    for (d=0; d<3; d++)
                                        sum[d] += mag * rv[d];
                                }
                            }
                        }
                    }
                }
            }

#ifdef USE_TTM
            mass = sysp->EAM.elements[element_i].mass;
            for (d=0; d<3; d++)
                item1_p->F[d] = -sum[d] + ttm_lattice_aux[ttm_index].xi *
                    mass * (item1_p->P.v[d] - ttm_lattice_aux[ttm_index].v_cm[d]);
            if (ttm_useSuction)
                if (item1_p->P.x[0] < ttm_suctionWidth)
                    item1_p->F[0] -= mass * ttm_suctionIntensity;
            if (ttm_pxx_compute)
            {
                dx = item1_p->P.x[0] - ttm_pxx_pos;
                pxx += item1_p->F[0] * ((dx > 0) - (dx < 0));
            }
#else
            for (d=0; d<3; d++)
                item1_p->F[d] = -sum[d];
#endif
        }
    }
#ifdef USE_TTM
    ttm_pxx_local[1] += pxx;
#endif
    potEnergy = 0.5 * potEnergy + FembSum;
    MPI_Allreduce(&potEnergy, &sysp->totalPotentialEnergy, 1, MPI_DOUBLE, MPI_SUM, sysp->MD_comm);
}

static void computeLJ(fmd_sys_t *sysp)
{
    int jc[3], kc[3];
    int d;
    TParticleListItem *item1_p, *item2_p;
    double r2, rv[3];
    int ic0, ic1, ic2;
    double potEnergy = 0.0;
    double F[3];
    potpair_t **pottable = sysp->potsys.pottable;

    // iterate over all cells(lists)
    #pragma omp parallel for private(ic0,ic1,ic2,item1_p,d,kc,jc,item2_p,rv,r2,F) \
      shared(sysp,pottable) default(none) collapse(3) reduction(+:potEnergy) schedule(static,1)
    for (ic0 = sysp->subDomain.ic_start[0]; ic0 < sysp->subDomain.ic_stop[0]; ic0++)
        for (ic1 = sysp->subDomain.ic_start[1]; ic1 < sysp->subDomain.ic_stop[1]; ic1++)
            for (ic2 = sysp->subDomain.ic_start[2]; ic2 < sysp->subDomain.ic_stop[2]; ic2++)
            {
                // iterate over all items in cell ic
                for (item1_p = sysp->subDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
                {
                    if (!(sysp->activeGroup == -1 || item1_p->P.groupID == sysp->activeGroup))
                        continue;

                    unsigned atomkind1 = item1_p->P.elementID;

                    for (d=0; d<3; d++)
                        F[d] = 0.0;

                    // iterate over neighbor cells of cell ic
                    for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
                    {
                        SET_jc_IN_DIRECTION(0)
                        for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                        {
                            SET_jc_IN_DIRECTION(1)
                            for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                            {
                                SET_jc_IN_DIRECTION(2)
                                // iterate over all items in cell jc
                                for (item2_p = sysp->subDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
                                {
                                    if (!(sysp->activeGroup == -1 || item2_p->P.groupID == sysp->activeGroup))
                                        continue;

                                    if (item1_p != item2_p)
                                    {
                                        COMPUTE_rv_AND_r2;

                                        unsigned atomkind2 = item2_p->P.elementID;
                                        LJ_6_12_t *lj = (LJ_6_12_t *)pottable[atomkind1][atomkind2].data;

                                        if (r2 < lj->cutoff_sqr)
                                        {
                                            double inv_r2, inv_rs2, inv_rs6, inv_rs12;

                                            // force, F = -(d/dr)U
                                            inv_r2 = 1.0/r2;
                                            inv_rs2 = SQR(lj->sig) * inv_r2;
                                            inv_rs6 = inv_rs2 * inv_rs2 * inv_rs2;
                                            inv_rs12 = SQR(inv_rs6);
                                            double factor = lj->eps * inv_r2 * (inv_rs12 - 0.5*inv_rs6);
                                            for (d=0; d<3; d++)
                                                F[d] += rv[d] * factor;

                                            // potential energy, U = 4*eps*( (sig/r)^12 - (sig/r)^6 )
                                            potEnergy += lj->eps * (inv_rs12 - inv_rs6);

                                            // pressure
                                            /*for (d=0; d<3; d++) {
                                                pressure += -F[d]*rv[d]
                                            }*/
                                        }
                                    }
                                }
                            }
                        }
                    }

                    for (d=0; d<3; d++)
                        item1_p->F[d] = 48.0 * F[d];
                }
            }

    potEnergy *= 2.0;
    MPI_Allreduce(&potEnergy, &sysp->totalPotentialEnergy, 1, MPI_DOUBLE, MPI_SUM, sysp->MD_comm);
}

static void computeMorse(fmd_sys_t *sysp)
{
    int jc[3], kc[3];
    int d;
    TParticleListItem *item1_p, *item2_p;
    double r2, rv[3];
    int ic0, ic1, ic2;
    double potEnergy = 0.0;
    double F[3];
    potpair_t **pottable = sysp->potsys.pottable;

    // iterate over all cells(lists)
#pragma omp parallel for private(ic0,ic1,ic2,item1_p,d,kc,jc,item2_p,rv,r2,F) \
      shared(sysp,pottable) default(none) collapse(3) reduction(+:potEnergy) schedule(static,1)
    for (ic0 = sysp->subDomain.ic_start[0]; ic0 < sysp->subDomain.ic_stop[0]; ic0++)
        for (ic1 = sysp->subDomain.ic_start[1]; ic1 < sysp->subDomain.ic_stop[1]; ic1++)
            for (ic2 = sysp->subDomain.ic_start[2]; ic2 < sysp->subDomain.ic_stop[2]; ic2++)
            {
                // iterate over all items in cell ic
                for (item1_p = sysp->subDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
                {
                    if (!(sysp->activeGroup == -1 || item1_p->P.groupID == sysp->activeGroup))
                        continue;

                    unsigned atomkind1 = item1_p->P.elementID;

                    for (d=0; d<3; d++)
                        F[d] = 0.0;

                    // iterate over neighbor cells of cell ic
                    for (kc[0]=ic0-1; kc[0]<=ic0+1; kc[0]++)
                    {
                        SET_jc_IN_DIRECTION(0)
                        for (kc[1]=ic1-1; kc[1]<=ic1+1; kc[1]++)
                        {
                            SET_jc_IN_DIRECTION(1)
                            for (kc[2]=ic2-1; kc[2]<=ic2+1; kc[2]++)
                            {
                                SET_jc_IN_DIRECTION(2)
                                // iterate over all items in cell jc
                                for (item2_p = sysp->subDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
                                {
                                    if (!(sysp->activeGroup == -1 || item2_p->P.groupID == sysp->activeGroup))
                                        continue;

                                    if (item1_p != item2_p)
                                    {
                                        COMPUTE_rv_AND_r2;

                                        unsigned atomkind2 = item2_p->P.elementID;
                                        morse_t *morse = (morse_t *)pottable[atomkind1][atomkind2].data;

                                        if (r2 < morse->cutoff_sqr)
                                        {
                                            // force, F = -(d/dr)U
                                            double r = sqrt(r2);
                                            double inv_r = 1.0/r;
                                            double exp1 = exp( -morse->alpha * (r - morse->r0) );
                                            double exp2 = SQR(exp1);

                                            for (d=0; d<3; d++)
                                                F[d] +=  2.0 * morse->alpha * morse->D0 * rv[d] * inv_r * (exp2 - exp1);

                                            // potential energy, U = D0 * ( exp(-2*alpha*(r-r0)) - 2*exp(-alpha*(r-r0)) )
                                            potEnergy += morse->D0 * (exp2 - 2.0 * exp1);

                                            // pressure
                                            /*for (d=0; d<3; d++) {
                                                pressure += -F[d]*rv[d]
                                            }*/
                                        }
                                    }
                                }
                            }
                        }
                    }
                    for (d=0; d<3; d++)
                        item1_p->F[d] = F[d];
                }
            }

    potEnergy *= 0.5;  /*correct double-counting*/
    MPI_Allreduce(&potEnergy, &sysp->totalPotentialEnergy, 1, MPI_DOUBLE, MPI_SUM, sysp->MD_comm);
}


void fmd_dync_updateForces(fmd_sys_t *sysp)
{
    if (sysp->potsys.potkinds == NULL)  // just for one time
        fmd_pot_prepareForForceComp(sysp);

    fmd_ghostparticles_init(sysp);

    if (sysp->potsys.potkinds_num == 1)
    {
        potkind_t potkind = *(potkind_t *)(sysp->potsys.potkinds->data);

        switch (potkind)
        {
            case POTKIND_LJ_6_12:
                computeLJ(sysp);
                break;

            case POTKIND_MORSE:
                computeMorse(sysp);
                break;

            case POTKIND_EAM_ALLOY:
                if (sysp->iCompLocOrdParam) compLocOrdParam(sysp);
                double FembSum;
                computeEAM_pass1(sysp, &FembSum);
                fmd_ghostparticles_update_Femb(sysp);
                computeEAM_pass0(sysp, FembSum);
                break;
        }
    }

    fmd_ghostparticles_delete(sysp);
}
