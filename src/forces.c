/*
  forces.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani, Hossein Ghorbanfekr Kalashami

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

static void computeEAM_pass1(fmd_sys_t *sysp, double *FembSum_p)
{
/*
    int jc[3], kc[3];
    int d, ir2, irho, ir2_h, irho_h;
    TParticleListItem *item1_p, *item2_p;
    const double r_cutSqd = SQR(sysp->EAM.cutoff);
    double r2, rv[3];
    double rho_host, *rho, *rhoDD, *F, *F_DD;
    double a, b, h;
    int ic0, ic1, ic2;
    double sum=0;

    // iterate over all cells(lists)
    #pragma omp parallel for private(ic0,ic1,ic2,item1_p,rho_host,kc,jc,item2_p,d,rv,r2,h,ir2,ir2_h,a,b,rho,rhoDD,F, \
      F_DD,irho,irho_h) shared(sysp) default(none) collapse(3) reduction(+:sum) schedule(static,1)
    for (ic0 = sysp->subDomain.ic_start[0]; ic0 < sysp->subDomain.ic_stop[0]; ic0++)
    for (ic1 = sysp->subDomain.ic_start[1]; ic1 < sysp->subDomain.ic_stop[1]; ic1++)
    for (ic2 = sysp->subDomain.ic_start[2]; ic2 < sysp->subDomain.ic_stop[2]; ic2++)
    {
        // iterate over all items in cell ic
        for (item1_p = sysp->subDomain.grid[ic0][ic1][ic2]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            if (!(sysp->activeGroup == -1 || item1_p->P.groupID == sysp->activeGroup))
                continue;
            rho_host = 0.0;
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
                                for (d=0; d<3; d++)
                                {
                                    if (sysp->ns[d] == 1)
                                    {
                                        if (kc[d]==-1)
                                            rv[d] = item1_p->P.x[d] - item2_p->P.x[d] + sysp->l[d];
                                        else
                                            if (kc[d] == sysp->nc[d])
                                                rv[d] = item1_p->P.x[d] - item2_p->P.x[d] - sysp->l[d];
                                            else
                                                rv[d] = item1_p->P.x[d] - item2_p->P.x[d];
                                    }
                                    else
                                        rv[d] = item1_p->P.x[d] - item2_p->P.x[d];
                                }
                                r2 = SQR(rv[0])+SQR(rv[1])+SQR(rv[2]);
                                if (r2 < r_cutSqd)
                                {
                                    h = sysp->EAM.dr2;
                                    ir2 = (int)(r2 / h);
                                    ir2_h = ir2 + 1;
                                    a = ir2_h - r2/h;
                                    b=1-a;
                                    rho = sysp->EAM.elements[item2_p->P.elementID].rho;
#ifdef USE_CSPLINE
                                    rhoDD = sysp->EAM.elements[item2_p->P.elementID].rhoDD;
                                    rho_host += SPLINE_VAL(a,b,rho,ir2,ir2_h,rhoDD,h);
#else
                                    rho_host += rho[ir2]*a + rho[ir2_h]*b;
#endif
                                }
                            }
                        }
                    }
                }
            }
            h = sysp->EAM.drho;
            irho = (int)(rho_host / h);
            assert(irho < sysp->EAM.Nrho-1);
            irho_h = irho + 1;
            F = sysp->EAM.elements[item1_p->P.elementID].F;
#ifdef USE_CSPLINE
            F_DD = sysp->EAM.elements[item1_p->P.elementID].F_DD;
            a = irho_h - rho_host/h;
            b = 1-a;
            item1_p->FembPrime = SPLINE_DERIV(a,b,F,irho,irho_h,F_DD,h);
            sum += SPLINE_VAL(a,b,F,irho,irho_h,F_DD,h);
#else
            item1_p->FembPrime = (F[irho_h] - F[irho]) / h;
            sum += F[irho] + (rho_host - irho * h) * item1_p->FembPrime;
#endif
        }
    }
    *FembSum_p=sum;
*/
}

static void computeEAM_pass2(fmd_sys_t *sysp, double FembSum)
{
/*
    int jc[3], kc[3];
    int d, ir2, ir2_h;
    TParticleListItem *item1_p, *item2_p;
    const double r_cutSqd = SQR(sysp->EAM.cutoff);
    double r2, rv[3];
    double *rho_i, *rho_j, *phi;
    double *rho_iDD, *rho_jDD, *phiDD;
    double rho_ip, rho_jp;
    double mag;
    double sum[3];
    int element_i, element_j;
    double phi_deriv;
    double a, b, h;
    int ic0, ic1, ic2;
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
    #pragma omp parallel for private(ic0,ic1,ic2,item1_p,d,element_i,rho_i,rho_iDD,kc,jc,item2_p,rv,r2,h,ir2, \
      ir2_h,element_j,phi,phiDD,a,b,phi_deriv,rho_ip,rho_jp,rho_jDD,rho_j,mag,sum) \
      shared(sysp) default(none) collapse(3) reduction(+:potEnergy) schedule(static,1)
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
            element_i = item1_p->P.elementID;
            rho_i = sysp->EAM.elements[element_i].rho;
#ifdef USE_CSPLINE
            rho_iDD = sysp->EAM.elements[element_i].rhoDD;
#endif
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
                                for (d=0; d<3; d++)
                                {
                                    if (sysp->ns[d] == 1)
                                    {
                                        if (kc[d]==-1)
                                            rv[d] = item1_p->P.x[d] - item2_p->P.x[d] + sysp->l[d];
                                        else
                                            if (kc[d] == sysp->nc[d])
                                                rv[d] = item1_p->P.x[d] - item2_p->P.x[d] - sysp->l[d];
                                            else
                                                rv[d] = item1_p->P.x[d] - item2_p->P.x[d];
                                    }
                                    else
                                        rv[d] = item1_p->P.x[d] - item2_p->P.x[d];
                                }
                                r2 = SQR(rv[0])+SQR(rv[1])+SQR(rv[2]);
                                if (r2 < r_cutSqd)
                                {
                                    h = sysp->EAM.dr2;
                                    ir2 = (int)(r2 / h);
                                    ir2_h = ir2 + 1;
                                    element_j = item2_p->P.elementID;
                                    phi = sysp->EAM.elements[element_i].phi[element_j];
#ifdef USE_CSPLINE
                                    phiDD = sysp->EAM.elements[element_i].phiDD[element_j];
                                    a = ir2_h - r2/h;
                                    b = 1-a;
                                    phi_deriv = SPLINE_DERIV(a,b,phi,ir2,ir2_h,phiDD,h);
                                    rho_ip = SPLINE_DERIV(a,b,rho_i,ir2,ir2_h,rho_iDD,h);
                                    if (element_j == element_i)
                                        rho_jp = rho_ip;
                                    else
                                    {
                                        rho_j = sysp->EAM.elements[element_j].rho;
                                        rho_jDD = sysp->EAM.elements[element_j].rhoDD;
                                        rho_jp = SPLINE_DERIV(a,b,rho_j,ir2,ir2_h,rho_jDD,h);
                                    }

                                    mag = 2 * (item1_p->FembPrime * rho_jp +
                                          item2_p->FembPrime * rho_ip + phi_deriv);
                                    potEnergy += SPLINE_VAL(a,b,phi,ir2,ir2_h,phiDD,h);
#else
                                    rho_j = sysp->EAM.elements[element_j].rho;
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
*/
}

static void computeEAM(fmd_sys_t *sysp)
{
    double FembSum;

    fmd_ghostparticles_init(sysp);
    if (sysp->iCompLocOrdParam) compLocOrdParam(sysp);
    computeEAM_pass1(sysp, &FembSum);
    fmd_ghostparticles_update_Femb(sysp);
    computeEAM_pass2(sysp, FembSum);
    fmd_ghostparticles_delete(sysp);
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
                                        for (d=0; d<3; d++)
                                        {
                                            if (sysp->ns[d] == 1)
                                            {
                                                if (kc[d]==-1)
                                                    rv[d] = item1_p->P.x[d] - item2_p->P.x[d] + sysp->l[d];
                                                else
                                                if (kc[d] == sysp->nc[d])
                                                    rv[d] = item1_p->P.x[d] - item2_p->P.x[d] - sysp->l[d];
                                                else
                                                    rv[d] = item1_p->P.x[d] - item2_p->P.x[d];
                                            }
                                            else
                                                rv[d] = item1_p->P.x[d] - item2_p->P.x[d];
                                        }

                                        r2 = SQR(rv[0])+SQR(rv[1])+SQR(rv[2]);
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

void fmd_dync_updateForces(fmd_sys_t *sysp)
{
    if (sysp->potsys.potkinds == NULL)  // just for one time
        fmd_pot_potkinds_update(sysp);

    fmd_ghostparticles_init(sysp);

    if (sysp->potsys.potkinds_num == 1)
    {
        potkind_t potkind = *(potkind_t *)(sysp->potsys.potkinds->data);

        switch (potkind)
        {
            case POTKIND_LJ_6_12:
                computeLJ(sysp);
                break;
        }
    }

    fmd_ghostparticles_delete(sysp);
}
