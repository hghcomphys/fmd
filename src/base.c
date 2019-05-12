/*
  base.c: This file is part of Free Molecular Dynamics

  Copyright (C) 2019 Arham Amouye Foumani

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

#include "base.h"
#include "md_ghost.h"
#include "forces.h"
#ifdef USE_TTM
#include "ttm.h"
#endif
#include <stdarg.h>

const int threeZeros[3] = {0, 0, 0};

static void refreshGrid(fmd_sys_t *sysp, int reverse);

void cleanGridSegment(TCell ***grid, int ic_from[3], int ic_to[3])
{
    int ic[3];
    TParticleListItem *item_p, **item_pp;

    ITERATE(ic, ic_from, ic_to)
    {
        item_pp = &grid[ic[0]][ic[1]][ic[2]];
        item_p = *item_pp;
        while (item_p != NULL)
        {
            removeFromList(item_pp);
            free(item_p);
            item_p = *item_pp;
        }
    }
}

void compLocOrdParam(fmd_sys_t *sysp)
{
/*
    float latticeParameter = sysp->EAM.elements[0].latticeParameter;
    float rCutSqd = SQR(1.32 * latticeParameter);
    int ic[3], jc[3], kc[3];
    int d;
    TParticleListItem *item1_p, *item2_p;
    float arg, real, img;
    int Z;
    float r2, rv[3];
    float q[6][3] = {{1,0,0},{0,1,0},{0,0,1},{1,1,0},{0,1,1},{1,0,1}};
    int i;

    if (sysp->LOPiteration == 0)
        ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
            for (item1_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item1_p != NULL; item1_p = item1_p->next_p)
            {
                item1_p->P.LocOrdParam = 0.;
                for (d=0; d<3; d++)
                    item1_p->P.x_avgd[d] = 0.0;
            }

    (sysp->LOPiteration)++;

    for (i=0; i<6; i++)
        for (d=0; d<3; d++)
            q[i][d] *= 4.0 * M_PI / latticeParameter * q[i][d];

    ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
        for (item1_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            real = img = 0.;
            Z = 0;
            // iterate over neighbor cells of cell ic
            for (kc[0]=ic[0]-1; kc[0]<=ic[0]+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0)
                for (kc[1]=ic[1]-1; kc[1]<=ic[1]+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1)
                    for (kc[2]=ic[2]-1; kc[2]<=ic[2]+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2)
                        // iterate over all items in cell jc
                        for (item2_p = sysp->subDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
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
                                if (r2 < rCutSqd)
                                {
                                    Z++;
                                    for (i=0; i<6; i++)
                                    {
                                        arg = 0.;
                                        for (d=0; d<3; d++)
                                            arg += q[i][d] * rv[d];
                                        real += cosf(arg);
                                        img  += sinf(arg);
                                    }
                                }
                            }
                    }
                }
            }
            item1_p->P.LocOrdParam += (SQR(real)+SQR(img)) / SQR(Z);
            for (d=0; d<3; d++)
                item1_p->P.x_avgd[d] += item1_p->P.x[d];
        }

    if (sysp->LOPiteration != sysp->locOrdParamPeriod) return;

    fmd_ghostparticles_update_LocOrdParam(sysp);

    ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
        for (item1_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item1_p != NULL; item1_p = item1_p->next_p)
        {
            real = item1_p->P.LocOrdParam;
            Z = 0;
            // iterate over neighbor cells of cell ic
            for (kc[0]=ic[0]-1; kc[0]<=ic[0]+1; kc[0]++)
            {
                SET_jc_IN_DIRECTION(0)
                for (kc[1]=ic[1]-1; kc[1]<=ic[1]+1; kc[1]++)
                {
                    SET_jc_IN_DIRECTION(1)
                    for (kc[2]=ic[2]-1; kc[2]<=ic[2]+1; kc[2]++)
                    {
                        SET_jc_IN_DIRECTION(2)
                        // iterate over all items in cell jc
                        for (item2_p = sysp->subDomain.grid[jc[0]][jc[1]][jc[2]]; item2_p != NULL; item2_p = item2_p->next_p)
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
                                if (r2 < rCutSqd)
                                {
                                    Z++;
                                    real += item2_p->P.LocOrdParam;

                                }
                            }
                    }
                }
            }
            item1_p->LocOrdParamAvg = real / ((Z+1) * 36 * sysp->locOrdParamPeriod);
            for (d=0; d<3; d++)
                item1_p->P.x_avgd[d] /= sysp->locOrdParamPeriod;
        }

    sysp->LOPiteration = 0;
*/
}

void fmd_dync_velocityVerlet_takeFirstStep(fmd_sys_t *sysp, int useThermostat)
{
/*
    int ic[3];
    int d;
    TParticleListItem *item_p, **item_pp;
    double velocityScale, mass;
    double x;
    int itemDestroyed;

    if (useThermostat) velocityScale = sqrt(1 + sysp->delta_t / sysp->BerendsenThermostatParam *
                       (sysp->desiredTemperature / sysp->globalTemperature - 1));
    ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
    {
        // iterate over all items in cell ic
        item_pp = &sysp->subDomain.grid[ic[0]][ic[1]][ic[2]];
        item_p = *item_pp;
        while (item_p != NULL)
        {
            if (!(sysp->activeGroup == -1 || item_p->P.groupID == sysp->activeGroup))
            {
                item_pp = &item_p->next_p;
                item_p = *item_pp;
                continue;
            }

            itemDestroyed = 0;
            mass = sysp->EAM.elements[item_p->P.elementID].mass;

            for (d=0; d<3; d++)
            {
                if (sysp->useAutoStep)
                {
                    item_p->P.v_bak[d] = item_p->P.v[d];
                    item_p->P.x_bak[d] = item_p->P.x[d];
                }
                if (useThermostat) item_p->P.v[d] *= velocityScale;
                item_p->P.v[d] += sysp->delta_t * 0.5 / mass * item_p->F[d];
                x = item_p->P.x[d] + sysp->delta_t * item_p->P.v[d];

                if ( (sysp->ns[d] == 1) && ((x < 0.0) || (x >= sysp->l[d])) )
                {
                    if (!sysp->PBC[d])
                    {
                        removeFromList(item_pp);
                        free(item_p);
                        (sysp->subDomain.numberOfParticles)--;
                        itemDestroyed = 1;
                        break;
                    }
                    else
                        if (x < 0.0) x += sysp->l[d]; else x -= sysp->l[d];
                }
                item_p->P.x[d] = x;
            }

            if (!itemDestroyed)
                item_pp = &item_p->next_p;
            item_p = *item_pp;
        }
    }
    refreshGrid(sysp, 0);
*/
}

int fmd_dync_velocityVerlet_takeLastStep(fmd_sys_t *sysp)
{
/*
    int ic[3];
    int d;
    TParticleListItem *item_p;
    double m_vSqd_Sum = 0, m_vSqd_SumSum;
    double mass;
    int returnVal = 0;
    int particlesNum = 0;
    double momentumSum[3] = {0., 0., 0.};

    for (ic[0] = sysp->subDomain.ic_start[0]; ic[0] < sysp->subDomain.ic_stop[0]; ic[0]++)
    {
#ifdef USE_TTM  // for non-reflecting boundary condition
        int reduce = 0;
        double amount;

        if (ttm_useExtended)
        {
            int ic_global0 = ic[0] - sysp->subDomain.ic_start[0] + sysp->subDomain.ic_lower_global[0];
            if (ic_global0 >= ttm_nonrefl_bound_begin)
            {
                int ttm_index;
                double vcm;

                reduce = 1;
                ttm_index = ic[0] - sysp->subDomain.ic_start[0] + 1;
                vcm = ttm_lattice_aux[ttm_index].v_cm[0];
                amount = (ic_global0-ttm_nonrefl_bound_begin) * sysp->cellh[0] / ttm_nonrefl_bound_width;
                if (amount>0.95)
                    amount = vcm*0.99;
                else
                    amount = vcm*delta_t/delta_t_initial*1e-4*(3.*CUBE(amount)+1.5*amount+0.25);
            }
            else
                reduce = 0;
        }
#endif
        for (ic[1] = sysp->subDomain.ic_start[1]; ic[1] < sysp->subDomain.ic_stop[1]; ic[1]++)
            for (ic[2] = sysp->subDomain.ic_start[2]; ic[2] < sysp->subDomain.ic_stop[2]; ic[2]++)
                for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                {
                    if (!(sysp->activeGroup == -1 || item_p->P.groupID == sysp->activeGroup))
                        continue;
                    particlesNum++;
                    mass = sysp->EAM.elements[item_p->P.elementID].mass;
                    for (d=0; d<3; d++)
                    {
                        item_p->P.v[d] += sysp->delta_t * 0.5 / mass * item_p->F[d];
                        momentumSum[d] += mass * item_p->P.v[d];
                    }
#ifdef USE_TTM  // for non-reflecting boundary condition
                    if (reduce) item_p->P.v[0] -= amount;
#endif
                    m_vSqd_Sum += mass * ( SQR(item_p->P.v[0]) +
                                           SQR(item_p->P.v[1]) +
                                           SQR(item_p->P.v[2]) );
                }
    }
    MPI_Reduce(momentumSum, sysp->totalMomentum, 3, MPI_DOUBLE, MPI_SUM,
               MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
    MPI_Allreduce(&particlesNum, &(sysp->activeGroupParticlesNum), 1, MPI_INT, MPI_SUM, sysp->MD_comm);

    if (sysp->activeGroup == -1)
        sysp->totalNoOfParticles = sysp->activeGroupParticlesNum;

    MPI_Allreduce(&m_vSqd_Sum, &m_vSqd_SumSum, 1, MPI_DOUBLE, MPI_SUM, sysp->MD_comm);
    sysp->totalKineticEnergy = 0.5 * m_vSqd_SumSum;
    sysp->totalMDEnergy = sysp->totalKineticEnergy + sysp->totalPotentialEnergy;

    if (sysp->useAutoStep)
    {
        if (sysp->_oldTotalMDEnergy!=0. && fabs((sysp->totalMDEnergy-sysp->_oldTotalMDEnergy)/sysp->_oldTotalMDEnergy) > sysp->autoStepSensitivity)
        {
            if (sysp->_prevFailedMDEnergy!=0. && fabs((sysp->totalMDEnergy-sysp->_prevFailedMDEnergy)/sysp->_prevFailedMDEnergy) < sysp->autoStepSensitivity)
            {  // was jump in energy due to escape of some energetic particle(s) from simulation box?
                if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
                    printf("Maybe the jump was caused by departure of some energetic particle(s). Increasing time step...\n");
                returnVal = 2;
            }
            else
            {
                if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
                {
                    printf("current delta_t = %e\n", sysp->delta_t);
                    printf("Jump in total MD energy (old=%e new=%e)! Decreasing time step...\n", sysp->_oldTotalMDEnergy, sysp->totalMDEnergy);
                }
                sysp->_prevFailedMDEnergy = sysp->totalMDEnergy;
                return 1;
            }
        }
        sysp->_prevFailedMDEnergy = 0.;
        sysp->_oldTotalMDEnergy = sysp->totalMDEnergy;
    }
    sysp->globalTemperature = m_vSqd_SumSum / (3.0 * sysp->activeGroupParticlesNum * K_BOLTZMANN);

    return returnVal;
*/
}

// not correct under periodic boundary conditions
// see [J. Chem. Phys. 131, 154107 (2009)]
double compVirial_internal(fmd_sys_t *sysp)
{
    int ic[3];
    TParticleListItem *item_p;
    double virial = 0.0;
    double virial_global;

    ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
        for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
        {
            virial += item_p->P.x[0] * item_p->F[0] +
                      item_p->P.x[1] * item_p->F[1] +
                      item_p->P.x[2] * item_p->F[2];
        }
    MPI_Reduce(&virial, &virial_global, 1, MPI_DOUBLE, MPI_SUM,
      MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);

    return virial_global;
}

void createCommunicators(fmd_sys_t *sysp)
{
    int mdnum, i;
    MPI_Group world_group, MD_group;
    int *ranks;

    // create MD_comm
    mdnum = sysp->ns[0] * sysp->ns[1] * sysp->ns[2];
    ranks = (int *)malloc(mdnum * sizeof(int));
    for (i=0; i<mdnum; i++)
        ranks[i] = i;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    MPI_Group_incl(world_group, mdnum, ranks, &MD_group);
    MPI_Comm_create(MPI_COMM_WORLD, MD_group, &sysp->MD_comm);
    MPI_Group_free(&world_group);
    MPI_Group_free(&MD_group);
    free(ranks);
}

TCell ***createGrid(int cell_num[3])
{
    TCell ***grid;
    int i, j, k;

    grid = (TCell ***)malloc(cell_num[0]*sizeof(TCell **));
    for (i=0; i < cell_num[0]; i++)
    {
        grid[i] = (TCell **)malloc(cell_num[1]*sizeof(TCell *));
        for (j=0; j < cell_num[1]; j++)
        {
            grid[i][j] = (TCell *)malloc(cell_num[2]*sizeof(TCell));
            for (k=0; k < cell_num[2]; k++)
                grid[i][j][k] = NULL;
        }
    }
    return grid;
}

void findLimits(fmd_sys_t *sysp, double lowerLimit[3], double upperLimit[3])
{
    TParticleListItem *item_p;
    int ic[3];
    int d;
    double localLower[3], localUpper[3];

    localLower[0] = localLower[1] = localLower[2] = DBL_MAX;
    localUpper[0] = localUpper[1] = localUpper[2] = DBL_MIN;

    ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
        for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            for (d=0; d<3; d++)
            {
                if (item_p->P.x[d] < localLower[d])
                    localLower[d] = item_p->P.x[d];
                if (item_p->P.x[d] > localUpper[d])
                    localUpper[d] = item_p->P.x[d];
            }
    MPI_Allreduce(localLower, lowerLimit, 3, MPI_DOUBLE, MPI_MIN, sysp->MD_comm);
    MPI_Allreduce(localUpper, upperLimit, 3, MPI_DOUBLE, MPI_MAX, sysp->MD_comm);
}

void freeGrid(TCell ***grid, int *cell_num)
{
    int i, j;

    cleanGridSegment(grid, threeZeros, cell_num);
    for (i=0; i < cell_num[0]; i++)
    {
        for (j=0; j < cell_num[1]; j++)
            free(grid[i][j]);
        free(grid[i]);
    }
    free(grid);
}

void fmd_subd_free(fmd_sys_t *sysp)
{
    if (sysp->subDomain.grid != NULL)
    {
        freeGrid(sysp->subDomain.grid, sysp->subDomain.cell_num);
        sysp->subDomain.grid = NULL;
    }
}

int getListLength(TParticleListItem *root_p)
{
    int i = 0;

    for ( ; root_p != NULL; root_p = root_p->next_p)
        i++;
    return i;
}

void handleFileOpenError(FILE *fp, char *filename)
{
    if (fp == NULL)
    {
        fprintf(stderr, "ERROR: Unable to open %s!\n", filename);
        MPI_Abort(MPI_COMM_WORLD, ERROR_UNABLE_OPEN_FILE);
    }
}

void identifyProcess(fmd_sys_t *sysp)
{
    int mdnum;

    mdnum = sysp->ns[0] * sysp->ns[1] * sysp->ns[2];
    if (sysp->world_rank < mdnum)
        sysp->isMDprocess = 1;
    else
        sysp->isMDprocess = 0;
#ifdef USE_TTM
    if (ttm_useExtended && sysp->world_rank == mdnum)
        ttm_is_extended_process = 1;
    else
        ttm_is_extended_process = 0;
#endif
}

void fmd_matt_setActiveGroup(fmd_sys_t *sysp, int groupID)
{
    sysp->activeGroup = groupID;
}

void fmd_matt_addVelocity(fmd_sys_t *sysp, int groupID, double vx, double vy, double vz)
{
    TCell ***grid;
    int *start, *stop;
    int ic[3];
    TParticleListItem *item_p;

    if (sysp->particlesDistributed)
    {
        grid = sysp->subDomain.grid;
        start = sysp->subDomain.ic_start;
        stop = sysp->subDomain.ic_stop;
    }
    else
    {
        start = threeZeros;
        if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
        {
            grid = sysp->global_grid;
            stop = sysp->nc;
        }
        else
            stop = start;
    }

    ITERATE(ic, start, stop)
        for (item_p = grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            if (groupID == -1 || groupID == item_p->P.groupID)
            {
                item_p->P.v[0] += vx;
                item_p->P.v[1] += vy;
                item_p->P.v[2] += vz;
            }
}

void fmd_matt_distribute(fmd_sys_t *sysp)
{
/*
    TParticleListItem *item_p;
    int i, k, d, nct, sum_length;
    int ic[3], *ic_length;
    TParticle *is_particles;

    if (sysp->subDomain.grid == NULL) fmd_subd_init(sysp);

    if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
    {
        int r, w;
        int is[3], global_icstart[3], global_icstop[3];
        TParticleListItem **item_pp;
        double m_vSqd_Sum = 0.0;
        double mass;

#ifdef USE_TTM
        ttm_comp_min_atomsNo(global_grid, s_p);
#endif

        ITERATE(ic, threeZeros, sysp->nc)
            for (item_p=sysp->global_grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            {
                mass = sysp->EAM.elements[item_p->P.elementID].mass;
                for (d=0; d<3; d++)
                    m_vSqd_Sum += mass * SQR(item_p->P.v[d]);
            }
        sysp->globalTemperature = m_vSqd_Sum / (3.0 * sysp->totalNoOfParticles * K_BOLTZMANN);

        for (i=0; i < MAINPROCESS(sysp->subDomain.numprocs); i++)
        {
            INVERSEINDEX(i, sysp->ns, is);
            nct = 1;
            for (d=0; d<3; d++)
            {
                r = sysp->nc[d] % sysp->ns[d];
                w = sysp->nc[d] / sysp->ns[d];
                if (is[d] < r)
                {
                    global_icstart[d] = is[d] * (w + 1);
                    global_icstop[d] = global_icstart[d] + w + 1;
                    nct *= w + 1;
                }
                else
                {
                    global_icstart[d] = is[d] * w + r;
                    global_icstop[d] = global_icstart[d] + w;
                    nct *= w;
                }
            }
            ic_length = (int *)malloc((nct+1) * sizeof(int));
            k = sum_length = 0;
            ITERATE(ic, global_icstart, global_icstop)
            {
                ic_length[k] = getListLength(sysp->global_grid[ic[0]][ic[1]][ic[2]]);
                sum_length += ic_length[k++];
            }
            ic_length[nct] = sum_length;
            MPI_Send(ic_length, nct+1, MPI_INT, i, 50, sysp->MD_comm);
            free(ic_length);
            sum_length *= sizeof(TParticle);
            is_particles = (TParticle *)malloc(sum_length);
            k = 0;
            ITERATE(ic, global_icstart, global_icstop)
                for (item_p=sysp->global_grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                    is_particles[k++] = item_p->P;
            MPI_Send(is_particles, sum_length, MPI_CHAR, i, 51, sysp->MD_comm);
            free(is_particles);
        }

        sysp->subDomain.numberOfParticles = 0;
        ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
        {
            item_pp = &sysp->global_grid[ ic[0] - sysp->subDomain.ic_start[0] + sysp->subDomain.ic_global_firstcell[0] ]
                                        [ ic[1] - sysp->subDomain.ic_start[1] + sysp->subDomain.ic_global_firstcell[1] ]
                                        [ ic[2] - sysp->subDomain.ic_start[2] + sysp->subDomain.ic_global_firstcell[2] ];
            item_p = *item_pp;
            while (item_p != NULL)
            {
                removeFromList(item_pp);
                insertInList(&sysp->subDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
                ++(sysp->subDomain.numberOfParticles);
                item_p = *item_pp;
            }
        }
        freeGrid(sysp->global_grid, sysp->nc);
    }
    else
    {
        MPI_Status status;
        int kreceive;

#ifdef USE_TTM
        ttm_comp_min_atomsNo(NULL, s_p);
#endif
        nct = 1;
        for (d=0; d<3; d++)
            nct *= sysp->subDomain.ic_stop[d] - sysp->subDomain.ic_start[d];
        ic_length = (int *)malloc((nct+1) * sizeof(int));
        MPI_Recv(ic_length, nct+1, MPI_INT, MAINPROCESS(sysp->subDomain.numprocs),
                 50, sysp->MD_comm, &status);
        sysp->subDomain.numberOfParticles = sum_length = ic_length[nct];
        sum_length *= sizeof(TParticle);
        is_particles = (TParticle *)malloc(sum_length);
        MPI_Recv(is_particles, sum_length, MPI_CHAR,
                 MAINPROCESS(sysp->subDomain.numprocs), 51, sysp->MD_comm, &status);
        kreceive = k = 0;
        ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
        {
            for (i=0; i<ic_length[kreceive]; i++)
            {
                item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
                item_p->P = is_particles[k++];
                insertInList(&sysp->subDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
            }
            kreceive++;
        }
        free(ic_length);
        free(is_particles);
    }

    MPI_Bcast(&sysp->totalNoOfParticles, 1, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs),
              sysp->MD_comm);
    MPI_Bcast(&sysp->globalTemperature, 1, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs),
              sysp->MD_comm);

    sysp->totalKineticEnergy = 3.0/2.0 * sysp->totalNoOfParticles * K_BOLTZMANN * sysp->globalTemperature;
    sysp->globalGridExists = 0;
    sysp->particlesDistributed = 1;
*/
}

void fmd_subd_init(fmd_sys_t *sysp)
{
    int d;

    // initialize is
    INVERSEINDEX(sysp->subDomain.myrank, sysp->ns, sysp->subDomain.is);
    // initialize rank_of_lower_subd and rank_of_upper_subd (neighbor processes)
    int istemp[3];
    for (d=0; d<3; d++)
        istemp[d] = sysp->subDomain.is[d];
    for (d=0; d<3; d++)
    {
        istemp[d] = (sysp->subDomain.is[d] - 1 + sysp->ns[d]) % sysp->ns[d];
        sysp->subDomain.rank_of_lower_subd[d] = INDEX(istemp, sysp->ns);
        istemp[d] = (sysp->subDomain.is[d] + 1) % sysp->ns[d];
        sysp->subDomain.rank_of_upper_subd[d] = INDEX(istemp, sysp->ns);
        istemp[d] = sysp->subDomain.is[d];
    }
    //
    for (d=0; d<3; d++)
    {
        int r, w;

        if (sysp->ns[d] == 1) sysp->subDomain.ic_start[d] = 0; else sysp->subDomain.ic_start[d] = 1;
        r = sysp->nc[d] % sysp->ns[d];
        w = sysp->nc[d] / sysp->ns[d];
        if (sysp->subDomain.is[d] < r)
        {
            sysp->subDomain.ic_stop[d] = sysp->subDomain.ic_start[d] + w + 1;
            sysp->subDomain.ic_global_firstcell[d] = sysp->subDomain.is[d] * (w + 1);
        }
        else
        {
            sysp->subDomain.ic_stop[d] = sysp->subDomain.ic_start[d] + w;
            sysp->subDomain.ic_global_firstcell[d] = sysp->subDomain.is[d] * w + r;
        }
        sysp->subDomain.cell_num[d] = sysp->subDomain.ic_stop[d] + sysp->subDomain.ic_start[d];
    }

    sysp->subDomain.grid = createGrid(sysp->subDomain.cell_num);
}

void insertInList(TParticleListItem **root_pp, TParticleListItem *item_p)
{
    item_p->next_p = *root_pp;
    *root_pp = item_p;
}

void fmd_io_loadState(fmd_sys_t *sysp, char *file, int useTime)
{
/*
    TParticleListItem *item_p;
    FILE *fp;
    char name[3];
    int i, j, d;
    int ic[3];
    double stateFileTime;
    int particlesNum;
    double l0, l1, l2;
    int PBC0, PBC1, PBC2;

    if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
    {
        fp = fopen(file, "r");
        handleFileOpenError(fp, file);
        fscanf(fp, "%lf", &stateFileTime);
        if (useTime)
            sysp->mdTime = stateFileTime;
        fscanf(fp, "%d\n", &particlesNum);
        sysp->totalNoOfParticles += particlesNum;
        fscanf(fp, "%lf%lf%lf", &l0, &l1, &l2);
        fscanf(fp, "%d%d%d", &PBC0, &PBC1, &PBC2);
    }

    if (!sysp->boxSizeDetermined)
    {
        if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
        {
            sysp->l[0] = l0;
            sysp->l[1] = l1;
            sysp->l[2] = l2;
        }
        MPI_Bcast(&sysp->l, 3, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        sysp->boxSizeDetermined = 1;
    }

    if (!sysp->PBCdetermined)
    {
        if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
        {
            sysp->PBC[0] = PBC0;
            sysp->PBC[1] = PBC1;
            sysp->PBC[2] = PBC2;
        }
        MPI_Bcast(&sysp->PBC, 3, MPI_INT, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        sysp->PBCdetermined = 1;
    }

    if (!sysp->globalGridExists)
        fmd_box_createGrid(sysp, sysp->cutoffRadius);

    if (useTime)
        MPI_Bcast(&sysp->mdTime, 1, MPI_DOUBLE, MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);

    if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
    {
        for (i=0; i < particlesNum; i++)
        {
            item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
            fscanf(fp, "%s%d", name, &item_p->P.groupID);
            for (j=0; j < sysp->EAM.elementsNo; j++)
                if (strcmp(name, sysp->EAM.elements[j].name) == 0)
                {
                    item_p->P.elementID = j;
                    break;
                }
            fscanf(fp, "%lf%lf%lf", &item_p->P.x[0], &item_p->P.x[1], &item_p->P.x[2]);
            fscanf(fp, "%lf%lf%lf", &item_p->P.v[0], &item_p->P.v[1], &item_p->P.v[2]);
            for (d=0; d<3; d++)
                ic[d] = (int)floor(item_p->P.x[d] / sysp->cellh[d]);
            insertInList(&(sysp->global_grid[ic[0]][ic[1]][ic[2]]), item_p);
        }
        fclose(fp);
    }
*/
}

static void refreshGrid(fmd_sys_t *sysp, int reverse)
{
    int ic[3], jc[3];
    int d;
    TParticleListItem **item_pp;
    TParticleListItem *item_p;

    // iterate over all cells(lists)
    ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
    {
        // iterate over all items in cell ic
        item_pp = &sysp->subDomain.grid[ic[0]][ic[1]][ic[2]];
        item_p = *item_pp;
        while (item_p != NULL)
        {
            for (d=0; d<3; d++)
            {
                jc[d] = (int)floor(item_p->P.x[d] / sysp->cellh[d]) - sysp->subDomain.ic_global_firstcell[d] + sysp->subDomain.ic_start[d];
                if (jc[d] < 0)
                {
                    if (reverse && sysp->PBC[d] && sysp->ns[d] > 1 &&
                        sysp->subDomain.is[d]==sysp->ns[d]-1 && ic[d] >= sysp->subDomain.ic_stop[d]-sysp->subDomain.ic_start[d])
                    {
                        item_p->P.x[d] += sysp->l[d];
                        jc[d] = ic[d] + sysp->subDomain.ic_start[d];
                    }
                    else
                    {
                        fprintf(stderr,"ERROR: Unexpected particle position!\n");
                        MPI_Abort(MPI_COMM_WORLD, ERROR_UNEXPECTED_PARTICLE_POSITION);
                    }
                }
                else
                    if (jc[d] >= sysp->subDomain.cell_num[d])
                    {
                        if (reverse && sysp->PBC[d] && sysp->ns[d] > 1 &&
                            sysp->subDomain.is[d]==0 && ic[d] < 2*sysp->subDomain.ic_start[d])
                        {
                            item_p->P.x[d] -= sysp->l[d];
                            jc[d] = ic[d] - sysp->subDomain.ic_start[d];
                        }
                        else
                        {
                            fprintf(stderr, "ERROR: Unexpected particle position!\n");
                            MPI_Abort(MPI_COMM_WORLD, ERROR_UNEXPECTED_PARTICLE_POSITION);
                        }
                    }
            }

            if ((ic[0] != jc[0]) || (ic[1] != jc[1]) || (ic[2] != jc[2]))
            {
                removeFromList(item_pp);
                insertInList(&sysp->subDomain.grid[jc[0]][jc[1]][jc[2]], item_p);
            }
            else
                item_pp = &item_p->next_p;

            item_p = *item_pp;
        }
    }

    // now particles in ghost cells migrate to neighbour subdomains
    fmd_particles_migrate(sysp);
}

void removeFromList(TParticleListItem **item_pp)
{
    *item_pp = (*item_pp)->next_p;
}

void rescaleVelocities(fmd_sys_t *sysp)
{
    int ic[3];
    int d;
    TParticleListItem *item_p;
    double scale;

    scale = sqrt(sysp->desiredTemperature / sysp->globalTemperature);

    ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
        for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            for (d=0; d<3; d++)
                item_p->P.v[d] *= scale;
    sysp->globalTemperature = sysp->desiredTemperature;
}

void restoreBackups(fmd_sys_t *sysp)
{
    int ic[3], d;
    TParticleListItem *item_p;

    ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
        for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
        {
            for (d=0; d<3; d++)
            {
                item_p->P.x[d] = item_p->P.x_bak[d];
                item_p->P.v[d] = item_p->P.v_bak[d];
            }
        }
    refreshGrid(sysp, 1);

#ifdef USE_TTM
    ttm_restoreBackups();
#endif
}

void fmd_matt_saveConfiguration(fmd_sys_t *sysp)
{
/*
    int ic[3];
    TParticleListItem *item_p;
    TXYZ_Struct *localData, *globalData;
    int *nums, *recvcounts, *displs;
    int k;

    localData = (TXYZ_Struct *)malloc(sysp->subDomain.numberOfParticles * sizeof(TXYZ_Struct));
    k=0;
    ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
        for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
        {
            if (sysp->iCompLocOrdParam)
            {
                localData[k].x[0] = (float)item_p->P.x_avgd[0];
                localData[k].x[1] = (float)item_p->P.x_avgd[1];
                localData[k].x[2] = (float)item_p->P.x_avgd[2];
            }
            else
            {
                localData[k].x[0] = (float)item_p->P.x[0];
                localData[k].x[1] = (float)item_p->P.x[1];
                localData[k].x[2] = (float)item_p->P.x[2];
            }
            localData[k].var = item_p->P.groupID;
            localData[k].elementID = item_p->P.elementID;
            k++;
        }

    nums = (int *)malloc(sysp->subDomain.numprocs * sizeof(int));
    MPI_Allgather(&sysp->subDomain.numberOfParticles, 1, MPI_INT, nums, 1, MPI_INT,
        sysp->MD_comm);
    sysp->totalNoOfParticles = 0;
    for (k=0; k < sysp->subDomain.numprocs; k++)
        sysp->totalNoOfParticles += nums[k];

    if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
    {
        int displ = 0;

        globalData = (TXYZ_Struct *)malloc(sysp->totalNoOfParticles * sizeof(TXYZ_Struct));
        recvcounts = (int *)malloc(sysp->subDomain.numprocs * sizeof(int));
        displs     = (int *)malloc(sysp->subDomain.numprocs * sizeof(int));
        for (k=0; k < sysp->subDomain.numprocs; k++)
        {
            recvcounts[k] = nums[k] * sizeof(TXYZ_Struct);
            displs[k] = displ;
            displ += recvcounts[k];
        }
    }
    free(nums);
    MPI_Gatherv(localData, sysp->subDomain.numberOfParticles * sizeof(TXYZ_Struct), MPI_CHAR,
        globalData, recvcounts, displs, MPI_CHAR, MAINPROCESS(sysp->subDomain.numprocs),
        sysp->MD_comm);
    free(localData);

    if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
    {
        char configPath[MAX_PATH_LENGTH];
        char *elementName;
        int i;

        free(recvcounts);
        free(displs);

        switch (sysp->saveConfigMode)
        {
            case scmXYZParticlesNum:
                if (sysp->totalNoOfParticles != sysp->_oldNumberOfParticles)
                {
                    if (sysp->_oldNumberOfParticles != -1) fclose(sysp->configFilep);
                    sprintf(configPath, "%s%d.xyz", sysp->saveDirectory, sysp->totalNoOfParticles);
                    sysp->configFilep = fopen(configPath, "w");
                    handleFileOpenError(sysp->configFilep, configPath);
                    sysp->_oldNumberOfParticles = sysp->totalNoOfParticles;
                }
                break;

            case scmXYZSeparate:
                sprintf(configPath, "%s%05d.xyz", sysp->saveDirectory, sysp->_fileIndex++);
                sysp->configFilep = fopen(configPath, "w");
                handleFileOpenError(sysp->configFilep, configPath);
                break;

            case scmCSV:
                sprintf(configPath, "%s%05d.csv", sysp->saveDirectory, sysp->_fileIndex++);
                sysp->configFilep = fopen(configPath, "w");
                handleFileOpenError(sysp->configFilep, configPath);
                for (i=0; i < sysp->EAM.elementsNo; i++)
                {
                    for (k=0; k < sysp->totalNoOfParticles; k++)
                        if (globalData[k].elementID == i)
                            fprintf(sysp->configFilep, "%.2f, %.2f, %.2f, %d, %.4f\n",
                                globalData[k].x[0], globalData[k].x[1],
                                globalData[k].x[2], i, globalData[k].var);
                }
                break;

            case scmVTF:
            {
                int atomID = 0;
                sprintf(configPath, "%s%05d.vtf", sysp->saveDirectory, sysp->_fileIndex++);
                sysp->configFilep = fopen(configPath, "w");
                handleFileOpenError(sysp->configFilep, configPath);
                for (i=0; i < sysp->EAM.elementsNo; i++)
                {
                    int count = 0;
                    for (k=0; k < sysp->totalNoOfParticles; k++)
                        if (globalData[k].elementID == i)
                            count++;
                    fprintf(sysp->configFilep, "atom %d:%d\tname %s\n", atomID,
                     atomID+count-1, sysp->EAM.elements[i].name);
                    atomID += count;
                }
                atomID = 0;
                for (i=0; i < sysp->EAM.elementsNo; i++)
                {
                    for (k=0; k < sysp->totalNoOfParticles; k++)
                        if (globalData[k].elementID == i)
                        {
                            fprintf(sysp->configFilep, "%d\tbeta %.4f\n", atomID,
                             globalData[k].var);
                            atomID++;
                        }
                }
                fprintf(sysp->configFilep, "timestep\n");
                for (i=0; i < sysp->EAM.elementsNo; i++)
                {
                    for (k=0; k < sysp->totalNoOfParticles; k++)
                        if (globalData[k].elementID == i)
                            fprintf(sysp->configFilep, "%.2f\t%.2f\t%.2f\n",
                                globalData[k].x[0], globalData[k].x[1],
                                globalData[k].x[2]);
                }
                break;
            }
        }

        if (sysp->saveConfigMode == scmXYZSeparate || sysp->saveConfigMode == scmXYZParticlesNum)
        {
            fprintf(sysp->configFilep, "%d\n\n", sysp->totalNoOfParticles);
            for (i=0; i < sysp->EAM.elementsNo; i++)
            {
                elementName = sysp->EAM.elements[i].name;
                for (k=0; k < sysp->totalNoOfParticles; k++)
                    if (globalData[k].elementID == i)
                        fprintf(sysp->configFilep, "%s\t%.2f\t%.2f\t%.2f\n",
                            elementName, globalData[k].x[0], globalData[k].x[1],
                            globalData[k].x[2]);
            }
        }
        free(globalData);

        if (sysp->saveConfigMode == scmXYZSeparate || sysp->saveConfigMode == scmCSV ||
         sysp->saveConfigMode == scmVTF)
            fclose(sysp->configFilep);
    }
*/
}

void fmd_io_saveState(fmd_sys_t *sysp, char *filename)
{
/*
    TParticle *is_particles, *P_p;
    int *nums;
    int i, k, ic[3];
    TParticleListItem *item_p;
    char stateFilePath[MAX_PATH_LENGTH];
    FILE *fp;
    MPI_Status status;

    if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
    {
        nums = (int *)malloc(sysp->subDomain.numprocs * sizeof(int));
        MPI_Gather(&sysp->subDomain.numberOfParticles, 1, MPI_INT, nums, 1, MPI_INT,
            MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        sysp->totalNoOfParticles = 0;
        for (k=0; k < sysp->subDomain.numprocs; k++)
            sysp->totalNoOfParticles += nums[k];
        sprintf(stateFilePath, "%s%s", sysp->saveDirectory, filename);
        fp = fopen(stateFilePath, "w");
        handleFileOpenError(fp, stateFilePath);
        fprintf(fp, "%.16e\n", sysp->mdTime);
        fprintf(fp, "%d\n", sysp->totalNoOfParticles);
        fprintf(fp, "%.16e\t%.16e\t%.16e\n", sysp->l[0], sysp->l[1], sysp->l[2]);
        fprintf(fp, "%d %d %d\n", sysp->PBC[0], sysp->PBC[1], sysp->PBC[2]);
        ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            {
                fprintf(fp, "%s %d\n", sysp->EAM.elements[item_p->P.elementID].name, item_p->P.groupID);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", item_p->P.x[0], item_p->P.x[1], item_p->P.x[2]);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", item_p->P.v[0], item_p->P.v[1], item_p->P.v[2]);
            }
        for (i=0; i < MAINPROCESS(sysp->subDomain.numprocs); i++)
        {
            is_particles = (TParticle *)malloc(nums[i] * sizeof(TParticle));
            MPI_Recv(is_particles, nums[i] * sizeof(TParticle), MPI_CHAR, i, 150,
                sysp->MD_comm, &status);
            for (k=0; k < nums[i]; k++)
            {
                P_p = is_particles + k;
                fprintf(fp, "%s %d\n", sysp->EAM.elements[P_p->elementID].name, P_p->groupID);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", P_p->x[0], P_p->x[1], P_p->x[2]);
                fprintf(fp, "%.16e\t%.16e\t%.16e\n", P_p->v[0], P_p->v[1], P_p->v[2]);
            }
            free(is_particles);
        }
        free(nums);
        fclose(fp);
    }
    else
    {
        MPI_Gather(&sysp->subDomain.numberOfParticles, 1, MPI_INT, nums, 1, MPI_INT,
            MAINPROCESS(sysp->subDomain.numprocs), sysp->MD_comm);
        is_particles = (TParticle *)malloc(sysp->subDomain.numberOfParticles * sizeof(TParticle));
        k = 0;
        ITERATE(ic, sysp->subDomain.ic_start, sysp->subDomain.ic_stop)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                is_particles[k++] = item_p->P;
        MPI_Send(is_particles, sysp->subDomain.numberOfParticles * sizeof(TParticle), MPI_CHAR,
            MAINPROCESS(sysp->subDomain.numprocs), 150, sysp->MD_comm);
        free(is_particles);
    }
*/
}

void fmd_matt_setDesiredTemperature(fmd_sys_t *sysp, double desiredTemperature)
{
    sysp->desiredTemperature = desiredTemperature;
}

void fmd_box_setPBC(fmd_sys_t *sysp, int PBCx, int PBCy, int PBCz)
{
    sysp->PBC[0] = PBCx;
    sysp->PBC[1] = PBCy;
    sysp->PBC[2] = PBCz;
    sysp->PBCdetermined = 1;
}

void fmd_box_setSubDomains(fmd_sys_t *sysp, int dimx, int dimy, int dimz)
{
    sysp->ns[0] = dimx;
    sysp->ns[1] = dimy;
    sysp->ns[2] = dimz;
    identifyProcess(sysp);
    createCommunicators(sysp);
    if (sysp->isMDprocess)
    {
        MPI_Comm_size(sysp->MD_comm, &sysp->subDomain.numprocs);
        MPI_Comm_rank(sysp->MD_comm, &sysp->subDomain.myrank);
    }
}

fmd_sys_t *fmd_sys_create()
{
    fmd_sys_t *sysp = (fmd_sys_t *)malloc(sizeof(fmd_sys_t));

    int isMPIInitialized;
    MPI_Initialized(&isMPIInitialized);
    if (!isMPIInitialized) MPI_Init(NULL, NULL);

    omp_set_num_threads(1);

    MPI_Comm_size(MPI_COMM_WORLD, &(sysp->world_numprocs));
    MPI_Comm_rank(MPI_COMM_WORLD, &(sysp->world_rank));
    sysp->LOPiteration = 0;
    sysp->useAutoStep = 0;
    sysp->mdTime = 0.0;
    sysp->saveDirectory[0] = '\0';
    sysp->iCompLocOrdParam = 0;
    sysp->subDomain.grid = NULL;
    sysp->totalNoOfParticles = 0;
    sysp->activeGroup = -1;             // all groups are active by default
    sysp->particlesDistributed = 0;
    sysp->globalGridExists = 0;
    sysp->boxSizeDetermined = 0;
    sysp->PBCdetermined = 0;
    sysp->_oldNumberOfParticles = -1;
    sysp->_fileIndex = 0;
    sysp->_oldTotalMDEnergy = 0.0;
    sysp->_prevFailedMDEnergy = 0.0;
    fmd_pot_init(sysp);

    // this must be the last statement before return
    sysp->wallTimeOrigin = MPI_Wtime();

    return sysp;
}

void fmd_box_setSize(fmd_sys_t *sysp, double sx, double sy, double sz)
{
    if (!sysp->globalGridExists)
    {
        sysp->l[0] = sx;
        sysp->l[1] = sy;
        sysp->l[2] = sz;
        sysp->boxSizeDetermined = 1;
    }
}

double fmd_proc_getWallTime(fmd_sys_t *sysp)
{
    return (MPI_Wtime() - sysp->wallTimeOrigin);
}

int fmd_proc_isMD(fmd_sys_t *sysp)
{
    return sysp->isMDprocess;
}

void fmd_box_createGrid(fmd_sys_t *sysp, double cutoff)
{
    int d;

    for (d=0; d<3; d++)
    {
        sysp->nc[d] = (int)(sysp->l[d] / cutoff);
        if ((sysp->nc[d] < 3) && sysp->PBC[d])
        {
            fprintf(stderr, "ERROR: nc[%d] = %d. Under PBC, this must be greater than 2!\n", d, sysp->nc[d]);
            MPI_Abort(MPI_COMM_WORLD, ERROR_NC_TOO_SMALL);
        }
        sysp->cellh[d] = sysp->l[d] / sysp->nc[d];
    }

    if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
        sysp->global_grid = createGrid(sysp->nc);
    sysp->globalGridExists = 1;
    sysp->cutoffRadius = cutoff;
}

void fmd_matt_makeCuboidFCC(fmd_sys_t *sysp, double x, double y, double z,
  int dimx, int dimy, int dimz, double latticeParameter, int elementID, int groupID)
{
/*
    if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
    {
        double mass, stdDevVelocity;
        int crystalCell[3], ic[3];
        double rFCC[4][3] = {{0.0, 0.0, 0.0}, {0.0, 0.5, 0.5},
                             {0.5, 0.0, 0.5}, {0.5, 0.5, 0.0}};
        double momentumSum[3] = {0.0, 0.0, 0.0};
        TParticleListItem *item_p;
        int dims[3] = {dimx, dimy, dimz};
        double r0[3] = {x, y, z};
        int i, d;

        gsl_rng *random_fast;
        random_fast = gsl_rng_alloc(gsl_rng_rand);
        gsl_rng_set(random_fast, time(NULL));

        int atomsNum = 4 * dimx * dimy * dimz;
        mass = sysp->EAM.elements[elementID].mass;
        stdDevVelocity = sqrt(K_BOLTZMANN * sysp->desiredTemperature / mass);
        ITERATE(crystalCell, threeZeros, dims)
            for (i=0; i<4; i++)
            {
                item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
                item_p->P.elementID = elementID;
                item_p->P.groupID = -2;
                for (d=0; d<3; d++)
                {
                    item_p->P.x[d] = r0[d] + (crystalCell[d] + .25 + rFCC[i][d]) * latticeParameter;
                    item_p->P.v[d] = gsl_ran_gaussian_ziggurat(random_fast, stdDevVelocity);
                    momentumSum[d] += mass * item_p->P.v[d];
                    ic[d] = (int)floor(item_p->P.x[d] / sysp->cellh[d]);
                }
                insertInList(&sysp->global_grid[ic[0]][ic[1]][ic[2]], item_p);
            }

        gsl_rng_free(random_fast);

        ITERATE(ic, threeZeros, sysp->nc)
            for (item_p=sysp->global_grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            {
                if (item_p->P.groupID == -2)
                {
                    item_p->P.groupID = groupID;
                    for (d=0; d<3; d++)
                        item_p->P.v[d] -= momentumSum[d] / (atomsNum * mass);
                }
            }

        sysp->totalNoOfParticles += atomsNum;
    }
*/
}

void fmd_io_setSaveDirectory(fmd_sys_t *sysp, char *directory)
{
    strcpy(sysp->saveDirectory, directory);
}

void fmd_io_setSaveConfigMode(fmd_sys_t *sysp, fmd_SaveConfigMode_t mode)
{
    sysp->saveConfigMode = mode;
}

void fmd_dync_setTimeStep(fmd_sys_t *sysp, double timeStep)
{
    sysp->delta_t = timeStep;
}

double fmd_dync_getTimeStep(fmd_sys_t *sysp)
{
    return sysp->delta_t;
}

double fmd_dync_getTime(fmd_sys_t *sysp)
{
    return sysp->mdTime;
}

void fmd_dync_incTime(fmd_sys_t *sysp)
{
    sysp->mdTime += sysp->delta_t;
}

void fmd_dync_equilibrate(fmd_sys_t *sysp, int groupID, double duration, double strength)
{
    sysp->mdTime = 0.0;
    sysp->globalTemperature = sysp->desiredTemperature;
    fmd_dync_setBerendsenThermostatParameter(sysp, strength);
    fmd_matt_setActiveGroup(sysp, groupID);

    // compute forces for the first time
    fmd_dync_updateForces(sysp);

    while (sysp->mdTime < duration)
    {
        // take first step of velocity Verlet integrator
        fmd_dync_velocityVerlet_takeFirstStep(sysp, 1);

        // compute forces
        fmd_dync_updateForces(sysp);

        // take last step of velocity Verlet integrator
        fmd_dync_velocityVerlet_takeLastStep(sysp);

        sysp->mdTime += sysp->delta_t;
    }
    // end of the time loop
    sysp->mdTime = 0.0;
}

void fmd_io_printf(fmd_sys_t *sysp, const char * restrict format, ...)
{
    if (sysp->isMDprocess && sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
    {
        va_list argptr;

        va_start(argptr, format);
        vprintf(format, argptr);
        va_end(argptr);
    }
}

double fmd_matt_getTotalEnergy(fmd_sys_t *sysp)
{
    return sysp->totalKineticEnergy + sysp->totalPotentialEnergy;
}

void fmd_matt_giveTemperature(fmd_sys_t *sysp, int groupID)
{
/*
    TCell ***grid;
    int *start, *stop;
    int ic[3];
    TParticleListItem *item_p;

    if (sysp->particlesDistributed)
    {
        grid = sysp->subDomain.grid;
        start = sysp->subDomain.ic_start;
        stop = sysp->subDomain.ic_stop;
    }
    else
    {
        start = threeZeros;
        if (sysp->subDomain.myrank == MAINPROCESS(sysp->subDomain.numprocs))
        {
            grid = sysp->global_grid;
            stop = sysp->nc;
        }
        else
            stop = start;
    }

    gsl_rng *rng;
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng, time(NULL));
    int d;

    ITERATE(ic, start, stop)
        for (item_p = grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
        {
            if (groupID == -1 || groupID == item_p->P.groupID)
            {
                double mass = sysp->EAM.elements[item_p->P.elementID].mass;
                double stdDevVelocity = sqrt(K_BOLTZMANN * sysp->desiredTemperature / mass);

                for (d=0; d<3; d++)
                    item_p->P.v[d] = gsl_ran_gaussian_ziggurat(rng, stdDevVelocity);
            }
        }

    gsl_rng_free(rng);
*/
}

double fmd_matt_getGlobalTemperature(fmd_sys_t *sysp)
{
    return sysp->globalTemperature;
}

void fmd_dync_setBerendsenThermostatParameter(fmd_sys_t *sysp, double parameter)
{
    sysp->BerendsenThermostatParam = parameter;
}

void fmd_sys_free(fmd_sys_t *sysp, int finalizeMPI)
{
    fmd_subd_free(sysp);
    fmd_pot_free(sysp);
    free(sysp);
    if (finalizeMPI)
    {
        int isMPIFinalized;
        MPI_Finalized(&isMPIFinalized);
        if (!isMPIFinalized) MPI_Finalize();
    }
}
