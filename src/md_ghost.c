/*
  md_ghost.c: This file is part of Free Molecular Dynamics

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

// functions for working with ghost cells of MD subdomains, including
// communication of data from or to them

#include "base.h"
#include "md_ghost.h"

static void particles_migrate_in_direction_d(
    fmd_sys_t *sysp, int d, int ic_start_send_lower[3],
    int ic_stop_send_lower[3], int ic_start_receive_lower[3],
    int ic_stop_receive_lower[3], int ic_start_send_upper[3], int ic_stop_send_upper[3],
    int ic_start_receive_upper[3], int ic_stop_receive_upper[3])
{
    MPI_Status status;
    MPI_Request request;
    int sum_length_send, sum_length_receive;
    int k, kreceive, cells_num;
    int *cells_length_send, *cells_length_receive, ic[3];
    TParticle *particles_send, *particles_receive;
    int dd;
    TParticleListItem *item_p, **item_pp;
    int cc;

    if ( ((sysp->subDomain.is[d] == 0) || (sysp->subDomain.is[d] == sysp->ns[d]-1)) && !sysp->PBC[d] )
    {
        if (sysp->subDomain.is[d] == 0)
        {
            // receiving from upper process
            cells_num = 1;
            for (dd=0; dd<3; dd++)
                cells_num *= ic_stop_receive_upper[dd] - ic_start_receive_upper[dd];
            cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
            MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d],
                     1, sysp->MD_comm, &status);
            sysp->subDomain.numberOfParticles += cells_length_receive[cells_num];
            sum_length_receive = cells_length_receive[cells_num] * sizeof(TParticle);
            particles_receive = (TParticle *)malloc(sum_length_receive);
            MPI_Recv(particles_receive, sum_length_receive, MPI_CHAR,
                     sysp->subDomain.rank_of_upper_subd[d], 2, sysp->MD_comm, &status);
            kreceive = k = 0;
            ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
            {
                for (cc=0; cc<cells_length_receive[kreceive]; cc++)
                {
                    item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
                    item_p->P = particles_receive[k++];
                    insertInList(&sysp->subDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
                }
                kreceive++;
            }
            free(cells_length_receive);
            free(particles_receive);
            //
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            {
                item_pp = &sysp->subDomain.grid[ic[0]][ic[1]][ic[2]];
                item_p = *item_pp;
                while (item_p != NULL)
                {
                    removeFromList(item_pp);
                    --(sysp->subDomain.numberOfParticles);
                    free(item_p);
                    item_p = *item_pp;
                }
            }
            // sending to upper process
            cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
            k = sum_length_send = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            {
                cells_length_send[k] = getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
                sum_length_send += cells_length_send[k++];
            }
            cells_length_send[cells_num] = sum_length_send;
            MPI_Send(cells_length_send, cells_num+1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 3,
                     sysp->MD_comm);
            free(cells_length_send);
            sysp->subDomain.numberOfParticles -= sum_length_send;
            sum_length_send *= sizeof(TParticle);
            particles_send = (TParticle *)malloc(sum_length_send);
            k = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            {
                item_pp = &sysp->subDomain.grid[ic[0]][ic[1]][ic[2]];
                item_p = *item_pp;
                while (item_p != NULL)
                {
                    removeFromList(item_pp);
                    particles_send[k++] = item_p->P;
                    free(item_p);
                    item_p = *item_pp;

                }
            }
            MPI_Send(particles_send, sum_length_send, MPI_CHAR,
                     sysp->subDomain.rank_of_upper_subd[d], 4, sysp->MD_comm);
            free(particles_send);
        }
        else // if (sysp->subDomain.is[d] == sysp->ns[d]-1)
        {
            // sending to lower process
            cells_num = 1;
            for (dd=0; dd<3; dd++)
                cells_num *= ic_stop_send_lower[dd] - ic_start_send_lower[dd];
            cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
            k = sum_length_send = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            {
                cells_length_send[k] = getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
                sum_length_send += cells_length_send[k++];
            }
            cells_length_send[cells_num] = sum_length_send;
            MPI_Send(cells_length_send, cells_num+1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 1,
                     sysp->MD_comm);
            free(cells_length_send);
            sysp->subDomain.numberOfParticles -= sum_length_send;
            sum_length_send *= sizeof(TParticle);
            particles_send = (TParticle *)malloc(sum_length_send);
            k = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            {
                item_pp = &sysp->subDomain.grid[ic[0]][ic[1]][ic[2]];
                item_p = *item_pp;
                while (item_p != NULL)
                {
                    removeFromList(item_pp);
                    particles_send[k++] = item_p->P;
                    free(item_p);
                    item_p = *item_pp;
                }
            }
            MPI_Send(particles_send, sum_length_send, MPI_CHAR,
                     sysp->subDomain.rank_of_lower_subd[d], 2, sysp->MD_comm);
            free(particles_send);
            // receiving from lower process
            cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
            MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 3,
                     sysp->MD_comm, &status);
            sysp->subDomain.numberOfParticles += cells_length_receive[cells_num];
            sum_length_receive = cells_length_receive[cells_num] * sizeof(TParticle);
            particles_receive = (TParticle *)malloc(sum_length_receive);
            MPI_Recv(particles_receive, sum_length_receive, MPI_CHAR,
                     sysp->subDomain.rank_of_lower_subd[d], 4, sysp->MD_comm, &status);
            kreceive = k = 0;
            ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
            {
                for (cc=0; cc<cells_length_receive[kreceive]; cc++)
                {
                    item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
                    item_p->P = particles_receive[k++];
                    insertInList(&sysp->subDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
                }
                kreceive++;
            }
            free(cells_length_receive);
            free(particles_receive);
            //
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            {
                item_pp = &sysp->subDomain.grid[ic[0]][ic[1]][ic[2]];
                item_p = *item_pp;
                while (item_p != NULL)
                {
                    removeFromList(item_pp);
                    --(sysp->subDomain.numberOfParticles);
                    free(item_p);
                    item_p = *item_pp;
                }
            }
        }
    }
    else
    {
        // sending to lower process, receiving from upper process
        cells_num = 1;
        for (dd=0; dd<3; dd++)
            cells_num *= ic_stop_send_lower[dd] - ic_start_send_lower[dd];
        cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
        cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
        k = sum_length_send = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
        {
            cells_length_send[k] = getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
            sum_length_send += cells_length_send[k++];
        }
        cells_length_send[cells_num] = sum_length_send;
        MPI_Isend(cells_length_send, cells_num+1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 1,
                  sysp->MD_comm, &request);
        MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 1,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        sysp->subDomain.numberOfParticles += cells_length_receive[cells_num] - sum_length_send;
        sum_length_send *= sizeof(TParticle);
        particles_send = (TParticle *)malloc(sum_length_send);
        sum_length_receive = cells_length_receive[cells_num] * sizeof(TParticle);
        particles_receive = (TParticle *)malloc(sum_length_receive);
        k = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
        {
            item_pp = &sysp->subDomain.grid[ic[0]][ic[1]][ic[2]];
            item_p = *item_pp;
            while (item_p != NULL)
            {
                removeFromList(item_pp);
                particles_send[k++] = item_p->P;
                free(item_p);
                item_p = *item_pp;
            }
        }
        MPI_Isend(particles_send, sum_length_send, MPI_CHAR,
                  sysp->subDomain.rank_of_lower_subd[d], 2, sysp->MD_comm, &request);
        MPI_Recv(particles_receive, sum_length_receive, MPI_CHAR,
                 sysp->subDomain.rank_of_upper_subd[d], 2, sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(particles_send);
        kreceive = k = 0;
        ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
        {
            for (cc=0; cc<cells_length_receive[kreceive]; cc++)
            {
                item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
                item_p->P = particles_receive[k++];
                if (sysp->subDomain.is[d] == sysp->ns[d] - 1)
                    item_p->P.x[d] += sysp->l[d];
                insertInList(&sysp->subDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
            }
            kreceive++;
        }
        free(particles_receive);
        // sending to upper process, receiving from lower process
        k = sum_length_send = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
        {
            cells_length_send[k] = getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
            sum_length_send += cells_length_send[k++];
        }
        cells_length_send[cells_num] = sum_length_send;
        MPI_Isend(cells_length_send, cells_num+1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 3,
                  sysp->MD_comm, &request);
        MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 3,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(cells_length_send);
        sysp->subDomain.numberOfParticles += cells_length_receive[cells_num] - sum_length_send;
        sum_length_send *= sizeof(TParticle);
        particles_send = (TParticle *)malloc(sum_length_send);
        sum_length_receive = cells_length_receive[cells_num] * sizeof(TParticle);
        particles_receive = (TParticle *)malloc(sum_length_receive);
        k = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
        {
            item_pp = &sysp->subDomain.grid[ic[0]][ic[1]][ic[2]];
            item_p = *item_pp;
            while (item_p != NULL)
            {
                removeFromList(item_pp);
                particles_send[k++] = item_p->P;
                free(item_p);
                item_p = *item_pp;
            }
        }
        MPI_Isend(particles_send, sum_length_send, MPI_CHAR,
                  sysp->subDomain.rank_of_upper_subd[d], 4, sysp->MD_comm, &request);
        MPI_Recv(particles_receive, sum_length_receive, MPI_CHAR,
                 sysp->subDomain.rank_of_lower_subd[d], 4, sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(particles_send);
        kreceive = k = 0;
        ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
        {
            for (cc=0; cc<cells_length_receive[kreceive]; cc++)
            {
                item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
                item_p->P = particles_receive[k++];
                if (sysp->subDomain.is[d] == 0)
                    item_p->P.x[d] -= sysp->l[d];
                insertInList(&sysp->subDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
            }
            kreceive++;
        }
        free(cells_length_receive);
        free(particles_receive);
    }
}

static void ghostparticles_init_in_direction_d(
    fmd_sys_t *sysp, int d, int ic_start_send_lower[3], int ic_stop_send_lower[3],
    int ic_start_receive_lower[3], int ic_stop_receive_lower[3], int ic_start_send_upper[3],
    int ic_stop_send_upper[3], int ic_start_receive_upper[3], int ic_stop_receive_upper[3])
{
    MPI_Status status;
    MPI_Request request;
    int sum_length_send, sum_length_receive;
    int k, kreceive, cells_num;
    int *cells_length_send, *cells_length_receive, ic[3];
    TPosition_Struct *data_send, *data_receive;
    int dd;
    TParticleListItem *item_p;
    int cc;

    if ( ((sysp->subDomain.is[d] == 0) || (sysp->subDomain.is[d] == sysp->ns[d]-1)) && !sysp->PBC[d] )
    {
        if (sysp->subDomain.is[d] == 0)
        {
            // receiving from upper process
            cells_num = 1;
            for (dd=0; dd<3; dd++)
                cells_num *= ic_stop_receive_upper[dd] - ic_start_receive_upper[dd];
            cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
            MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 5,
                     sysp->MD_comm, &status);
            sum_length_receive = sizeof(TPosition_Struct) * cells_length_receive[cells_num];
            data_receive = (TPosition_Struct *)malloc(sum_length_receive);
            MPI_Recv(data_receive, sum_length_receive, MPI_CHAR,
                     sysp->subDomain.rank_of_upper_subd[d], 6, sysp->MD_comm, &status);
            kreceive = 0;
            k = -1;
            ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
            {
                k += cells_length_receive[kreceive];
                for (cc=0; cc<cells_length_receive[kreceive]; cc++)
                {
                    item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
                    for (dd=0; dd<3; dd++)
                        item_p->P.x[dd] = data_receive[k].x[dd];
                    item_p->P.elementID = data_receive[k].elementID;
                    item_p->P.groupID = data_receive[k--].groupID;
                    insertInList(&sysp->subDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
                }
                k += cells_length_receive[kreceive];
                kreceive++;
            }
            free(cells_length_receive);
            free(data_receive);
            // sending to upper process
            cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
            k = sum_length_send = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            {
                cells_length_send[k] = getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
                sum_length_send += cells_length_send[k++];
            }
            cells_length_send[cells_num] = sum_length_send;
            MPI_Send(cells_length_send, cells_num+1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 7,
                     sysp->MD_comm);
            free(cells_length_send);
            sum_length_send *= sizeof(TPosition_Struct);
            data_send = (TPosition_Struct *)malloc(sum_length_send);
            k = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
                for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                {
                    for (dd=0; dd<3; dd++)
                        data_send[k].x[dd] = item_p->P.x[dd];
                    data_send[k].elementID = item_p->P.elementID;
                    data_send[k++].groupID = item_p->P.groupID;
                }
            MPI_Send(data_send, sum_length_send, MPI_CHAR, sysp->subDomain.rank_of_upper_subd[d],
                     8, sysp->MD_comm);
            free(data_send);
        }
        else // if (sysp->subDomain.is[d] == sysp->ns[d]-1)
        {
            // sending to lower process
            cells_num = 1;
            for (dd=0; dd<3; dd++)
                cells_num *= ic_stop_send_lower[dd] - ic_start_send_lower[dd];
            cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
            k = sum_length_send = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            {
                cells_length_send[k] = getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
                sum_length_send += cells_length_send[k++];
            }
            cells_length_send[cells_num] = sum_length_send;
            MPI_Send(cells_length_send, cells_num+1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 5,
                     sysp->MD_comm);
            free(cells_length_send);
            sum_length_send *= sizeof(TPosition_Struct);
            data_send = (TPosition_Struct *)malloc(sum_length_send);
            k = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
                for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                {
                    for (dd=0; dd<3; dd++)
                        data_send[k].x[dd] = item_p->P.x[dd];
                    data_send[k].elementID = item_p->P.elementID;
                    data_send[k++].groupID = item_p->P.groupID;
                }
            MPI_Send(data_send, sum_length_send, MPI_CHAR, sysp->subDomain.rank_of_lower_subd[d],
                     6, sysp->MD_comm);
            free(data_send);
            // receiving from lower process
            cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
            MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 7,
                     sysp->MD_comm, &status);
            sum_length_receive = sizeof(TPosition_Struct) * cells_length_receive[cells_num];
            data_receive = (TPosition_Struct *)malloc(sum_length_receive);
            MPI_Recv(data_receive, sum_length_receive, MPI_CHAR,
                     sysp->subDomain.rank_of_lower_subd[d], 8, sysp->MD_comm, &status);
            kreceive = 0;
            k = -1;
            ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
            {
                k += cells_length_receive[kreceive];
                for (cc=0; cc<cells_length_receive[kreceive]; cc++)
                {
                    item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
                    for (dd=0; dd<3; dd++)
                        item_p->P.x[dd] = data_receive[k].x[dd];
                    item_p->P.elementID = data_receive[k].elementID;
                    item_p->P.groupID = data_receive[k--].groupID;
                    insertInList(&sysp->subDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
                }
                k += cells_length_receive[kreceive];
                kreceive++;
            }
            free(cells_length_receive);
            free(data_receive);
        }
    }
    else
    {
        // sending to lower process, receiving from upper process
        cells_num = 1;
        for (dd=0; dd<3; dd++)
            cells_num *= ic_stop_send_lower[dd] - ic_start_send_lower[dd];
        cells_length_send = (int *)malloc((cells_num+1) * sizeof(int));
        cells_length_receive = (int *)malloc((cells_num+1) * sizeof(int));
        k = sum_length_send = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
        {
            cells_length_send[k] = getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
            sum_length_send += cells_length_send[k++];
        }
        cells_length_send[cells_num] = sum_length_send;
        MPI_Isend(cells_length_send, cells_num+1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 5,
                  sysp->MD_comm, &request);
        MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 5,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        sum_length_send *= sizeof(TPosition_Struct);
        data_send = (TPosition_Struct *)malloc(sum_length_send);
        sum_length_receive = sizeof(TPosition_Struct) * cells_length_receive[cells_num];
        data_receive = (TPosition_Struct *)malloc(sum_length_receive);
        k = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            {
                for (dd=0; dd<3; dd++)
                    data_send[k].x[dd] = item_p->P.x[dd];
                data_send[k].elementID = item_p->P.elementID;
                data_send[k++].groupID = item_p->P.groupID;
            }
        MPI_Isend(data_send, sum_length_send, MPI_CHAR, sysp->subDomain.rank_of_lower_subd[d], 6,
                  sysp->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, MPI_CHAR,
                 sysp->subDomain.rank_of_upper_subd[d], 6, sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);
        kreceive = 0;
        k = -1;
        ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
        {
            k += cells_length_receive[kreceive];
            for (cc=0; cc<cells_length_receive[kreceive]; cc++)
            {
                item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
                for (dd=0; dd<3; dd++)
                    item_p->P.x[dd] = data_receive[k].x[dd];
                item_p->P.elementID = data_receive[k].elementID;
                item_p->P.groupID = data_receive[k--].groupID;
                if (sysp->subDomain.is[d] == sysp->ns[d] - 1)
                    item_p->P.x[d] += sysp->l[d];
                insertInList(&sysp->subDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
            }
            k += cells_length_receive[kreceive];
            kreceive++;
        }
        free(data_receive);
        // sending to upper process, receiving from lower process
        k = sum_length_send = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
        {
            cells_length_send[k] = getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
            sum_length_send += cells_length_send[k++];
        }
        cells_length_send[cells_num] = sum_length_send;
        MPI_Isend(cells_length_send, cells_num+1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 7,
                  sysp->MD_comm, &request);
        MPI_Recv(cells_length_receive, cells_num+1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 7,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(cells_length_send);
        sum_length_send *= sizeof(TPosition_Struct);
        data_send = (TPosition_Struct *)malloc(sum_length_send);
        sum_length_receive = sizeof(TPosition_Struct) * cells_length_receive[cells_num];
        data_receive = (TPosition_Struct *)malloc(sum_length_receive);
        k = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
            {
                for (dd=0; dd<3; dd++)
                    data_send[k].x[dd] = item_p->P.x[dd];
                data_send[k].elementID = item_p->P.elementID;
                data_send[k++].groupID = item_p->P.groupID;
            }
        MPI_Isend(data_send, sum_length_send, MPI_CHAR, sysp->subDomain.rank_of_upper_subd[d], 8,
                  sysp->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, MPI_CHAR,
                 sysp->subDomain.rank_of_lower_subd[d], 8, sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);
        kreceive = 0;
        k = -1;
        ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
        {
            k += cells_length_receive[kreceive];
            for (cc=0; cc<cells_length_receive[kreceive]; cc++)
            {
                item_p = (TParticleListItem *)malloc(sizeof(TParticleListItem));
                for (dd=0; dd<3; dd++)
                    item_p->P.x[dd] = data_receive[k].x[dd];
                item_p->P.elementID = data_receive[k].elementID;
                item_p->P.groupID = data_receive[k--].groupID;
                if (sysp->subDomain.is[d] == 0)
                    item_p->P.x[d] -= sysp->l[d];
                insertInList(&sysp->subDomain.grid[ic[0]][ic[1]][ic[2]], item_p);
            }
            k += cells_length_receive[kreceive];
            kreceive++;
        }
        free(cells_length_receive);
        free(data_receive);
    }
}

static void ghostparticles_update_Fprime_in_direction_d(
    fmd_sys_t *sysp, int d, int ic_start_send_lower[3],
    int ic_stop_send_lower[3], int ic_start_receive_lower[3],
    int ic_stop_receive_lower[3], int ic_start_send_upper[3], int ic_stop_send_upper[3],
    int ic_start_receive_upper[3], int ic_stop_receive_upper[3])
{
    MPI_Status status;
    MPI_Request request;
    int sum_length_send, sum_length_receive;
    int ic[3];
    double *data_receive, *data_send;
    TParticleListItem *item_p;
    int k;

    if ( ((sysp->subDomain.is[d] == 0) || (sysp->subDomain.is[d] == sysp->ns[d]-1)) && !sysp->PBC[d] )
    {
        if (sysp->subDomain.is[d] == 0)
        {
            // receiving from upper process
            MPI_Recv(&sum_length_receive, 1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 100,
                     sysp->MD_comm, &status);
            data_receive = (double *)malloc(sum_length_receive * sizeof(double));
            MPI_Recv(data_receive, sum_length_receive, MPI_DOUBLE,
                     sysp->subDomain.rank_of_upper_subd[d], 101, sysp->MD_comm, &status);
            k = 0;
            ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
                for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                    item_p->FembPrime = data_receive[k++];
            free(data_receive);
            // sending to upper process
            sum_length_send = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
                sum_length_send += getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
            MPI_Send(&sum_length_send, 1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 102,
                     sysp->MD_comm);
            data_send = (double *)malloc(sum_length_send * sizeof(double));
            k = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
                for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                    data_send[k++] = item_p->FembPrime;
            MPI_Send(data_send, sum_length_send, MPI_DOUBLE, sysp->subDomain.rank_of_upper_subd[d],
                     103, sysp->MD_comm);
            free(data_send);
        }
        else // if (sysp->subDomain.is[d] == sysp->ns[d]-1)
        {
            // sending to lower process
            sum_length_send = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
                sum_length_send += getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
            MPI_Send(&sum_length_send, 1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 100,
                     sysp->MD_comm);
            data_send = (double *)malloc(sum_length_send * sizeof(double));
            k = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
                for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                    data_send[k++] = item_p->FembPrime;
            MPI_Send(data_send, sum_length_send, MPI_DOUBLE, sysp->subDomain.rank_of_lower_subd[d],
                     101, sysp->MD_comm);
            free(data_send);
            // receiving from lower process
            MPI_Recv(&sum_length_receive, 1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 102,
                     sysp->MD_comm, &status);
            data_receive = (double *)malloc(sum_length_receive * sizeof(double));
            MPI_Recv(data_receive, sum_length_receive, MPI_DOUBLE,
                     sysp->subDomain.rank_of_lower_subd[d], 103, sysp->MD_comm, &status);
            k = 0;
            ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
                for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                    item_p->FembPrime = data_receive[k++];
            free(data_receive);
        }
    }
    else
    {
        // sending to lower process, receiving from upper process
        sum_length_send = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            sum_length_send += getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
        MPI_Isend(&sum_length_send, 1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 100,
                  sysp->MD_comm, &request);
        MPI_Recv(&sum_length_receive, 1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 100,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        data_send = (double *)malloc(sum_length_send * sizeof(double));
        data_receive = (double *)malloc(sum_length_receive * sizeof(double));
        k = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                data_send[k++] = item_p->FembPrime;
        MPI_Isend(data_send, sum_length_send, MPI_DOUBLE, sysp->subDomain.rank_of_lower_subd[d], 101,
                  sysp->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, MPI_DOUBLE, sysp->subDomain.rank_of_upper_subd[d], 101,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);
        k = 0;
        ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                item_p->FembPrime = data_receive[k++];
        free(data_receive);
        // sending to upper process, receiving from lower process
        sum_length_send = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            sum_length_send += getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
        MPI_Isend(&sum_length_send, 1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 102,
                  sysp->MD_comm, &request);
        MPI_Recv(&sum_length_receive, 1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 102,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        data_send = (double *)malloc(sum_length_send * sizeof(double));
        data_receive = (double *)malloc(sum_length_receive * sizeof(double));
        k = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                data_send[k++] = item_p->FembPrime;
        MPI_Isend(data_send, sum_length_send, MPI_DOUBLE, sysp->subDomain.rank_of_upper_subd[d], 103,
                  sysp->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, MPI_DOUBLE, sysp->subDomain.rank_of_lower_subd[d], 103,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);
        k = 0;
        ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                item_p->FembPrime = data_receive[k++];
        free(data_receive);
    }
}

static void ghostparticles_update_LocOrdParam_in_direction_d(
    fmd_sys_t *sysp, int d, int ic_start_send_lower[3],
    int ic_stop_send_lower[3], int ic_start_receive_lower[3],
    int ic_stop_receive_lower[3], int ic_start_send_upper[3], int ic_stop_send_upper[3],
    int ic_start_receive_upper[3], int ic_stop_receive_upper[3])
{
    MPI_Status status;
    MPI_Request request;
    int sum_length_send, sum_length_receive;
    int ic[3];
    float *data_receive, *data_send;
    TParticleListItem *item_p;
    int k;

    if ( ((sysp->subDomain.is[d] == 0) || (sysp->subDomain.is[d] == sysp->ns[d]-1)) && !sysp->PBC[d] )
    {
        if (sysp->subDomain.is[d] == 0)
        {
            // receiving from upper process
            MPI_Recv(&sum_length_receive, 1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 130,
                     sysp->MD_comm, &status);
            data_receive = (float *)malloc(sum_length_receive * sizeof(float));
            MPI_Recv(data_receive, sum_length_receive, MPI_FLOAT,
                     sysp->subDomain.rank_of_upper_subd[d], 131, sysp->MD_comm, &status);
            k = 0;
            ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
                for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                    item_p->P.LocOrdParam = data_receive[k++];
            free(data_receive);
            // sending to upper process
            sum_length_send = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
                sum_length_send += getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
            MPI_Send(&sum_length_send, 1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 132,
                     sysp->MD_comm);
            data_send = (float *)malloc(sum_length_send * sizeof(float));
            k = 0;
            ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
                for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                    data_send[k++] = item_p->P.LocOrdParam;
            MPI_Send(data_send, sum_length_send, MPI_FLOAT, sysp->subDomain.rank_of_upper_subd[d],
                     133, sysp->MD_comm);
            free(data_send);
        }
        else // if (sysp->subDomain.is[d] == sysp->ns[d]-1)
        {
            // sending to lower process
            sum_length_send = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
                sum_length_send += getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
            MPI_Send(&sum_length_send, 1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 130,
                     sysp->MD_comm);
            data_send = (float *)malloc(sum_length_send * sizeof(float));
            k = 0;
            ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
                for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                    data_send[k++] = item_p->P.LocOrdParam;
            MPI_Send(data_send, sum_length_send, MPI_FLOAT, sysp->subDomain.rank_of_lower_subd[d],
                     131, sysp->MD_comm);
            free(data_send);
            // receiving from lower process
            MPI_Recv(&sum_length_receive, 1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 132,
                     sysp->MD_comm, &status);
            data_receive = (float *)malloc(sum_length_receive * sizeof(float));
            MPI_Recv(data_receive, sum_length_receive, MPI_FLOAT,
                     sysp->subDomain.rank_of_lower_subd[d], 133, sysp->MD_comm, &status);
            k = 0;
            ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
                for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                    item_p->P.LocOrdParam = data_receive[k++];
            free(data_receive);
        }
    }
    else
    {
        // sending to lower process, receiving from upper process
        sum_length_send = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            sum_length_send += getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
        MPI_Isend(&sum_length_send, 1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 130,
                  sysp->MD_comm, &request);
        MPI_Recv(&sum_length_receive, 1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 130,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        data_send = (float *)malloc(sum_length_send * sizeof(float));
        data_receive = (float *)malloc(sum_length_receive * sizeof(float));
        k = 0;
        ITERATE(ic, ic_start_send_lower, ic_stop_send_lower)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                data_send[k++] = item_p->P.LocOrdParam;
        MPI_Isend(data_send, sum_length_send, MPI_FLOAT, sysp->subDomain.rank_of_lower_subd[d], 131,
                  sysp->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, MPI_FLOAT, sysp->subDomain.rank_of_upper_subd[d], 131,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);
        k = 0;
        ITERATE(ic, ic_start_receive_upper, ic_stop_receive_upper)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                item_p->P.LocOrdParam = data_receive[k++];
        free(data_receive);
        // sending to upper process, receiving from lower process
        sum_length_send = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            sum_length_send += getListLength(sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]);
        MPI_Isend(&sum_length_send, 1, MPI_INT, sysp->subDomain.rank_of_upper_subd[d], 132,
                  sysp->MD_comm, &request);
        MPI_Recv(&sum_length_receive, 1, MPI_INT, sysp->subDomain.rank_of_lower_subd[d], 132,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        data_send = (float *)malloc(sum_length_send * sizeof(float));
        data_receive = (float *)malloc(sum_length_receive * sizeof(float));
        k = 0;
        ITERATE(ic, ic_start_send_upper, ic_stop_send_upper)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                data_send[k++] = item_p->P.LocOrdParam;
        MPI_Isend(data_send, sum_length_send, MPI_FLOAT, sysp->subDomain.rank_of_upper_subd[d], 133,
                  sysp->MD_comm, &request);
        MPI_Recv(data_receive, sum_length_receive, MPI_FLOAT, sysp->subDomain.rank_of_lower_subd[d], 133,
                 sysp->MD_comm, &status);
        MPI_Wait(&request, &status);
        free(data_send);
        k = 0;
        ITERATE(ic, ic_start_receive_lower, ic_stop_receive_lower)
            for (item_p = sysp->subDomain.grid[ic[0]][ic[1]][ic[2]]; item_p != NULL; item_p = item_p->next_p)
                item_p->P.LocOrdParam = data_receive[k++];
        free(data_receive);
    }
}

static void particles_prepare_migration_in_direction_d(
    TSubDomain *s_p, int d, int ic_start_send_lower[3],
    int ic_stop_send_lower[3], int ic_start_receive_lower[3],
    int ic_stop_receive_lower[3], int ic_start_send_upper[3], int ic_stop_send_upper[3],
    int ic_start_receive_upper[3], int ic_stop_receive_upper[3])
{
    int dd;

    for (dd=0; dd<3; dd++)
    {
        if (dd == d) // only ghost
        {
            ic_start_send_lower[dd] = s_p->ic_start[dd];
            ic_stop_send_lower[dd] = ic_start_send_lower[dd] + s_p->ic_start[dd];
            ic_start_receive_lower[dd] = 0;
            ic_stop_receive_lower[dd] = s_p->ic_start[dd];
            ic_stop_send_upper[dd] = s_p->ic_stop[dd];
            ic_start_send_upper[dd] = ic_stop_send_upper[dd] - s_p->ic_start[dd];
            ic_start_receive_upper[dd] = s_p->ic_stop[dd];
            ic_stop_receive_upper[dd] = s_p->cell_num[dd];
        }
        else if (dd > d) // including ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd] = 0;
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd] = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd] = s_p->cell_num[dd];
        }
        else // excluding ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd] = s_p->ic_start[dd];
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd] = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd] = s_p->ic_stop[dd];
        }
    }
}

static void ghostparticles_prepare_init_update_in_direction_d(
    fmd_sys_t *sysp, int d, int ic_start_send_lower[3],
    int ic_stop_send_lower[3], int ic_start_receive_lower[3],
    int ic_stop_receive_lower[3], int ic_start_send_upper[3], int ic_stop_send_upper[3],
    int ic_start_receive_upper[3], int ic_stop_receive_upper[3])
{
    int dd;

    for (dd=0; dd<3; dd++)
    {
        if (dd == d) // only ghost
        {
            ic_start_send_lower[dd] = sysp->subDomain.ic_start[dd];
            ic_stop_send_lower[dd] = ic_start_send_lower[dd] + 1;
            ic_stop_receive_lower[dd] = sysp->subDomain.ic_start[dd];
            ic_start_receive_lower[dd] = ic_stop_receive_lower[dd] - 1;
            ic_stop_send_upper[dd] = sysp->subDomain.ic_stop[dd];
            ic_start_send_upper[dd] = ic_stop_send_upper[dd] - 1;
            ic_start_receive_upper[dd] = sysp->subDomain.ic_stop[dd];
            ic_stop_receive_upper[dd] = ic_start_receive_upper[dd] + 1;
        }
        else if (dd > d) // including ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd]
                                       = sysp->subDomain.ic_start[dd] - (sysp->ns[dd] == 1 ? 0 : 1);
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd]
                                       = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd]
                                       = sysp->subDomain.ic_stop[dd] + (sysp->ns[dd] == 1 ? 0 : 1);
        }
        else // excluding ghost
        {
            ic_start_receive_lower[dd] = ic_start_send_lower[dd]
                                       = ic_start_receive_upper[dd]
                                       = ic_start_send_upper[dd]
                                       = sysp->subDomain.ic_start[dd];
            ic_stop_receive_lower[dd]  = ic_stop_send_lower[dd]
                                       = ic_stop_receive_upper[dd]
                                       = ic_stop_send_upper[dd]
                                       = sysp->subDomain.ic_stop[dd];
        }
    }
}

void fmd_ghostparticles_update_LocOrdParam(fmd_sys_t *sysp)
{
    int d;
    int ic_start_send_lower[3], ic_stop_send_lower[3];
    int ic_start_send_upper[3], ic_stop_send_upper[3];
    int ic_start_receive_lower[3], ic_stop_receive_lower[3];
    int ic_start_receive_upper[3], ic_stop_receive_upper[3];

    for (d = 3-1; d >= 0; d--)
    {
        if (sysp->ns[d] != 1)
        {
            ghostparticles_prepare_init_update_in_direction_d(
                sysp, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
            ghostparticles_update_LocOrdParam_in_direction_d(
                sysp, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
        }
    }
}

void fmd_ghostparticles_update_Femb(fmd_sys_t *sysp)
{
    int d;
    int ic_start_send_lower[3], ic_stop_send_lower[3];
    int ic_start_send_upper[3], ic_stop_send_upper[3];
    int ic_start_receive_lower[3], ic_stop_receive_lower[3];
    int ic_start_receive_upper[3], ic_stop_receive_upper[3];

    for (d = 3-1; d >= 0; d--)
    {
        if (sysp->ns[d] != 1)
        {
            ghostparticles_prepare_init_update_in_direction_d(
                sysp, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
            ghostparticles_update_Fprime_in_direction_d(
                sysp, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
        }
    }
}

void fmd_ghostparticles_init(fmd_sys_t *sysp)
{
    int d;
    int ic_start_send_lower[3], ic_stop_send_lower[3];
    int ic_start_send_upper[3], ic_stop_send_upper[3];
    int ic_start_receive_lower[3], ic_stop_receive_lower[3];
    int ic_start_receive_upper[3], ic_stop_receive_upper[3];

    for (d = 3-1; d >= 0; d--)
    {
        if (sysp->ns[d] != 1)
        {
            ghostparticles_prepare_init_update_in_direction_d(
                sysp, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
            ghostparticles_init_in_direction_d(
                sysp, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
        }
    }
}

void fmd_ghostparticles_delete(fmd_sys_t *sysp)
{
    int d;
    int ic_from[3], ic_to[3];
    int jc[3];

    for (d=0; d<3; d++)
    {
        ic_from[d] = sysp->subDomain.ic_start[d] - (sysp->ns[d] == 1 ? 0 : 1);
        jc[d] = ic_to[d] = sysp->subDomain.ic_stop[d] + (sysp->ns[d] == 1 ? 0 : 1);
    }

    for (d=0; d<3; d++)
    {
        jc[d] = sysp->subDomain.ic_start[d];
        cleanGridSegment(sysp->subDomain.grid, ic_from, jc);
        ic_from[d] = sysp->subDomain.ic_stop[d];
        cleanGridSegment(sysp->subDomain.grid, ic_from, ic_to);
        ic_from[d] = sysp->subDomain.ic_start[d];
        jc[d] = ic_to[d] = sysp->subDomain.ic_stop[d];
    }
}

void fmd_particles_migrate(fmd_sys_t *sysp)
{
    int ic_start_send_lower[3], ic_stop_send_lower[3];
    int ic_start_send_upper[3], ic_stop_send_upper[3];
    int ic_start_receive_lower[3], ic_stop_receive_lower[3];
    int ic_start_receive_upper[3], ic_stop_receive_upper[3];
    int d;

    for (d = 0; d < 3; d++)
    {
        if (sysp->ns[d] != 1)
        {
            particles_prepare_migration_in_direction_d(
                &sysp->subDomain, d,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_lower, ic_stop_send_lower, ic_start_receive_upper,
                ic_stop_receive_upper, ic_start_send_upper, ic_stop_send_upper);
            particles_migrate_in_direction_d(
                sysp, d, ic_start_send_lower, ic_stop_send_lower,
                ic_start_receive_lower, ic_stop_receive_lower,
                ic_start_send_upper, ic_stop_send_upper, ic_start_receive_upper,
                ic_stop_receive_upper);
        }
    }
}
