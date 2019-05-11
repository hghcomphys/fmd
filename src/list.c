/*
  list.c: This file is part of Free Molecular Dynamics

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

// this module owes some ideas to the doubly-linked list of GLib

#include "list.h"
#include <stdlib.h>

list_t *fmd_list_prepend(list_t *list, void *data)
{
    list_t *item = malloc(sizeof(list_t));
    if (item != NULL)
    {
        item->data = data;
        item->next = list;
        if (list == NULL)
            item->prev = NULL;
        else
        {
            item->prev = list->prev;
            if (list->prev != NULL)
                list->prev->next = item;
            list->prev = item;
        }
    }
    return item;
}

list_t *fmd_list_item_remove(list_t *list, list_t *item)
{
}
