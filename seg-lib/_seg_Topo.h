#pragma once

#include "_seg_common.h"

void ConnectComp(int * Old, int * New, int dimensions[3],int varin);

void Dillate(bool * Image, int kernel, int dimensions[3],int verbose);

void Dillate_const(bool * Image,
                   bool * Const,
                   int kernel,
                   int dimensions[3],
                   int direction);
