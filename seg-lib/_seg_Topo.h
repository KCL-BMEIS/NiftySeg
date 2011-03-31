#pragma once

#include "_seg_common.h"

void ConnectComp(int * Old, int * New, int dimensions[3]);

void Dillate(bool * Image, int kernel, int dimensions[3]);

void Dillate_const(bool * Image,
                   bool * Const,
                   int kernel,
                   int dimensions[3]);
