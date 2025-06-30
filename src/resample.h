#ifndef PRONAME_RESAMPLE_H
#define PRONAME_RESAMPLE_H
#include <arrayfire.h>
#include "tool.h"
using namespace af;

array resample(const array& x, int p, int q);

#endif // PRONAME_RESAMPLE_H
