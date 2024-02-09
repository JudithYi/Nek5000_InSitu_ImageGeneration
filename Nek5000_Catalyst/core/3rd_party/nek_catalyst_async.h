
#ifndef nek_catalyst_h
#define nek_catalyst_h

#include <vtkIdTypeArray.h>
#include <vtkUnsignedCharArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkUnstructuredGrid.h>
#include <vtkCellArray.h>
#include "vtkDoubleArray.h"
#include "vtkCPDataDescription.h"
#include "vtkCPInputDataDescription.h"
#include "vtkCPProcessor.h"
#include "CPythonAdaptorAPI.h"
#include "CAdaptorAPI.h"
#include "myCPPythonAdaptorAPI.h"


void catalyst_usrpipe(char* name, int length=9);

void creategrid(const double *x, const double *y, const double *z,
			    const int *lx1, const int *ly1, const int *lz1, 
			    const int *lelt, const int *dim);

void add_scalar_field(double *data, char *name) ;

void add_vector_field(double *xdata, double *ydata, double *zdata, 
				  int *dim, char *name);

#endif
// HeaderTest-Exclude: nek_catalyst.h

