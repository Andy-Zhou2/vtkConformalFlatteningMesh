#ifndef vtkConformalFlatteningMeshFilter_h
#define vtkConformalFlatteningMeshFilter_h

#include "vtkPolyData.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkVector.h"
#include "vtkSmartPointer.h"
#include "vtkTriangle.h"
#include "vtkCellArray.h"

#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkInformationVector.h"
#include "vtkDataObject.h"
#include "vtkInformation.h"

// vnl headers
#include "vnl/vnl_sparse_matrix.h"
#include "vnl/vnl_vector.h"


namespace vtk
{

  class ConformalFlatteningMeshFilter : public vtkPolyDataAlgorithm
  {
  public:
    /** Standard class typedefs. */
    typedef ConformalFlatteningMeshFilter Self;

    //typedef vtkPolyDataAlgorithm Superclass;

    /** Select the cell that will be used as reference for the flattening.
     * This value must be the identifier of a cell existing in the input Mesh.
     * A point of this cell will be mapped to infinity on the plane, or it
     * will be mapped to the north-pole on the sphere. It is recommended to
     * select a cell whose curvature is relatively flat. */
    void SetPolarCellIdentifier(vtkIdType cellId);

    /** Define the scale of the mapping. The largest coordinates of the
     * furthest point in the plane is m_MapScale. */
    void SetScale(double);

    /** Define that the input surface will be mapped to a sphere */
    void MapToSphere();

    /** Define that the input surface will be mapped to a plane.
     *  This skips the steps of the stereographic projection. */
    void MapToPlane();
    //vtkTypeMacro(ConformalFlatteningMeshFilter, vtkPolyDataAlgorithm);
    void PrintSelf(std::ostream & os, vtkIndent indent);
    vtkTypeMacro(ConformalFlatteningMeshFilter, vtkPolyDataAlgorithm);
    static ConformalFlatteningMeshFilter *New();

  protected:
    ConformalFlatteningMeshFilter();
    ~ConformalFlatteningMeshFilter() {}

    /** Generate Requested Data */
    int RequestData(vtkInformation *, vtkInformationVector **, vtkInformationVector *);


  private:
    ConformalFlatteningMeshFilter(const ConformalFlatteningMeshFilter&);

    typedef vnl_vector< double >        VectorCoordType;
    typedef vnl_sparse_matrix< double > SparseMatrixCoordType;

    /** Cell Id  in which the point P, which is used
     * to define the mapping, lies in. */
    unsigned int m_PolarCellIdentifier;

    /** Whether the result is sphere or plane.  */
    bool m_MapToSphere;

    /** The scale when mapping to the plane.
     *  Determines how far the farthest point goes. */
    double m_MapScale;
		
    void operator=(const ConformalFlatteningMeshFilter&);  // Not implemented.
  };
} // end namespace vtk

#endif
