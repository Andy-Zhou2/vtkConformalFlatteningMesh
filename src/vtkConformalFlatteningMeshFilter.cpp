#include "vtkConformalFlatteningMeshFilter.h"

// #ifndef vtkConformalFlatteningMeshFilter_hxx
// #define vtkConformalFlatteningMeshFilter_hxx
#include "vtkObjectFactory.h"

#include <cfloat>  // for DBL_MIN

//#include <vtkXMLPolyDataWriter.h>
namespace vtk
{
	vtkStandardNewMacro(ConformalFlatteningMeshFilter);

	ConformalFlatteningMeshFilter ::ConformalFlatteningMeshFilter()
	{
		this->SetNumberOfInputPorts(1); 
		this->SetNumberOfOutputPorts(1);
		
		this->m_PolarCellIdentifier = 0;
		this->m_MapToSphere = false;

		this->m_MapScale = -1.0;
		// If during the stage when this parameter is used and it is still
		// -1.0, then it indicates that the user doesn't assign a scale
		// factor. Then automatically calculate it s.t. after doing the
		// stereo-graphic projection, upper and lower hemi-sphere will have
		// same number of vertics.
	}

	/**
	 * Set the triangle used to define the boundary of the flattened region.
	 */

	void
		ConformalFlatteningMeshFilter ::SetPolarCellIdentifier(vtkIdType cellId)
	{
		this->m_PolarCellIdentifier = cellId;
	}

	/**
	 * Define the scale of the mapping. The largest coordinates of the
	 * furthest point in the plane is m_MapScale.
	 */

	void ConformalFlatteningMeshFilter ::SetScale(double scale)
	{
		this->m_MapScale = scale;
	}

	/**
	 * Define that the input surface will be mapped to a sphere
	 */
	void ConformalFlatteningMeshFilter ::MapToSphere(void)
	{
		this->m_MapToSphere = true;
	}

	/** Define that the input surface will be mapped to a plane.
	 *  This skips the steps of the stereographic projection.
	 */
	
	void ConformalFlatteningMeshFilter ::MapToPlane(void)
	{
		this->m_MapToSphere = false;
	}



	void ConformalFlatteningMeshFilter ::PrintSelf(std::ostream & os, vtkIndent indent)
	{
		Superclass::PrintSelf(os, indent);
	}

	/**
	 * This method causes the filter to generate its output.
	 */

	int ConformalFlatteningMeshFilter::RequestData(vtkInformation *vtkNotUsed(request),
			vtkInformationVector **inputVector,
			vtkInformationVector *outputVector)
	{
		vtkDebugMacro(<< 3);
		vtkDebugMacro(<< "Into function RequestData" << endl);
		// get the input and output
		vtkSmartPointer<vtkPolyData> inputMesh = vtkPolyData::GetData(inputVector[0], 0);
		vtkSmartPointer<vtkPolyData> output = vtkPolyData::GetData(outputVector, 0);

		vtkSmartPointer<vtkPolyData> outputMesh = vtkSmartPointer<vtkPolyData>::New(); // Create an empty mesh and copy it to output at last.

		if (!inputMesh)
		{
			vtkDebugMacro(<< "Missing Input Mesh");
		}

		if (!outputMesh)
		{
			vtkDebugMacro(<< "Missing Output Mesh");
		}

		vtkSmartPointer<vtkPoints> outPoints = vtkSmartPointer<vtkPoints>::New();

		const unsigned int numberOfPoints = inputMesh->GetNumberOfPoints();


		unsigned int i;

		SparseMatrixCoordType D(numberOfPoints, numberOfPoints);

		VectorCoordType bx(numberOfPoints, 0.0);
		VectorCoordType by(numberOfPoints, 0.0);

		vtkDebugMacro(<< "m_PolarCellIdentifier " << this->m_PolarCellIdentifier << endl);

		vtkSmartPointer<vtkCellArray> cells = inputMesh->GetPolys();

		vtkSmartPointer<vtkIdList> idList = vtkSmartPointer<vtkIdList>::New();
		cells->GetCell(m_PolarCellIdentifier, idList);

		unsigned int cellNumberOfPoints = idList->GetNumberOfIds();

		if (cellNumberOfPoints != 3)
		{
			vtkDebugMacro(<< "Polar cell has " << cellNumberOfPoints << " points"
				"\nThis filter can only process triangle meshes. "
				"Use vtkTriangleFilter to convert your mesh to a triangle mesh." << endl);
			return 0;
		}

		vtkIdType boundaryId0 = idList->GetId(0);
		vtkIdType boundaryId1 = idList->GetId(1);
		vtkIdType boundaryId2 = idList->GetId(2);

		double ptA[3] = { 0.0, 0.0, 0.0 };
		double ptB[3] = { 0.0, 0.0, 0.0 };
		double ptC[3] = { 0.0, 0.0, 0.0 };

		inputMesh->GetPoint(boundaryId0, ptA); 
		inputMesh->GetPoint(boundaryId1, ptB);
		inputMesh->GetPoint(boundaryId2, ptC);

		vtkDebugMacro(<< "Stage 1 finish" << endl);

		double AB[3];
		double BC[3];
		double CA[3];

		double normAB2;
		double normBC2;
		double normCA2;

		double normBC;
		double normCA;

		double prodABBC;
		double prodBCCA;
		double prodCAAB;

		AB[0] = ptB[0] - ptA[0];
		AB[1] = ptB[1] - ptA[1];
		AB[2] = ptB[2] - ptA[2];

		BC[0] = ptC[0] - ptB[0];
		BC[1] = ptC[1] - ptB[1];
		BC[2] = ptC[2] - ptB[2];

		CA[0] = ptA[0] - ptC[0];
		CA[1] = ptA[1] - ptC[1];
		CA[2] = ptA[2] - ptC[2];

		normAB2 = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];

		if (normAB2 < 1e-10)
		{
			vtkDebugMacro(<< "||AB||^2 = " << normAB2
				<< "\nRisk of division by zero" << endl);
			return 0;
		}

		double E[3];
		double CE[3];

		prodCAAB = CA[0] * AB[0] + CA[1] * AB[1] + CA[2] * AB[2];

		// E = projection of C onto AB orthogonal to AB.
		// t = Parameter to find E = A + t * ( B - A ).
		//
		// If t = 0.0,  E = A.
		// If t = 1.0,  E = B.
		//
		// |AC| * cos(alpha) = t * |AB|
		// AB * AC = |AB| |AC| cos(alpha)
		//
		// t = (|AC| / |AB|)  *  ((AB * AC) / (|AB| * |AC|))
		//   = (AB * AC) / (|AB| * |AB|)

		double t = -prodCAAB / normAB2;

		E[0] = ptA[0] + t * AB[0];
		E[1] = ptA[1] + t * AB[1];
		E[2] = ptA[2] + t * AB[2];

		CE[0] = ptC[0] - E[0];
		CE[1] = ptC[1] - E[1];
		CE[2] = ptC[2] - E[2];

		double normCE2 = CE[0] * CE[0] + CE[1] * CE[1] + CE[2] * CE[2];

		double normAB = std::sqrt(normAB2);
		double normCE = std::sqrt(normCE2);

		vtkDebugMacro(<< "scale " << this->m_MapScale << endl);

		double tmp = 2.0 / normAB;
		//double factor = normAB / normCE;

		bx(boundaryId0) = -tmp;  // -t * factor;
		bx(boundaryId1) = tmp;   // (1.0 - t) * factor;

		double tmp2 = 2.0 / normCE;

		by(boundaryId0) = tmp2 * (1.0 - t); // 0.0;
		by(boundaryId1) = tmp2 * t;           // 0.0;
		by(boundaryId2) = -tmp2;              // 1.0;

		int cellIterCount = 0;
		int numberOfCells = inputMesh->GetNumberOfCells();

		int ptIdA;
		int ptIdB;
		int ptIdC;

		double cosABC;
		double cosBCA;
		double cosCAB;

		double sinABC;
		double sinBCA;
		double sinCAB;

		double cotgABC;
		double cotgBCA;
		double cotgCAB;

		vtkDebugMacro(<< "Stage 2 finish" << endl);

		while(cellIterCount < numberOfCells)
		{
			vtkSmartPointer<vtkCell> aCell;
			aCell = inputMesh->GetCell(cellIterCount);
			unsigned int aCellNumberOfPoints = aCell->GetNumberOfPoints();



			if (aCellNumberOfPoints > 3)
			{
				vtkDebugMacro(<< "cell has " << aCellNumberOfPoints << " points\n" <<
					"This filter can only process triangle meshes." << endl);
				return 0;
			}

			while (aCellNumberOfPoints < 3) // leave the edges and points untouched
			{
				cellIterCount++;
				if (cellIterCount < numberOfCells)
				{
					aCell = inputMesh->GetCell(cellIterCount);
					aCellNumberOfPoints = aCell->GetNumberOfPoints();
				}
			}
			if (cellIterCount == numberOfCells) { break; } 

			ptIdA = aCell->GetPointId(0);
			ptIdB = aCell->GetPointId(1);
			ptIdC = aCell->GetPointId(2);

			inputMesh->GetPoint(ptIdA, ptA);
			inputMesh->GetPoint(ptIdB, ptB);
			inputMesh->GetPoint(ptIdC, ptC);

			AB[0] = ptB[0] - ptA[0];
			AB[1] = ptB[1] - ptA[1];
			AB[2] = ptB[2] - ptA[2];

			BC[0] = ptC[0] - ptB[0];
			BC[1] = ptC[1] - ptB[1];
			BC[2] = ptC[2] - ptB[2];

			CA[0] = ptA[0] - ptC[0];
			CA[1] = ptA[1] - ptC[1];
			CA[2] = ptA[2] - ptC[2];

			normAB2 = AB[0] * AB[0] + AB[1] * AB[1] + AB[2] * AB[2];
			normBC2 = BC[0] * BC[0] + BC[1] * BC[1] + BC[2] * BC[2];
			normCA2 = CA[0] * CA[0] + CA[1] * CA[1] + CA[2] * CA[2];

			if (normAB2 < 1e-10)
			{
				vtkDebugMacro(<< "normAB2 " << normAB2 << endl);
				return 0;
			}

			if (normBC2 < 1e-10)
			{
				vtkDebugMacro(<< "normBC2 " << normBC2 << endl);
				return 0;
			}

			if (normCA2 < 1e-10)
			{
				vtkDebugMacro(<< "normCA2 " << normCA2 << endl);
				return 0;
			}

			normAB = std::sqrt(normAB2);
			normBC = std::sqrt(normBC2);
			normCA = std::sqrt(normCA2);

			prodABBC = AB[0] * BC[0] + AB[1] * BC[1] + AB[2] * BC[2];
			prodBCCA = BC[0] * CA[0] + BC[1] * CA[1] + BC[2] * CA[2];
			prodCAAB = CA[0] * AB[0] + CA[1] * AB[1] + CA[2] * AB[2];

			cosABC = -prodABBC / (normAB * normBC);
			cosBCA = -prodBCCA / (normBC * normCA);
			cosCAB = -prodCAAB / (normCA * normAB);

			if (cosABC <= -1.0 || cosABC >= 1.0)
			{
				vtkDebugMacro(<< "cosABC= " << cosABC << endl);
				return 0;
			}

			if (cosBCA <= -1.0 || cosBCA >= 1.0)
			{
				vtkDebugMacro(<< "cosBCA= " << cosBCA << endl);
				return 0;
			}

			if (cosCAB <= -1.0 || cosCAB >= 1.0)
			{
				vtkDebugMacro(<< "cosCAB= " << cosCAB << endl);
				return 0;
			}

			sinABC = std::sqrt(1.0 - cosABC * cosABC);
			sinBCA = std::sqrt(1.0 - cosBCA * cosBCA);
			sinCAB = std::sqrt(1.0 - cosCAB * cosCAB);

			if (sinABC < 1e-10)
			{
				vtkDebugMacro(<< "sinABC= " << sinABC << endl);
				return 0;
			}

			if (sinBCA < 1e-10)
			{
				vtkDebugMacro(<< "sinBCA= " << sinBCA << endl);
				return 0;
			}

			if (sinCAB < 1e-10)
			{
				vtkDebugMacro(<< "sinCAB= " << sinCAB << endl);
				return 0;
			}

			cotgABC = cosABC / sinABC;
			cotgBCA = cosBCA / sinBCA;
			cotgCAB = cosCAB / sinCAB;

			D(ptIdA, ptIdA) += cotgABC + cotgBCA;
			D(ptIdA, ptIdB) -= cotgBCA;
			D(ptIdA, ptIdC) -= cotgABC;

			D(ptIdB, ptIdB) += cotgBCA + cotgCAB;
			D(ptIdB, ptIdA) -= cotgBCA;
			D(ptIdB, ptIdC) -= cotgCAB;

			D(ptIdC, ptIdC) += cotgCAB + cotgABC;
			D(ptIdC, ptIdB) -= cotgCAB;
			D(ptIdC, ptIdA) -= cotgABC;

			cellIterCount++; 
		}

		vtkDebugMacro(<< "Stage 3 finish" << endl);

		VectorCoordType x(numberOfPoints, 0.0);
		VectorCoordType y(numberOfPoints, 0.0);
		{
			// solving Ax = b (D x = bx)
			VectorCoordType rx = bx;
			VectorCoordType zx(numberOfPoints);

			VectorCoordType ry = by;
			VectorCoordType zy(numberOfPoints);

			// Jacobi preconditioner
			VectorCoordType Dinv(numberOfPoints);
			for (unsigned int ip = 0; ip < numberOfPoints; ++ip) //PointIdentifier -> int
			{
				Dinv[ip] = 1.0 / (D(ip, ip) + DBL_MIN);

				zx[ip] = rx[ip] * Dinv[ip];
				zy[ip] = ry[ip] * Dinv[ip];
			}

			VectorCoordType dx = zx;
			VectorCoordType dy = zy;

			unsigned int numIter = bx.size();
			if (bx.size() != numberOfPoints)
			{
				// check for safe
				std::cerr << "bx.size() != numberOfPoints\n";
			}
			numIter += numIter / 10; // let the iteration times a little more than the
									 // dimension

			double tol = 1e-10;

			for (i = 0; i <= numIter; ++i)
			{
				VectorCoordType Dxd;
				D.pre_mult(dx, Dxd);
				VectorCoordType Dyd;
				D.pre_mult(dy, Dyd);

				double dDxd = inner_product(dx, Dxd);
				double dDyd = inner_product(dy, Dyd);

				double zxTrx = inner_product(zx, rx);
				double zyTry = inner_product(zy, ry);

				double alphax = zxTrx / (dDxd + DBL_MIN);
				double alphay = zyTry / (dDyd + DBL_MIN);

				x += alphax * dx;
				y += alphay * dy;

				rx -= alphax * Dxd;
				ry -= alphay * Dyd;

				double rxTrx = inner_product(rx, rx);
				double ryTry = inner_product(ry, ry);
				if (rxTrx < tol && ryTry < tol)
				{
					//      std::cout<<"out from here when i = "<<i<<std::endl;
					break;
				}

				for (unsigned int id = 0; id < numberOfPoints; ++id)
				{
					zx[id] = rx[id] * Dinv[id];
					zy[id] = ry[id] * Dinv[id];
				}

				double betaX = inner_product(zx, rx) / (zxTrx + DBL_MIN);
				double betaY = inner_product(zy, ry) / (zyTry + DBL_MIN);

				dx = zx + betaX * dx;
				dy = zy + betaY * dy;
			}
		}

		vtkDebugMacro(<< "Stage 4 finish" << endl);

		int outputPointIterator = 0;

		double point[3] = { 0.0, 0.0, 0.0 }; 

		double bounds[6];

		bounds[0] = std::numeric_limits< double >::max();
		bounds[1] = -std::numeric_limits< double >::max();

		bounds[2] = std::numeric_limits< double >::max();
		bounds[3] = -std::numeric_limits< double >::max();

		bounds[4] = std::numeric_limits< double >::max();
		bounds[5] = -std::numeric_limits< double >::max();

		if (this->m_MapToSphere)
		{
			vtkDebugMacro(<< "Map to sphere." << endl);
			if (m_MapScale < 0)
			{
				// < 0 means user doesn't explicitly assign it. Then
				// automatically calculate it s.t. after doing the
				// stereo-graphic projection, upper and lower hemi-sphere will have
				// same number of vertics.

				std::vector< double >           v_r2(numberOfPoints);
				std::vector< double >::iterator itv_r2 = v_r2.begin();

				for (i = 0; i < numberOfPoints; ++i, ++itv_r2)
				{
					*itv_r2 = x(i) * x(i) + y(i) * y(i);
				}

				std::sort(v_r2.begin(), v_r2.end());
				unsigned int uiMidPointIdx = 0;
				if (numberOfPoints % 2)
				{
					uiMidPointIdx = (numberOfPoints - 1) / 2;
				}
				else
				{
					uiMidPointIdx = numberOfPoints / 2;
				}
				this->m_MapScale = 1.0 / std::sqrt(v_r2[uiMidPointIdx]);
			}

			i = 0;
			while (outputPointIterator != numberOfPoints)
			{
				double xx = (this->m_MapScale) * x(i);
				double yy = (this->m_MapScale) * y(i);

				double radius2 = xx * xx + yy * yy;

				point[0] = 2.0 * xx / (1.0 + radius2);
				point[1] = 2.0 * yy / (1.0 + radius2);
				point[2] = 2.0 * radius2 / (1.0 + radius2) - 1.0;

				if (point[0] < bounds[0]) { bounds[0] = point[0]; }
				if (point[0] > bounds[1]) { bounds[1] = point[0]; }

				if (point[1] < bounds[2]) { bounds[2] = point[1]; }
				if (point[1] > bounds[3]) { bounds[3] = point[1]; }

				if (point[2] < bounds[4]) { bounds[4] = point[2]; }
				if (point[2] > bounds[5]) { bounds[5] = point[2]; }

				outPoints->InsertNextPoint(point);
				outputPointIterator++;
				i++;
			}
		}
		else
		{
			vtkDebugMacro(<< "Map to plane." << endl);
			i = 0;
			while (outputPointIterator != numberOfPoints) 
			{
				point[0] = x(i);
				point[1] = y(i);

				if (point[0] < bounds[0]) { bounds[0] = point[0]; }
				if (point[0] > bounds[1]) { bounds[1] = point[0]; }

				if (point[1] < bounds[2]) { bounds[2] = point[1]; }
				if (point[1] > bounds[3]) { bounds[3] = point[1]; }

				if (point[2] < bounds[4]) { bounds[4] = point[2]; }
				if (point[2] > bounds[5]) { bounds[5] = point[2]; }

				outPoints->InsertNextPoint(point);
				outputPointIterator++;
				i++;
			}
		}

		vtkDebugMacro(<< "Stage 5 finish" << endl);

		vtkDebugMacro(<< "bounds"
			<< " " << bounds[0] << " " << bounds[1]
			<< " " << bounds[2] << " " << bounds[3]
			<< " " << bounds[4] << " " << bounds[5] << endl);

		//Create duplicate references to the rest of data on the mesh


		//this->CopyInputMeshToOutputMeshPointData();
		//this->CopyInputMeshToOutputMeshCellLinks();
		//this->CopyInputMeshToOutputMeshCells();
		//this->CopyInputMeshToOutputMeshCellData();

		outputMesh->SetPoints(outPoints);
		outputMesh->SetPolys(inputMesh->GetPolys());

		output->DeepCopy(outputMesh);

		vtkDebugMacro(<< "Stage 6 finish" << endl);

		return 1;
	}
} // end namespace vtk

//#endif
