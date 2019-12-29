#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredGrid.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <string>
#include <vtkPolyDataWriter.h>
#include "vtkConformalFlatteningMeshFilter.h"


int main(int argc, char *argv[])
{
	 //Ensure a filename was specified
	if (argc != 4)
	{
		std::cerr << "Usage: " << argv[0] << " InputFilename OutputFilename Operation" << endl;
		return EXIT_FAILURE;
	}

	// Get the filename from the command line
	std::string inputFilename = argv[1];
	std::string outputFilename = argv[2];
	std::string operation = argv[3];

	// Get all data from the file
	
	vtkSmartPointer<vtkGenericDataObjectReader> reader =
		vtkSmartPointer<vtkGenericDataObjectReader>::New();
	reader->SetFileName(inputFilename.c_str()); 
	reader->Update();

	vtkPolyData* inputMesh = reader->GetPolyDataOutput();
	cout << "Get input file" << endl;

	vtkSmartPointer<vtkPolyDataWriter> writer =
		vtkSmartPointer<vtkPolyDataWriter>::New();
	writer->SetFileName(outputFilename.c_str());

	vtkSmartPointer<vtk::ConformalFlatteningMeshFilter> filter =
		vtkSmartPointer<vtk::ConformalFlatteningMeshFilter>::New();
	if (operation == "0")
		filter->MapToPlane();
	else if (operation == "1")
		filter->MapToSphere();
	else {
		std::cerr << "Unknown operation type: " << operation << endl;
		std::cerr << "Try 1 for sphere or 0 for plane." << endl;
		return EXIT_FAILURE;
	}
	filter->DebugOn();
	filter->SetInputData(inputMesh);
	filter->Update();
	
	cout << "Update success." << endl;

	writer->SetInputData(filter->GetOutput());

	// Optional - set the mode. The default is binary.
	//writer->SetDataModeToAscii();

	writer->Write();
	cout << "Writing finish." << endl;

	return EXIT_SUCCESS;
}