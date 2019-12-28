# vtkConformalFlatteningMesh

**Abstract**

This paper describes the Visualization Toolkit (VTK) Conformal Flattening Filter: $ConformalFlatteningFilter$. This VTK Polydata Algorithm is an implementation of a paper by Sigurd Angenent, et al., “On the Laplace-Beltrami Operator and Brain Surface Flattening” [1]. This filter performs an angle preserving map of any genus zero (i.e. no handles) surface to the sphere or, alternatively, to the plane. In this paper, we describe the code and provide the user with enough details to reproduce the results which we present in this paper. This filter has a variety of applications including the flattening of brain surfaces, which was the initial motivation for this work.



**User’s Guide**

The conformal flattening filter takes an VTK PolyData as input and will generate another VTK PolyData as output. The usage is basically the same as other PolyDataAlgorithms in VTK.

 

**Basic usage**

The filter is instantiated by:

 

*vtkSmartPointer filter =* *vtkSmartPointer::New();*

 

Then the input can be set and results can be obtained by:

 

*filter**->**SetInputData(inputMesh);*

and

*newPolyData = filter->GetOutput();*

 

**2.2 More about APIs**

The filter has several APIs for further manipulation of the output.


1. $setPointP$ function. On the right side of equation (1), the  $\delta p$  function depends on the location of the point $p$. Basically, this point will be mapped to infinity on the plane and the north-pole of the sphere. Hence the selection of the point $p$ determines which patch on the original mesh is mapped up to the north pole. 

 

The API setPointP takes an integer as input indicating the cell number in which the point $p$  lies. It’s a good choice to set the point $p$  where the local surface is relatively flat, i.e., having a small local curvature. If setting the point $p$ at some flat area is crucial, we suggest that user first use $vtkCurvatures$ to obtain the number of cells having low curvatures and then call this function using one of the cells with a low curvature.

2. The switch functions $mapToPlane$ and $mapToSphere$ determine the output to be either a plane or a sphere, the sphere being the default. The difference between the two mappings is simply a stereographic projection from the plane to the sphere. Simply by

 

*filter->mapToSphere( );*

or

*filter->mapToPlane( );*

 

users can switch between two different outputs.

![img](file:///C:/Users/GIGABYTE/AppData/Local/Temp/msohtmlclip1/01/clip_image006.jpg)

3. $setScale$ function. The mapping, calculated from the equation (1), is dependent on the number of the nodes within the mesh. Given a mesh of a large number of nodes and cells, the image of the flattening mapping, is constrained in a small range around origin. To make it cover a sufficient area of the plane and further get a reasonable result from stereographic projection, re-scale of the flattened plane is needed. This function is used to set the scale factor, by:

 

*filter->setScale( scale factor );*

 

For a mesh of around several thousands of nodes, a factor of 100 is likely to get a good result. The factor should grow as the number of nodes grows.
