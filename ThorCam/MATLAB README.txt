MATLAB README

For MATLAB support for all Thorlabs DCC- and DCU-series cameras, please see the DCx Camera Support folder typically installed at C:\Program Files\Thorlabs\Scientific Imaging\DCx Camera Support. MATLAB support is provided through the .NET interface.

For MATLAB support for all Thorlabs Scientific Cameras (340, 1500, 1501, 4070, 8050, 8051, CS2100M-USB, CS235MU, CS235CU, CS505MU, CS505CU, CS505MUP, CS895MU, and CS895CU), please use the .NET camera interface by following these directions:

1. Install ThorCam and the appropriate drivers using the installer CD that came with a camera or download an installer from the appropriate product page at Thorlabs.com.

2. Copy the managed DLLs from

 Scientific Camera Interfaces\SDK\DotNet Toolkit\dlls\Managed_64_lib\*.dll 
 
 to the same folder with your MATLAB .m files.

3. See the following guides found in the Documentation folder (usually C:\Program Files\Thorlabs\Scientific Imaging\Documentation\Scientific Camera Documents):
 
 TSI_Camera_MATLAB_Interface_Guide.pdf
 TSI_Camera_DotNET-LabVIEW-MATLAB_Programming_Guide.chm
 