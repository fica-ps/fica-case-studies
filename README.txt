

##Case Studies##

1)	To run the Case studies install Jupyter and Julia, follow the instructions:
	https://datatofish.com/add-julia-to-jupyter/

2)	Clone or download the files in https://github.com/fica-ps/fica-case-studies

3)	Open Jupyter and navigate to the src folder, open the desired Case Study

4)	The first block of each Case Study is commented, uncomment(remove the '#') and run it once, 
	this will add the necessary packages for each case study. After running it once
	you can don't need to run it again for that device.
	
##FICA##

1) Prerequisites
	The prerequisites to build fica are the following:

	cmake: version 13.0 or above.
	C++ compiler:
		GCC: version 4.8 or newer for Linux/MacOS.
		Microsoft Visual C++: 2012 or newer for Windows

2) Getting Started
	To build fica, first, clone the repository to a directory (https://github.com/fica-ps/fica)

3) Generating build configurations
	Open a terminal and use the following command:

	cmake PATH/TO/FICA -DCMAKE_BUILD_TYPE=Release

	doing so in a directory different from where fica sources were cloned is recommended

4.1) Building on Linux/MacOS
	To build the libfica Shared Object on Linux/MacOS just use the make command on the directory where cmake created the build configuration.

4.2)Building on Windows
	To build the fica Dynamic-link library on windows, open the solution on visual studio (fica.sln), then change the configuration to Release mode instead of Debug.

	Right click on the fica project on Solution Explorer and click on Build.

	This generates a subdirectory on the directory that contains fica.sln named Release, there you'll find fica.dll

More information and documentation can be found in the project Wiki: https://github.com/fica-ps/fica/wiki