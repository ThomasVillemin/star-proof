# Star-4V/S
This program statistically computes "An invariant property of diffuse random
walks" (S. Blanco and R. Fournier - Europhysics Letter 2003) on an user defined
geometry. The submitted geometry must define a set of polygonal meshes saved
with respect to the Alias Wavefront Obj fileformat. The resulting surfaces must
define a set of closed volumes whose normals point inward the volumes. The
program assumes that the shape surfaces are fully emissive and the diffuse
property of their volume is constant.

## How to build

Star-4V/S relies on [CMake](http://www.cmake.org) and the
[RCMake](https://gitlab.com/vaplv/rcmake/#tab-readme) package to build. It also
depends on the
[RSys](https://gitlab.com/vaplv/rsys/),
[Star-3D](https://gitlab.com/meso-star/star-3d/),
[Star-3DAW](https://gitlab.com/meso-star/star-3daw/),
[Star-MC](https://gitlab.com/meso-star/star-mc/) and the
[Star-SP](https://gitlab.com/meso-star/star-sp/) libraries.

First ensure that CMake is installed on your system. Then install the RCMake
package as well as all the aforementioned prerequisites. Then generate the
project from the `cmake/CMakeLists.txt` file by appending to the
`CMAKE_PREFIX_PATH` variable the install directory of its dependencies.

## Release notes

### Version 0.3

- Make the code compatible with Star-MC 0.4.

### Version 0.2

- Make the code compatible with Star-3D 0.4 and Star-SP 0.4.

## License

Star-4V/S is Copyright (C) |Meso|Star> 2015-2018 (<contact@meso-star.com>).
It is a free software released under the [OSI](http://opensource.org)-approved
CeCILL license. You are welcome to redistribute it under certain conditions;
refer to the COPYING files for details.

