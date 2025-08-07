# MathLib
MathLib/
├── include/         # Vector3D.h
├── src/             # Vector3D.cpp
├── tests/           # test_vector.cpp
├── CMakeLists.txt
└── README.md
# From root folder
mkdir build && cd build
cmake ..
cmake --build .
ctest  # or run MathLibTests executable
