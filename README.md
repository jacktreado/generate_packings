# generate_packings

Code base to create packings compressed to a certain packing fraction $\phi$. 

## run test file
Compile/run `residue_main_test.cpp` with the following commands. From the directory with the `src/` and `examples/`directories:

`$ cd examples/test_main_files`

`$ g++ -I ../../src/ residue_pack_test.cpp ../../src/rigidbody.cpp ../../src/packing*.cpp ../../src/Quaternion.cpp -o -test.o`

`$ ./test.o`

Makefiles are on their way!
