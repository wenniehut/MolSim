MolSim GROUP D
===



## Instructions
### Build and Run
`mkdir build`<br>
`cd build`<br>
`cmake .. [-DBUILD_DOCUMENTATION=ON|OFF]`<br>
`make`<br>
`./MolSim <path_to_input_file>`

### Generate Doxygen Documentations
`make doc_doxygen`<br>

Note: Only possible if the project was built with the option
`-DBUILD_DOCUMENTATION=ON`