# Simple c++ rsa algo realization. 
# Requirement  
Require boost cpp_int and boost random. Can be easily changed to different long ints and random
# Build
## Linux
`g++ rsa.cpp -lboost_random -lpthread`
## Windows
For MinGW: `-lboost_random-mgw12-mt-x64-1_81`
Note that it should be your boost version.