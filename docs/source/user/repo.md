# Structure of the repository

In order to facilitate a quicker comprehension of the project, a structure of
the repository is provided:

```
pas-cg                       Repository
├── cg                       Source code
|   ├── data                 Data provided along with the program
|   |   ├── 1ubq             1ubq example files
|   |   ├── 9aac             9aac example files
|   |   ├── default          Default parameters
|   |   ├── extras           Extra files, for the most part M-J matrices
|   |   └── glut             glut example files
|   ├── ext                  External libraries
|   ├── include              Headers
|   |   └── cg
|   ├── scripts
|   ├── src                  Sources
|   └── CMakeLists.txt
├── docs                     Sources for the documentation
├── tests
├── reference                Reference (Fortran 77) implementation
|   ├── cg                   Code for the reference implementation
|   └── papers               Some research papers on the model
|       └── CPC14.pdf        "Main" (?) research paper
├── .clang-format            Settings for clang-format used in the sources
├── CMakeLists.txt
├── LICENSE
└── README.md
```

It is important to take a look at the inputfiles, such
as `cg/data/default/inputfile.yml` and other associated files, in case the
explanations in this documentation were unclear. 