# Overview

Asimov-ccs is a modular codebase where different simulation applications can be built up by putting together a set of components. This flexibilty requires the support of a sophisticated build system. 

The general approach is to compile all components available in ccs individually and link a subset of the components together to form the desired application. Not all combinations of components are valid. 

When a component can have multiple implementations, these implementations are put into different submodules. A particular implementation is selected at link time by linking in the compiled object of the desired submodule.

Note: the environment variable `CCS_DIR` needs to be defined and point to the root directory of the `asimov_ccs` project.

# Compilation 

All components can be compiled individually without concern of linking incompatabilites in the main application. Of course, the components themselves have some dependencies. These dependencies are found by `makedepf90` which reads all source files and parses Fortran `use` statements and module/submodule relationships. 

Some extra scripts are used to work around inflexibilities in Makefiles: 

- `filter-out.py` : filters out filenames from a list without taking the path into account

Additional information can be provided to the build system using doxygen build tags. For example, `!> @build mpi` or `!> @build petsc`. Multiple tags can be put on a single line: `!> @build mpi petsc`. This information could be used to know when to add compile flags for required libraries/headers, though this has not been added to the Makefile yet. 
Procesisng tags is handled by `process_build_tags.py`, which outputs a variable that can be used in a Makefile rule, e.g. `PETSC_OBJ = file1.o file2.o`.


# Linking

Once all object files are built, they need to be linked together according to a user specified configuration. This is mainly done by `generate_link_deps.py` that outputs a single Makefile rule that looks like `ccs_app: files.o to.o link.o`. 

First, some linking files are populated automatically. `process_dependencies.py` can read the build dependencies produced by `makedepf90` and process them into a graph. Nodes represent files and directed edges represent `use` relationships (e.g. if A contains `use B` then "A <- B" will appear in the graph). The graph can be visualised and analysed. The main analysis that is used in the link rule generation is finding the "internal" nodes of the graph which are here defined as nodes that have one or more outgoing edges, i.e. it the files are used by some other file. This set represents all nodes that are not submodules, programs and lone nodes. This set of files can be automatically added to the link rule because they do not have multiple possible implementations. 

Second, a config file is processed to add user chosen component implementations; this is done in `generate_link_deps.py`. The user config specifies a "main" (which is matched directly with a filename containing the Fortran `program`) and "options". Options can either be groupings of component implementation choices (e.g. base: mpi_petsc) or individual components to override earlier choices. The `config_mapping.yaml` file provides the mapping from user-facing config keywords to implementation filenames. Components that do not appear in `config_mapping.yaml` can be specified with `extra:` and will be added to the link line directly; single items or lists can be specified.

It is possible for a user to provide a configuration that cannot be linked or run. Such cases could be caught automatically at configure/link time given comprehensive compatability rules, but this is not implemented at the moment.


# Testing

Test building follows the same basic pattern as building ccs applications. All ccs files are built independently and then linked to provide the functionality that the test requires. However, an application would be specific in all components where there is a choice of implementation but a test would not. A test exercises the interfaces of components and has to be successful with all implementations of the components. Thus there is additionaly machinery to repeat test cases with all relevant submodule options.

The test suites are executed using the LLVM LIT framework.

## Running tests
`lit` is required to run tests. The easiest method is to run `pip install lit`.

Tests need to be run from the `src` directory in order to pass through the compilation flags that are used to build CCS. Execute `CMP=gnu make tests` (substitute gnu for your chosen compilation environment) to run the test suite. 
The expected output is compilation messages of CCS files (if not built already) followed by "-- Testing: n tests, m workers --" (where n and m are integers). Passing tests are listed briefly and failed tests provide their entire output. Further debugging can be done by examining temporary files in `Output` directories within the test subdirectories under `testing`.

## How tests are written
Test cases are written as single-file Fortran programs plus one or more configuration `yaml` files. 
The configuration files define every submodule combination that has to be tested for the given test program. 
The `main` in the configuration has to be equal to the name of the test program file (minus the `.f90` extension).

The configuration file needs to have one or more comments starting with `RUN:` that say how the test should be executed. 
The test case fails if any of the RUN commands fail.
`# RUN: %build_test %s %t1 mytest.f90` will compile the test case (identified with the macro `%s`) and store the executable in `t1`.
Note that the program file name (`mytest.f90` in this case) has to be specified as well; multiple files can be specified if there are, for example, utilities shared between multiple test programs.
This can then be executed with, for example `# RUN: %mpirun -n 4 %t1` (note, '%mpirun' is a macro that gets resolved to the mpi invocation command set in the used Makefile arch file).
If a compiled test case returns non-zero, the test case has failed.

Test cases can `use testing_lib` to included the testing library which contains various utility functions for test initialisation, finalisation, stopping and assertions (not everything has been implemented yet). 

## Testing framework details
The `lit.cfg` file configures the LIT testing framework. See therein or refer to the LIT user guide for details.

Tests are built using the `test_builder.py` script (this is what the `%build_test` macro invokes). 
It does the following actions:
- Parse the `yaml` file that defines the test case to find the test's source filename.
- Compile the test case object.
- Evaluate dependencies on CCS.
- Generate a link rule to the required CCS objects.
- Link the test case to create an executable.

## More resources
- LIT user guides
  - https://llvm.org/docs/CommandGuide/lit.html
  - https://llvm.org/docs/TestingGuide.html
- LIT python source
  - https://github.com/llvm/llvm-project/tree/main/llvm/utils/lit/lit
- LIT usage examples in LLVM
  - https://github.com/llvm/llvm-project/tree/main/llvm/test
  - https://github.com/llvm/llvm-project/blob/main/llvm/test/lit.cfg.py


# Linting
The code is linted primarily using Flint. Install via pip: `pip install flinter`.

Invoke the linter via the utility script using: `bash build_tools/lint.sh ACTION TARGET` where ACTION is "lint" or "score" and TARGET is the file or directory to lint. Specifying "lint" will print out detailed information on every discovered issue while "score" will return an overall numberic score (10 corresponds to "no issues").

Whitespace errors can largely be fixed automatically using fprettify. Install via pip: `pip install fprettify`.
The config file `.fprettify.rc` is set up to match as closely as possible the configuration of flint, and should be picked up automatically if `fprettify` is run from the root directory of ccs.

Use `fprettify FILE` to fix whitespace in FILE. To see the changes in stdout add `--stdout` or to see the diff add `--diff`; both flags prevent FILE from being changed.


