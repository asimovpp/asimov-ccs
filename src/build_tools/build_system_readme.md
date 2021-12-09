# Overview

Asimov-ccs is a modular codebase where different simulation applications can be built up by putting together a set of components. This flexibilty requires the support of a sophisticated build system. 

The general approach is to compile all components available in ccs individually and link a subset of the components together to form the desired application. Not all combinations of components are valid. 

When a component can have multiple implementations, these implementations are put into different submodules. A particular implementation is selected at link time by linking in the compiled object of the desired submodule.


# Compilation 

All components can be compiled individually without concern of linking incompatabilites in the main application. Of course, the components themselves have some dependencies. These dependencies are found by `makedepf90` which reads all source files and parses Fortran `use` statements and module/submodule relationships. 

Some extra scripts are used to work around inflexibilities in Makefiles: 

- `filter-out.py` : filters out filenames from a list without taking the path into account
- `fix_makefile.py` : brute force fix to adding the compile command to all file dependency rules 

Additional information can be provided to the build system using doxygen build tags. For example, `!> @build mpi` or `!> @build petsc`. Multiple tags can be put on a single line: `!> @build mpi petsc`. This information could be used to know when to add compile flags for required libraries/headers, though this has not been added to the Makefile yet. 
Procesisng tags is handled by `process_build_tags.py`, which outputs a variable that can be used in a Makefile rule, e.g. `PETSC_OBJ = file1.o file2.o`.


# Linking

Once all object files are built, they need to be linked together according to a user specified configuration. This is mainly done by `generate_link_deps.py` that outputs a single Makefile rule that looks like `ccs_app: files.o to.o link.o`. 

First, some linking files are populated automatically. `process_dependencies.py` can read the build dependencies produced by `makedepf90` and process them into a graph. Nodes represent files and directed edges represent `use` relationships (e.g. if A contains `use B` then "A <- B" will appear in the graph). The graph can be visualised and analysed. The main analysis that is used in the link rule generation is finding the "internal" nodes of the graph which are here defined as nodes that have one or more outgoing edges, i.e. it the files are used by some other file. This set represents all nodes that are not submodules, programs and lone nodes. This set of files can be automatically added to the link rule because they do not have multiple possible implementations. 

Second, a config file is processed to add user chosen component implementations; this is done in `generate_link_deps.py`. The user config specifies a "main" (which is matched directly with a filename containing the Fortran `program`) and "options". Options can either be groupings of component implementation choices (e.g. parallel: poisson) or individual components to override earlier choices. The `config_mapping.yaml` file provides the mapping from user-facing config keywords to implementation filenames.

It is possible for a user to provide a configuration that cannot be linked or run. Such cases could be caught automatically at configure/link time given comprehensive compatability rules, but this is not implemented at the moment.


# Testing

Test building follows the same basic pattern as building ccs applications. All ccs files are built independently and then linked to provide the functionality that the test requires. However, an application would be specific in all components where there is a choice of implementation but a test would not. A test exercises the interfaces of components and has to be successful with all implementations of the components. Thus there is additionaly machinery to repeat test cases with all submodule options.
