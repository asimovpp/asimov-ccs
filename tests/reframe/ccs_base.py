#!/usr/bin/env python
"""ReFrame base module for Asimov-CCS tests"""

import os

import reframe as rfm
import reframe.utility.sanity as sn


class SetupDependencies(rfm.RegressionMixin):
  """Sets up the dependencies for build and run"""
  # https://reframe-hpc.readthedocs.io/en/stable/regression_test_api.html#reframe.core.pipeline.RegressionTestPlugin

  @run_before("compile")
  def prepare_build(self):
    """Setup environment for build"""
    # General options
    self.build_system.max_concurrency = 8
    self.build_system.build_dir = f"{self.stagedir}/ccs_build"
    self.build_system.options = ["ccs_app"]

    # select compiler
    if self.current_environ == "PrgEnv-Cray":
      cmp = "cray"
    elif self.current_environ == "PrgEnv-gnu":
      cmp = "gnu"
    else:
      cmp = "cray"
    print(self.current_environ)
    self.env_vars["CMP"] = cmp

    # local_vars
    install_dir="/work/e609/e609/shared/epcc/libs"
    adios2_v="2.10.2"
    parhip_v="3.19"
    petsc_v="3.23.7"
    hdf5_v="1.14.6"
    fyamlc_v="0.2.5"

    # Env vars
    python_dir = install_dir + "/python-" + cmp
    self.env_vars["PYTHONPATH"] = os.getenv("PYTHONPATH") + ":" + python_dir
    self.env_vars["PATH"] = os.getenv("PATH") + ":" + python_dir + "/bin"
    self.env_vars["PYTHONUSERBASE"] = python_dir

    # petsc (consider adding single precision)
    petsc_dir = install_dir + "/petsc-" + cmp + "-v" + petsc_v
    self.env_vars["PETSC"] = petsc_dir
    self.env_vars["PETSC_PRECISION"] = "double"
    self.env_vars["PETSC_ROOT"] = petsc_dir
    self.env_vars["PETSC_DIR"] = petsc_dir

    # FYAMLC
    fyamlc_dir = install_dir + "/fyaml-c-" + cmp + "-v" + fyamlc_v
    self.env_vars["FYAMLC"] = fyamlc_dir

    # ADIOS2
    adios2_dir = install_dir + "/adios2-" + cmp + "-v" + adios2_v
    self.env_vars["ADIOS2"] = adios2_dir

    # HDF5
    hdf5_dir = install_dir + "/hdf5-" + cmp + "-v" + hdf5_v
    self.env_vars["HDF5"] = hdf5_dir

    # PARHIP
    parhip_dir = install_dir + "/parhip-" + cmp + "-v" + parhip_v
    self.env_vars["PARHIP"] = parhip_dir

    # PARMETIS
    parmetis_dir = install_dir + "/parmetis-" + cmp
    self.env_vars["PARMETIS"] = parmetis_dir

    # PARMETIS_32BIT
    parmetis_32bit_dir = install_dir + "/parmetis-" + cmp + "-32bit"
    self.env_vars["PARMETIS_32bit"] = parmetis_32bit_dir

    # RCMF90
    rcmf90_dir = install_dir + "/rcm-f90-" + cmp
    self.env_vars["RCMF90"] = rcmf90_dir

    # MAKEDEPF90
    makedepf90_dir = install_dir + "/makedepf90-" + cmp
    self.env_vars["MAKEDEPF90"] = makedepf90_dir

    # update LD_LIBRARY_PATH
    self.env_vars["LD_LIBRARY_PATH"] = (
            petsc_dir + "/lib" + ":" +
            fyamlc_dir + "/lib" + ":" +
            parhip_dir + "/lib" + ":" +
            parmetis_dir + "/lib" + ":" +
            parmetis_32bit_dir + "/lib" + ":" +
            rcmf90_dir + "/lib" + ":" +
            makedepf90_dir + "/lib" + ":" +
            os.getenv("LD_LIBRARY_PATH")
            )

    # PETSC


class BuildCCS(rfm.CompileOnlyRegressionTest, SetupDependencies):
  """
  Compiles CCS, making it available as a reframe fixture:
  https://reframe-hpc.readthedocs.io/en/stable/tutorial.html#test-fixtures
  """

  build_system = "Make"
  modules = ["cray-python", "cmake"]
  #sourcesdir = "git@git.ecdf.ed.ac.uk:asimov/asimov-ccs.git"
  sourcepath = "asimov-ccs"
  #sourcesdir = "../../"
  local = True
  build_locally = False

  valid_prog_environs = ["PrgEnv-gnu", "PrgEnv-cray"]

  prebuild_cmds = ["git clone git@git.ecdf.ed.ac.uk:asimov/asimov-ccs.git", "cd asimov-ccs"]
  # prebuild_cmds = [
  #         "module list > ml.txt",
  #         "python -m venv --system-site-packages ./ccs-venv",
  #         "source ./ccs-venv/bin/activate",
  #         "python -m pip install lit",
  #         ]

  @sanity_function
  def sanity_executable_exists(self):
    """Confirms that the executable was built"""
    build_dir = f"{self.stagedir}/ccs_build"
    return sn.path_exists(os.path.join(build_dir, "ccs_app"))


@rfm.simple_test
class RunCCS(rfm.RunOnlyRegressionTest, SetupDependencies):
  """Run CCS"""

  valid_systems = ["archer2:compute"]
  valid_prog_environs = ["PrgEnv-cray"]
  executable = "ccs_app"
  stream_binary = fixture(BuildCCS, scope="environment")

  keep_files = ["residulas.log", "timers.csv"]

  strict_check = True
  tags = {"full_application"}


  @run_after("setup")
  def set_executable(self):
    """sets up executable"""
    self.executable = os.path.join(self.stream_binary.build_system.builddir, "ccs_app")


  @sanity_function
  def assert_finished(self):
    """Sanity check that simulation finished successfully"""
    file_len = len(self.keep_files[0])
    return sn.assert_ge(file_len, 10)
