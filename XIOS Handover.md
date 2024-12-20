## Summary

#### Deliverable

A version of nextSIM-DG using the XIOS framework for model I/O.

#### Request

* Expertise in I/O and parallelisation.
	* Key aspects for eventually coupling neXtSIM to a climate model (a major deliverable of SASIP) and running it at high resolutions (which is of enormous scientific interest).
* Review PRs, discuss development strategies and decisions.

#### Requirements

1. XIOS must not become a compulsory dependency at this stage and, in the event that a user does not link against XIOS during the build, we revert back to the previous functionality.
2. XIOS must be able to create nextSIM-DG outputs – a diagnostic file and a restart file – via methods which leverage XIOS and match current outputs. Production of the log file will be unchanged.
3. XIOS must be enabled/disabled by a parameter in the config file under its own heading.
4. XIOS data must be managed by calls to methods provided in the `nextsim` namespace.

## Background

#### Current I/O approach

* Here we restrict attention to the `ParametricGrid` class is used to represent the mesh (other options are available).
* nextSIM-DG currently reads in information and outputs data to/from NetCDF files using the NetCDF library `netcdf-cxx4`.
	* Read occurs during model initialisation.
	* Write occur at the end of each dynamics step and at the end of the simulation.
* Core class for I/O: `ParaGridIO`, defined in `core/src/include/ParaGridIO.hpp` and `core/src/ParaGridIO/cpp`.
	* Reading driven by `getModelState` method of `ParaGridIO`.
	* Writing driven by `dumpModelState` method of `ParaGridIO`.
	* The `ParaGridIO` class also implements other methods supporting these, as well as for other I/O related functionality.
* ...

## Design

* **Configuration**
	* XIOS typically uses XML files to configure its objects and their attributes.
	* We would like to avoid this as far as possible and instead use a combination of nextSIM-DG's existing configuration files and calls to XIOS' API.
	* We currently still use a minimal XML file in `core/test/iodef.xml` for setting up the XIOS context.
		* This gets symlinked to the build directory as part of the build process.
* **XIOS handler class**
	* Implementing XIOS functionality in a separate I/O handler class allows the code's purpose to be clearly defined and prevents bloat.
	* Access to XIOS functionality should ordinarily be via the existing methods and a flag in the config file (or failure to build with the XIOS library) should determine if read/write is to be done via XIOS.
	* The XIOS handler class will be in `core`. (`dynamics` and `physics` should not depend on this class.)
	* The XIOS handler class will be created and maintained by the model (or an adjacent) class.
	* The XIOS handler class will need to be able to process information regarding the grid-processor distribution and creates a dependency on the MPI feature.
	* The XIOS handler class `Xios` is defined in `core/src/include/Xios.hpp` and `core/src/Xios.cpp`.
* **Integration**
	* Interfacing with the XIOS library and XIOS calls will be made by a handler class which will be owned by the grid. Sever initialisation will be handled during XIOS class creation.
	* XIOS state will be handled (approximately) when the `ModelState` is initialised.
	* XIOS state will need to be periodically updated.
* **Config**
	* Users will be able to enable XIOS via an `ENABLE_XIOS` flag in the config file, with the default being that it isn't included.
* **Calendar**
	* Various setters and getters have been created to handle the calendar properties. This removes the need for any pre-processing steps to ensure that the `iodef.xml` file matches the content of the neXtSIM-DG config file.
	* Before accessing any calendar information, it must be set either using these methods or within the `iodef.xml` file (or children thereof).
		* We are advocating for only using the methods of the `Xios` class.
	* Some of the methods in an earlier implementation raised XIOS C-exceptions, which were undocumented and would require debugging in a third-party library.
		* As such, it was decided to create bespoke converters for creating XIOS data types `cxios_date` and `cxios_duration` (which behave like structs).
			* The cost of maintaining this is not significant, since we would need to convert date-time standards between XIOS and neXtSIM anyway (XIOS drops the `T` and the `Z` from the [ISO 8601](https://www.iso.org/iso-8601-date-and-time-format.html) standard).
* **Calendar duration**
	* nextSIM-DG uses the [ISO 8601](https://www.iso.org/iso-8601-date-and-time-format.html) duration standard, except that it uses the convention `PY-DTh:m:s`, where `D` can run to 365 and `h`, `m`, and `s` are *any* integers, including those larger than 24, 60, and 60, respectively. No additional padding is required where an entry is zero, e.g. 10 hours can be `P0-0T10:0:0`. Use of the format `PY-M-DTh:m:s` is allowed but `1M` is interpreted as 30 days and is converted to seconds internally for all calculations.
	* We may wish to throw a warning in this scenario and we will need to translate all `M` entries to 30 days to advance duration in the XIOS calendar.
* **Tests**
	* We want tests which verify:
		1. that the XIOS library has been installed correctly by the user;
		2. that XIOS doesn't produce errors at runtime;
		3. that our usage of XIOS is and continues to be correct and performing functionality/behaviours as intended.
	* *Unit tests*
		* XIOS C++ demos run inside `doctest` as part of the nextSIM-DG test suite.
			* Does it run? Can the XIOS server be initialised with 1, 2, or 3 MPI processors?
			* Can we recreate a 'golden' output file / 'known good output (KGO)' for a demo grid? (Requirement 2)
			* Can we modify data and query it? (Requirement 4)
	* *Integration tests*
		* Feature control/access switch
			* Can we control usage of XIOS via switch(es) in the config file? Do outputs match in both cases? (Requirement 3)
			* Can we compile XIOS when not linked during build? (Requirement 1)
		* Dump model data
			* Can we run the model for each grid type and produce the 'golden' output file (`restart.nc`)? (Requirement 2)

## Work so far

* Background investigation (Alex) [point to issue on this]
* Implementation of `FindXios.cmake` file (Alex).
* Wrapping core XIOS functions in `core/src/include/xios_c_interface.hpp`, `core/src/Xios.hpp`, and `core/src/Xios.cpp` (Tom).
* Tests for core XIOS functionality in `core/test/XiosInit_test.cpp` (Tom).
* Separation of concerns of `core/test/XiosInit_test.cpp` into `core/test/Xios<Calendar/Axis/Domain/Grid/Field/File/Write>_test.cpp` (Joe).
* Wrapping XIOS reading functionality, `core/test/XiosRead_test.cpp` test (Joe).
* Wrapping other misc XIOS functionality (Joe).
* Use nextSIM-DG `TimePoint` and `Duration` classes in XIOS interface (Joe).
* Methods for converting between nextSIM-DG and XIOS' time and duration concepts (Joe).
* ...

## Still to do

* ...

## Building nextSIM-DG with XIOS support

* XIOS has not been added on the `main` branch of nextSIM-DG, so you need to branch off of `develop`.
* [Install XIOS2]
* [Install XIOS3]
* [Symlink]
* [Set `${xios_DIR}`]
* [Spack env for nextSIM-DG, although might want to ask Tim and Einar how they have things set up on Mac]

Make a copy of the `setup-joe.sh` script in the nextSIM-DG base directory:
```sh
cp setup-joe.sh setup.sh
```
Run
```sh
source setup.sh
```
and check that your prompt gets prepended by `[nextsim] (nextsim)`, denoting that the Spack and Python virtual environments have been activated, respectively.

Create a copy of the `build-joe.sh` script, too:
```sh
cp build-joe.sh build.sh
```
[TODO: Add help; debug mode vs release mode]
Run
```sh
./build.sh --help
```
to see the command line options that the build script accepts. Then try building the model with your desired options.

## Useful resources

* There is a [toy code](https://github.com/nextsimhub/xios_cpp_toy/blob/master/src/test_tp3.cpp) under the `nextsimhub` GitHub organisation, written by one of the XIOS developers. It demonstrates how to use some basic XIOS functionality in C++.
* [XIOS Mattermost + SASIP channel]
* [neXtSIM paper](www.the-cryosphere.net/10/1055/2016) for the model that nextSIM-DG is a re-write of.
* [nextSIM-DG preprint](https://doi.org/10.5194/egusphere-2023-391) for a paper on nextSIM-DG itself.
* [Robert's GPU port paper]
* [Original ICCS resource allocation](https://github.com/Cambridge-ICCS/ResourceAllocation/issues/28#event-12262163234)
* [neXtSIM-DG GitHub](https://github.com/nextsimhub/nextsimdg)
* [domain_decomp GitHub](https://github.com/nextsimhub/domain_decomp)
* [nextSIM-DG development practices]


## Tasks

#### Low-hanging fruit

* Change the number of MPI processors for some of the XIOS unit tests so that we cover more cases than just 2 MPI ranks. e.g., `XiosCalendar_test` could use 1 MPI rank and one of the others could be modified to run with 3 MPI ranks.

#### Major tasks

* ...

#### Other tasks

- Ensure the conversion of `1M` (1 month) to 30 days is documented on the XIOS pages and the nextSIM-DG time/config file pages