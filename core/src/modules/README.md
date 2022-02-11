# ModuleLoader

Defining implementations of an interface to load.

## Brief
Given multiple implementations of a given function, how can we choose at the start of run time which implementation to use? For example, there might be three different ways of calculating sea-ice albedo. Only one will be used in a particular model run, but we don't want to have to recompile the model just to change to a different implementation.

The alternative implementations inherit from an interface class. The specific implementation is then accessed through a base class `std::unique_ptr` and executes the chosen implementation using virtual functions.

The configuration comes from a JSON module configuration file. The selection of the implementation of each module is done on initialization of the module loader. Access to the loaded module is provided by templated `getImplementation<T>()` functions, where the function parameter `T` is the interface class.

## Usage
To use the ModuleLoader module loader, there must be a set of interfaces and implementations of those interfaces. These are coded in the usual way, but the interfaces must have header files named in the pattern `Iheader.hpp`, where `header` is replaced by the name of the interface, as specified below.

In the code that uses ModuleLoader, new instances of the chosen implementation of the interface are obtained by calling `getImplementation<T>`, where `T` is the name of the interface class, including the leading 'I', as in the header file name.

For example, if the interface is named "Albedo" in the configuration, then the header file will be called `IAlbedo.hpp` and an instance of the implementing class will be pointed to by the return value of `loader.getImplementation<IAlbedo>()`.

The configuration of the interfaces and implementations is specified in a JSON file. By default this is named `modules.json`, but this name can be changed in the `CMakeLists.txt` file by changing the value of the CMake variable `ModuleLoaderFile`.

The JSON file should specify an array. Each element is an object with two members: `name` and `implementations`. The value of `name` is that of the interface. In the example above this would be the string "Albedo". The member  `implementations` is an array of the names of the classes that implement this particular interface. These will be found in the header that is their name suffixed with `.hpp` (with no additional affixes). In that header the class is derived from the `Iinterface` class.

### Example
Using alternative implementations of albedo as an example, the JSON file looks like:

    [
        {
        "name": "Albedo",
        "implementations": [
            "SnowAlbedo",
            "IceAlbedo",
            ]
        },
    ]

The interface class is declared as

    class IAlbedo {
 in the header file `IAlbedo.hpp`. The two implementations would be declared as

    class SnowAlbedo: public IAlbedo {
in `SnowAlbedo.hpp` and

    class IceAlbedo: public IAlbedo {
in `IceAlbedo.hpp`, with these classes defined somewhere in the source tree that makes up the overall executable. The ModuleLoader does not care about implementation source file naming.

Within the code, implementations for the interfaces are selected using the `ModuleLoader::setImplementation()` function, which takes the names as string arguments.

```
    loader.setImplementation("IAlbedo", "IceAlbedo");
```
Once this is done, objects of the implementing class can be accessed using templated functions either as a static instance that is stored within the ModuleLoader class

```
    IAlbedo& alb = loader.getImplementation<IAlbedo>();
```

or as a new `std::unique_ptr`.

```
    std::unique_ptr<IAlbedo> pAlb = std::move(loader.getInstance<IAlbedo>();
```

### Building
The module file is passed to a Python script that will parse the JSON and produce the inclusion files that are required to implement the ModuleLoader class with the chosen sets of interfaces and implementations. The command line will look like

    python moduleloader_builder.py [--ipp path_prefix] [--hpp path_prefix] module_file
The `--ipp` option will apply the specified text directly in front of the defined `.ipp` file names. This will usually be a directory (including trailing separator) where the files are to be found, but could include a file name prefix if desired, and the names in `ModuleLoader.cpp` are edited to match.

The `--hpp` option will apply the specified text directly in front of the generated header file names in the interface and implementation classes include statements. This will usually be a directory (including trailing separator) where the files are to be found, but could include a file name prefix if desired, as long as the header files have the same prefix.

The module_file argument is the path to the file specifying the interface and implementation names, as specified above. If no name is given, the default file `modules.json` will be tried to be read.

### CMake integration
The builder system is designed to be integrated into a CMake build process. This is done by making the executable target depend on a custom target that runs the Python script. An example can be found in `modules/CMakeLists.txt` in this repository, but the custom target is based on information from a StackOverflow [answer](https://stackoverflow.com/a/49021383). This set-up allows the modules to be built with minimal intrusion into the main build process. The subsidiary module loader CMake file `modules/CMakeLists.txt` requires the CMake variable `ModuleLoaderIppTargetDirectory` to be defined. This should be any directory defined as an include directory for the executable, including a trailing directory separator.
