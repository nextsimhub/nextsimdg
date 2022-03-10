# Modules
The module system provides a way of selecting the configuration of NextSIM_DG at runtime, without having to re-compile the model to obtain different behaviour. It often operates along with the configuration system to allow selection of the particular implementation using the text config files. The module system provides a method of decoupling the selection of a particular implementation from its point of use.

## Interfaces and implementations
The module system works with classes defined using typical C++ polymorphic classes. An abstract interface class defines the visible behaviour of the module and one or more implementation classes derived from it. No further code is required in the interface or implementation classes to use the module system. However, and addition Module class is required to provide the necessary infrastructure.

## `Module` namespace and class
The Module infrastructure provides a pair of templated functions `getInstance<T>()` and `getImplementation<T>()` in the `Module` namespace. By calling specializations of these functions, an instance of the  the define implementation class can be obtained. For a new instance `getInstance<T>()` is called which returns a `std::unique_ptr` of type of the interface class, pointing at a new instance of the implementation class. The `getImplementation<T>()` function returns a reference to a static instance of the implementation class, which can be used when a new instance is not needed. This could be used for pseudo-static functions, for instance.

The implementation is set using the `setImplementation<T>()` function, which takes a string argument of the name of the implementaion to be used, allowing the implementation to be set from the text of a configuration file.

## Implementing a new `IInterfaceModule` class

In the following description, the interface class is named `IInterface` and the implementation class is `Implmentation`.

### The header
The header contains the declaration of the Module specialization for `IInterface`. It will contain an include of the `Module` header and of the `IInterface` header. Within the `Module` namespace there is the explicit specialization of the `functionMap` variable.
```
template <> Module<IInterface>::map Module<IInterface>::functionMap;
```

There is also the declaration of the Module class itself. This consists mostly of boilerplate code to set up the module class.
```
class IInterfaceModule : public Module<IInterface> {
    struct Constructor {
        Constructor();
    };
    static Constructor ctor;
};
```

### The source
The source file for the Module class consists of template specializations to relate the names of the implementations to the classes themselves. The code for the source file can be generated using the Python 3 script `module_maker.py`. This takes an `--interface` argument for the interface class name and then a list of the implementation classes. The first implementation in the list is used as the default which will be returned if no other implementation is selected. The text of the source file is then returned to standard output.

To implement the source of the module without using the python script, it is best to base it off one of the existing source files.

## Using the module
### Selecting an implementation
When the module is constructed a default implmentation is loaded and will be returned if no other implementation is selected. Otherwise, one of the other implementations will be selected when named as an argument to `setImplementation()`. If the argument to `setImplementation()` is not the name of an implementation of the interface then a `std::out_of_range` exception will be throw, based on the underlying `std::map`.

How the implementation selecting string is selected is entirely up to the developer. The `boost::program_options` based classes `Configurator` and `ConfiguredModule` may be of use.

### Using the interface
The usual way of using a module is to either acquire a fresh instance using `getInstance()` or by using the static instance that can be referenced through `getImplementation()`. The instance returned by `getInstance()` is referred to by the `std::unique_ptr` that is returned, which removes the need for disposing of it when its use is complete.

Either through the reference or the smart pointer, the implementation can be accessed using any of the functions defined in the interface class, just as with any other C++ polymorphic class.
