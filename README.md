[![Doxygen -> gh-pages](https://github.com/nextsimdg/nextsimdg/actions/workflows/doxygen.yml/badge.svg)](https://nextsimdg.github.io/nextsimdg/index.html) 
[![Documentation Status](https://readthedocs.org/projects/nextsim-dg/badge/?version=latest)](https://nextsim-dg.readthedocs.io/en/latest/?badge=latest)

# neXtSIM_DG : next generation sea-ice model with DG

Here you will find the sources for neXtSIM_DG sea-ice model.

Below are some conventions that every contributors to this model must apply in order to have the most efficient and consistent workflow for the model developments.

## Coding conventions

For neXtSIM_DG we use clang-format and the [Webkit style](https://webkit.org/code-style-guidelines/), with a 100 character line length limit.
  -  neXtSIM is written using ISO C++11
  -  All array operations should be done using std::vectors - not C-style arrays
  -  The use of C-style pointers, new, and delete is strongly discouraged
  -  Names and values of physical constants reside in src/include/constants.hpp
  -  Runtime modifiable model parameters are set in src/options.cpp and can be set using option files.


If your code editor does not have a clang formatting plugin, you can use  the git clang-format command before pushing your code (see Section 3 for git conventions) 

## Commenting conventions for a nice automatic documentation

The automatic generation of documentation is done with Doxygen and Sphinx. The update of the documentation each time the code is modified is done with readthedocs and github actions.

In order to have the most understandable and consistent documentation of the API, we expect that contributors to neXtSIM write some comments alongside their code for every C++ component (namespaces, types, methods, functions, and constants) they implement.

Every file must begin with :
  -  a description of what it contains
  -  the name of the authors
  -  the date of creation/modification

Example :

```
/*!
* @file   gridoutput.hpp
* @author Einar Olason <einar.olason@nersc.no>
* @date   Thu Aug  4 09:47:27 CEST 2016
*/

```

We ask every coders to follow the following commenting conventions for comments providing the automatic documentation:
  -  documentation blocks must appear before the component it describes, with the same indentation (classes and functions),
  -  comments that should appear in the automatic documentation start with ```/*!``` or ```//!```
  -  multi-line documentation block must begin with ```/*!``` and end with ```*/``` on their own line,
  -  single-line block begins with ```//!```,
  -  description of variables or constants on the same line they are declared must preceded by ```//!<```,
  -  documentation must use tags like ```@see``` or ```@param``` and Markdown for formatting.

Comments starting with ```//``` or ```/*``` won’t appear in the automatic documentation so they can be used for internal or development comments.

The structure of a documentation block for classes and functions should be organized as follow :
  -  short summary
  -  extended summary (optional)
  -  template parameters (for classes, methods and functions)
  -  function/method parameters (for methods and functions)
  -  returns (for methods and functions)
  -  throws (for methods and functions)
  -  exception safety (for methods and functions)
  -  helper functions (for functions)
  -  initializer declaration (for constants)
  -  see also (optional)
  -  notes (optional)
  -  references (optional)
  -  examples (optional)

Examples :

```
/*!
 * Sum numbers in a vector.
 *
 * @param values Container whose values are summed.
 * @return sum of `values`, or 0.0 if `values` is empty.
 */
double sum(std::vector<double> & const values) {
    ...
}
```

```
/*!
 * Compute an element-wise cosine.
 *
 * @see sin
 * @see tan
 * @see [numpy.vectorize](https://docs.scipy.org/doc/numpy/reference/generated/numpy.vectorize.html)
 */
vector<double> cos(vector<double> const & angles);
```

```
/*!
 * This is an amazing function! For example:
 *
 *     auto cosines = cos(angles);
 *
 * Comment explaining the second example.
 *
 *     auto cosines = cos(radians(angles));
 */
```

```
/*!
 * Supported coordinate systems for flux-conserving transformations.
 *
 * These values are used in arguments to various functions in this package.
 * Unless otherwise stated, these functions do not validate whether the data
 * set makes sense in the "from" coordinates.
 */
enum class CoordType
{
    PIXEL,    //!< Untransformed detector coordinates
    WARP_PIXEL,    //!< Distortion-corrected detector coordinates
    SKY_WCS    //!< Equatorial J2000.0 coordinates
};
```

```
/*!
 * Read an image from disk.
 *
 * @param fileName the file to read. Must be either absolute or relative to
 *     the program working directory.
 *
 * @return the image stored in `fileName`. If the image on disk does not
 *     have `double` pixels, they will be cast to `double`.
 *
 * @throws IoError Thrown if `fileName` does not exist or is not readable.
 *
 * @exceptsafe Strong exception guarantee.
 */
lsst::afw::image::Image<double> loadImage(std::string const & fileName);
```

The comments providing automatic documentation should be supplemented with ordinary explanatory and internal comments using the standard comment syntax of ```/*``` and ```//```. These will not be picked up by the automatic documentation and can take a more general format.

## Code managment on github


### Version numbering

We use [semantic versioning](https://semver.org/). In brief this means the main branch has a version number assigned (tagged) to each commit. The numbers are of the form major.minor.patch, where:

1. MAJOR version when you make incompatible API changes
2. MINOR version when you add functionality in a backwards-compatible manner, and
3. PATCH version when you make backwards-compatible bug fixes.

This is not directly applicable to our workflow, but changes in major numbers should be related to major user facing changes (different input or output format, for instance), while minor numbers should (mostly) relate to changes in functionality (new physics, for instance). The patch number is incremented for each hotfix (see below).

### Git branching and merging

For the code we adopt the following system for branching and merging:

1. **main** branch: numbered releases of the code. Never edited. Merged from *develop* and *hot fix* branches (see notes on workflow below). Long living.
2. **develop** branch: rather stable version of code under development. Never edited. Merged from topic specific issue branches. Long living.
3. **feature**```<NNN>_<short_heading>``` : feature branches to develop significant components. Related to one or more issues, which should be linked to the relevant draft pull request (NNN = draft pull request number). Long living. Branched from, regularly updated from and merged back into develop.
4. **issue**```<NNN>_<short-heading>``` : issue specific branches (NNN = issue number). Main working area. Short living. Branched from, and merged back into develop or a feature branch.
5. **hotfix**```<NNN>_<short-heading>``` : branches that are specific to a hotfix issue. Hotfixes are bugfixes on main that should be fixed as soon as possible.

Note :
1. Never edit code in the main or develop branch. Always make a new branch for your edits.
2. A new branch should be very specific to only one problem. It should be short lived.
3. Commit often.
4. Branch often.
5. Branch only from main (hotfixes), develop (features and issues) or feature branches (feature related-issues).
6. Create pull requests for your branches and always assign a reviewer to merge and delete the branch, and close the issue. Commit messages should be formatted as a short description, followed by an empty line, followed by a detailed description. Line wrapping at 80 characters is preferred.
7. Every commit must compile.
  - With one exception: If it is necessary to commit code that does not compile, for example as part of a large refactor, then the commit message must be prefixed by WIP (for Work In Progress). This will make it clear to other developers that it is a commit that does not compile.

  
### How to report and handle new issues (bugs, improvements, new features, etc.)

If you discover a bug check that no one else has reported the same issue on https://github.com/nextsimdg/nextsimdg/issues. If the bug has not been reported, create an issue and assign someone to fix it (possibly yourself). Please notify people in person if you're assigning issues to them.

If you would like to suggest an improvement or a new feature check that no one else has made a similar request on https://github.com/nextsimdg/nextsimdg/issues. If this is not the case create a new issue, assigning or mentioning anyone you think could be affected or interested by your suggestion.

If you have been assigned an issue on https://github.com/nextsimdg/nextsimdf/issues address it using the following steps. For issues requiring a hotfix (bugs on the main branch that should be fixed as soon as possible):

1. Branching off from main, create an issue branch on your local system named **hotfix**```<NNN>_<short-heading>``` where NNN is the issue number from GitHub. This will be the main (short living) working area.
2. Write the necessary code.
    a. Make sure you test your modifications well. 
    b. Feel free to commit and push your issue branch often
3. Once the issue is fixed merge the **main** branch back into your issue branch and resolve any conflicts.
4. Create a pull request on GitHub to merge the issue branch back into **main**,  always include at least one reviewer who will then merge and delete the issue branch.
5. Merge the **main** branch into **develop** in your local repository, resolve conflicts, test, and push the updated **develop** branch.
6. Tag the merge by incrementing the patch number of the version number.
    a. Do this on the command line with ```git tag LABEL```.
    b. See the existing tags with ```git tag -l```.
7. Close the issue.

For issues not requiring a hotfix (less urgent bug-fixes and feature requests):

1. Branching off from **develop**, create an issue branch on your local system named **issue**```<NNN>_<short-heading>``` where NNN is the issue number from GitHub. This will be the main (short living) working area.
2. Write the necessary code.
    a. Make sure you test your modifications well. 
    b. Feel free to commit and push your issue branch often
3. Once the issue is fixed merge the **develop** branch back into your issue branch and resolve any conflicts.
4. Create a pull request on GitHub to merge the issue branch back into **develop**. Always include at least one reviewer who will then merge and delete the issue branch, and close the issue.

