# Changelog

## [Unreleased](https://github.com/nextsimdg/nextsimdg/tree/HEAD)

[Full Changelog](https://github.com/nextsimdg/nextsimdg/compare/v0.2.0...HEAD)

**Fixed bugs:**

- Initial ice depends on a call to aoState.update\(\) [\#140](https://github.com/nextsimdg/nextsimdg/issues/140)

**Closed issues:**

- Integrating MEB to parametric mesh [\#128](https://github.com/nextsimdg/nextsimdg/issues/128)
- Implement logging [\#67](https://github.com/nextsimdg/nextsimdg/issues/67)

## [v0.2.0](https://github.com/nextsimdg/nextsimdg/tree/v0.2.0) (2022-09-05)

**Description:**
Last version including both meshes before merge


[Full Changelog](https://github.com/nextsimdg/nextsimdg/compare/v0.1.0...v0.2.0)

**Implemented enhancements:**

- Passing time values to parameterizations [\#133](https://github.com/nextsimdg/nextsimdg/issues/133)
- More complex output: ConfiguredOutput [\#124](https://github.com/nextsimdg/nextsimdg/issues/124)
- Land-sea mask [\#107](https://github.com/nextsimdg/nextsimdg/issues/107)
- Implement simple output [\#100](https://github.com/nextsimdg/nextsimdg/issues/100)
- Targetted C++ standard could be C++ 2017 [\#85](https://github.com/nextsimdg/nextsimdg/issues/85)
- Self-contained modules [\#82](https://github.com/nextsimdg/nextsimdg/issues/82)
- Modules should be separable [\#59](https://github.com/nextsimdg/nextsimdg/issues/59)
- Restart file structure [\#58](https://github.com/nextsimdg/nextsimdg/issues/58)
- Setup for clang-format [\#9](https://github.com/nextsimdg/nextsimdg/issues/9)

**Fixed bugs:**

- IConcentrationModel is no longer used [\#122](https://github.com/nextsimdg/nextsimdg/issues/122)

**Closed issues:**

- How does the Module system work with undefined member functions? [\#117](https://github.com/nextsimdg/nextsimdg/issues/117)
- Configuration dataflow [\#113](https://github.com/nextsimdg/nextsimdg/issues/113)
- ModelArrayRef should provide access to Eigen::Matrix data [\#111](https://github.com/nextsimdg/nextsimdg/issues/111)
- Time parsing [\#102](https://github.com/nextsimdg/nextsimdg/issues/102)
- Duplicate code in CI jobs [\#88](https://github.com/nextsimdg/nextsimdg/issues/88)
- Crash with no command line options [\#78](https://github.com/nextsimdg/nextsimdg/issues/78)
- A more flexible rectangular grid [\#72](https://github.com/nextsimdg/nextsimdg/issues/72)
- Use CMake imported targets when available [\#71](https://github.com/nextsimdg/nextsimdg/issues/71)
- Prognostic ice data types [\#25](https://github.com/nextsimdg/nextsimdg/issues/25)
- Github workflow for CI [\#15](https://github.com/nextsimdg/nextsimdg/issues/15)

## [v0.1.0](https://github.com/nextsimdg/nextsimdg/tree/v0.1.0) (2022-02-11)

**Description:**
Merge pull request #74 from nextsimdg/develop

Minimum Viable Model from develop.

[Full Changelog](https://github.com/nextsimdg/nextsimdg/compare/cb4f6ad730c95330c552fe27498e6d0a9fba40f5...v0.1.0)

**Implemented enhancements:**

- Initialize PrognosticData using a IPrognosticUpdater instance [\#65](https://github.com/nextsimdg/nextsimdg/issues/65)
- The number of ice levels should be defined by the physics [\#64](https://github.com/nextsimdg/nextsimdg/issues/64)
- Model logic for controlling \(physics\) timesteps [\#57](https://github.com/nextsimdg/nextsimdg/issues/57)
- Model logic for reading and writing restart files [\#56](https://github.com/nextsimdg/nextsimdg/issues/56)
- Create grid/mesh structure [\#55](https://github.com/nextsimdg/nextsimdg/issues/55)
- Restart file class\(es\) [\#54](https://github.com/nextsimdg/nextsimdg/issues/54)
- Minimum viable model [\#52](https://github.com/nextsimdg/nextsimdg/issues/52)
- Make arguments to module functions const [\#44](https://github.com/nextsimdg/nextsimdg/issues/44)
- Implement an actual Timer class [\#22](https://github.com/nextsimdg/nextsimdg/issues/22)
- Configuration [\#14](https://github.com/nextsimdg/nextsimdg/issues/14)
- Module loading [\#10](https://github.com/nextsimdg/nextsimdg/issues/10)
- Constants [\#6](https://github.com/nextsimdg/nextsimdg/issues/6)
- Implement ElementData [\#4](https://github.com/nextsimdg/nextsimdg/issues/4)
- Minimum viable model [\#53](https://github.com/nextsimdg/nextsimdg/pull/53) ([timspainNERSC](https://github.com/timspainNERSC))
- Merge the ElementData feature branch [\#5](https://github.com/nextsimdg/nextsimdg/pull/5) ([timspainNERSC](https://github.com/timspainNERSC))

**Fixed bugs:**

- NextsimPhysics should be a module [\#60](https://github.com/nextsimdg/nextsimdg/issues/60)
- ConfiguredModule should only configure modules that are defined. [\#48](https://github.com/nextsimdg/nextsimdg/issues/48)
- Modules configuration section [\#43](https://github.com/nextsimdg/nextsimdg/issues/43)

**Closed issues:**

- ConfiguredModule argument checking [\#50](https://github.com/nextsimdg/nextsimdg/issues/50)
- Inconsistent test results between mac and ubuntu [\#42](https://github.com/nextsimdg/nextsimdg/issues/42)
- Fix formating for clang format [\#20](https://github.com/nextsimdg/nextsimdg/issues/20)
- Change indentation from tabs to spaces [\#1](https://github.com/nextsimdg/nextsimdg/issues/1)

**Merged pull requests:**

- Docs [\#77](https://github.com/nextsimdg/nextsimdg/pull/77) ([timspainNERSC](https://github.com/timspainNERSC))
- merge develop into docs [\#76](https://github.com/nextsimdg/nextsimdg/pull/76) ([auraoupa](https://github.com/auraoupa))
- Develop [\#74](https://github.com/nextsimdg/nextsimdg/pull/74) ([timspainNERSC](https://github.com/timspainNERSC))
- Issue64 ice levels [\#66](https://github.com/nextsimdg/nextsimdg/pull/66) ([timspainNERSC](https://github.com/timspainNERSC))
- Issue60 Physics Module [\#61](https://github.com/nextsimdg/nextsimdg/pull/61) ([timspainNERSC](https://github.com/timspainNERSC))
- Issue48 configmod [\#51](https://github.com/nextsimdg/nextsimdg/pull/51) ([timspainNERSC](https://github.com/timspainNERSC))
- Doxygen comments for all\(?\) of the headers [\#49](https://github.com/nextsimdg/nextsimdg/pull/49) ([timspainNERSC](https://github.com/timspainNERSC))
- Put module configuration in a \[Modules\] configuration section. [\#47](https://github.com/nextsimdg/nextsimdg/pull/47) ([timspainNERSC](https://github.com/timspainNERSC))
- WIP: Physics subdirectory [\#45](https://github.com/nextsimdg/nextsimdg/pull/45) ([timspainNERSC](https://github.com/timspainNERSC))
- Develop [\#41](https://github.com/nextsimdg/nextsimdg/pull/41) ([auraoupa](https://github.com/auraoupa))
- Merge develop into the timer branch [\#40](https://github.com/nextsimdg/nextsimdg/pull/40) ([timspainNERSC](https://github.com/timspainNERSC))
- Update README.md [\#38](https://github.com/nextsimdg/nextsimdg/pull/38) ([auraoupa](https://github.com/auraoupa))
- Issue15 ci workflow [\#37](https://github.com/nextsimdg/nextsimdg/pull/37) ([auraoupa](https://github.com/auraoupa))
- Merge from develop [\#36](https://github.com/nextsimdg/nextsimdg/pull/36) ([timspainNERSC](https://github.com/timspainNERSC))
- Merge pull request \#34 from nextsimdg/issue14\_config [\#35](https://github.com/nextsimdg/nextsimdg/pull/35) ([timspainNERSC](https://github.com/timspainNERSC))
- Parse for command line for config files, and load them to Configurator. [\#34](https://github.com/nextsimdg/nextsimdg/pull/34) ([timspainNERSC](https://github.com/timspainNERSC))
- Issue22 timer [\#33](https://github.com/nextsimdg/nextsimdg/pull/33) ([timspainNERSC](https://github.com/timspainNERSC))
- Module loading into ElementData [\#32](https://github.com/nextsimdg/nextsimdg/pull/32) ([timspainNERSC](https://github.com/timspainNERSC))
- Issue15 ci workflow [\#31](https://github.com/nextsimdg/nextsimdg/pull/31) ([auraoupa](https://github.com/auraoupa))
- Configuration infrastructure into Module Loading [\#29](https://github.com/nextsimdg/nextsimdg/pull/29) ([timspainNERSC](https://github.com/timspainNERSC))
- Issue22: Timer infrastructure [\#28](https://github.com/nextsimdg/nextsimdg/pull/28) ([timspainNERSC](https://github.com/timspainNERSC))
- Develop [\#27](https://github.com/nextsimdg/nextsimdg/pull/27) ([auraoupa](https://github.com/auraoupa))
- Merger develop into timer issue branch [\#26](https://github.com/nextsimdg/nextsimdg/pull/26) ([timspainNERSC](https://github.com/timspainNERSC))
- Get develop up to date with main [\#24](https://github.com/nextsimdg/nextsimdg/pull/24) ([einola](https://github.com/einola))
- A format file for clang\_format [\#23](https://github.com/nextsimdg/nextsimdg/pull/23) ([einola](https://github.com/einola))
- Issue20 clang formatting [\#21](https://github.com/nextsimdg/nextsimdg/pull/21) ([timspainNERSC](https://github.com/timspainNERSC))
- Merge config changes into the elementdata feature branch [\#19](https://github.com/nextsimdg/nextsimdg/pull/19) ([timspainNERSC](https://github.com/timspainNERSC))
- clang format in the CI tests [\#18](https://github.com/nextsimdg/nextsimdg/pull/18) ([auraoupa](https://github.com/auraoupa))
- rebase from main [\#17](https://github.com/nextsimdg/nextsimdg/pull/17) ([auraoupa](https://github.com/auraoupa))
- Add ModuleLoader to the feature branch source tree with no current modules. [\#13](https://github.com/nextsimdg/nextsimdg/pull/13) ([timspainNERSC](https://github.com/timspainNERSC))
- Merge from develop [\#12](https://github.com/nextsimdg/nextsimdg/pull/12) ([timspainNERSC](https://github.com/timspainNERSC))
- Docs [\#8](https://github.com/nextsimdg/nextsimdg/pull/8) ([auraoupa](https://github.com/auraoupa))
- Create the namespaces, constants and conversion functions. [\#7](https://github.com/nextsimdg/nextsimdg/pull/7) ([timspainNERSC](https://github.com/timspainNERSC))
- Merge main into develop [\#3](https://github.com/nextsimdg/nextsimdg/pull/3) ([timspainNERSC](https://github.com/timspainNERSC))
- Spacify [\#2](https://github.com/nextsimdg/nextsimdg/pull/2) ([timspainNERSC](https://github.com/timspainNERSC))



\* *This Changelog was automatically generated by [github_changelog_generator](https://github.com/github-changelog-generator/github-changelog-generator)*
