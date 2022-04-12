.. Copyright (c) 2021, Nansen Environmental and Remote Sensing Center


neXtSIM_DG : the next generation sea-ice model with a Discontinuous Galerkin method

Introduction
------------

**neXtSIM is a Lagrangian dynamical-thermodynamical sea-ice model that faithfully represents large-scale sea ice thickness, concentration and drift patterns** [1]_ [2]_ **, statistical properties of pack ice deformation and drift** [3]_ [2]_ **, diffusion and dispersion** [4]_ [5]_ **; and additionally exhibits a strong capacity for short-term sea ice forecasting at Pan-Arctic scales** [6]_**. neXtSIM was developed to explore the role of sea-ice mechanics and dynamics in determining the large-scale behaviour of Arctic sea ice. As such, it has served as a testing and development ground for the brittle sea-ice rheologies developed over the past decade** [7]_ [8]_ [9]_ [2]_**.** The initial development of the model was focused on understanding and representing as well as possible the impact of the new physics developments, which was the main rationale for the Lagrangian advection scheme used. As the new physics developments matured it became evident that these could impact not only the ice, but also the ocean and atmosphere below and above it [10]_.

In order to better understand the role of sea-ice mechanics and dynamics in the coupled atmosphere-ocean-ice system we are developing a new version of neXtSIM, neXtSIM_DG. neXtSIM_DG is the Eulerian version of the neXtSIM model, using the discontinuous Galerkin method (DG) for advection, instead of the Lagrangian scheme used in previous neXtSIM versions. The DG method that has proven extremely effective and efficient in retaining high gradients in fluid dynamics [11]_ [12]_ , while the Lagrangian neXtSIM proved challenging to couple with other models and difficult to optimise for highly parallel environments.

The core idea behind neXtSIM_DG remains the same as that behind the original neXtSIM: to explore the role of sea-ice mechanics and dynamics in determining the large-scale behaviour of Arctic sea ice. With the new version we aim to extend the use of the model to fully coupled and global simulations, resulting in a model that can be used from regional short term forecasting applications to global climate simulations. To achieve this, the model infrastructure will be modular and flexible, so that it can easily be coupled to other models and systems and so that new parameterisations and physical routines can easily be added. The dynamical core will use the DG method to ensure the best possible preservation of high gradients in an Eulerian framework, while still maintaining good computational performance.

.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/ZLhtYUB_Xfo" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
  
.. raw:: html

    <div style="position: relative; padding-bottom: 56.25%; height: 0; overflow: hidden; max-width: 100%; height: auto;">
        <iframe src="https://www.youtube.com/embed/gQq7YYOlsHk" frameborder="0" allowfullscreen style="position: absolute; top: 0; left: 0; width: 100%; height: 100%;"></iframe>
    </div>
  

References
----------

.. [1] P. Rampal, S. Bouillon, E. Olason, and M. Morlighem. neXtSIM: a new Lagrangian sea ice model.  ́ The Cryosphere, 10:1055–1073, 2016.
.. [2] Olason, E., Boutin, G., Korosov, A., Rampal, P., Williams, T., Kimmritz, M., Dansereau V., Samaké, A., A new brittle rheology and numerical framework for large-scale sea-ice models. Earth and Space Science Open Archive (2021) doi:10.1002/essoar.10507977.2.
.. [3] P. Rampal, V. Dansereau, E. Olason, S. Bouillon, T. Williams, A. Korosov, and A. Samake. On the multi-fractal scaling properties of sea ice deformation. The Cryosphere, 13(9):2457–2474, 2019.
.. [4] M. Rabatel, P. Rampal, A. Carrassi, L. Bertino, and C. K. R. T. Jones. Impact of rheology on probabilistic forecasts of sea ice trajectories: application for search and rescue operations in the arctic. The Cryosphere, 12(3):935–953, 2018.
.. [5] P. Rampal, S. Bouillon, J. Bergh, and E. Olason. Arctic sea-ice diffusion from observed and simulated Lagrangian trajectories. The Cryosphere, 10(4):1513–1527, 2016.
.. [6] T. Williams, A. Korosov, P. Rampal, and E. Olason. Presentation and evaluation of the arctic sea ice forecasting system nextsim-f. The Cryosphere Discussions, 2019:1–31, 2019.
.. [7] Girard, L., Bouillon, S., Weiss, J., Amitrano, D., Fichefet, T., & Legat, V. (2011). A new modeling framework for sea-ice mechanics based on elasto-brittle rheology. Annals of Glaciology, 52, 123–132, doi:10.3189/172756411795931499
.. [8] Bouillon, S., & Rampal, P. (2015). Presentation of the dynamical core of neXtSIM, a new sea ice model. Ocean Modelling, 91, 23–37, doi:10.1016/j.ocemod.2015.04.005
.. [9] Dansereau, V., Weiss, J., Saramito, P., & Lattes, P. (2016). A Maxwell elasto-brittle rheology for sea ice modelling. Cryosphere, 10, 1339–1359, doi:10.5194/tc-10-1339-2016
.. [10] Ólason, E., Rampal, P., & Dansereau, V. (2021). On the statistical properties of sea-ice lead fraction and heat fluxes in the Arctic. The Cryosphere, 15, 1053–1064, doi:10.5194/tc-15-1053-2021
.. [11] E. Di Pietro and A. Ern. Mathematical Aspects of Discontinuous Galerkin Methods, volume 69 of Mathematiques et Applications. Springer-Verlag Berlin Heidelberg, 2012.
.. [12] P. Saramito. Efficient C++ finite element computing with Rheolef. CNRS-CCSD ed., 2013. http://cel.archives-ouvertes.fr/cel-00573970.




Licensing
---------

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.

You may obtain a copy of the License at
    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

See the License for the specific language governing permissions and
limitations under the License.


.. toctree::
   :caption: INSTALLATION
   :maxdepth: 1

   installation
   changelog

.. toctree::
   :caption: USAGE
   :maxdepth: 2

   getting_started


.. toctree::
   :caption: API Documentation
   :maxdepth: 2

   api/library_root
   api/index


All `Doxygen documentation`_

.. _Doxygen documentation : https://nextsimdg.github.io/nextsimdg/index.html
