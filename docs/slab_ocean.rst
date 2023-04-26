.. Copyright (c) 2023, Nansen Environmental and Remote Sensing Center

The neXtSIM_DG Slab Ocean Model
===============================

Introduction
------------

The slab ocean of neXtSIM_DG exists to ensure consistency between the
conditions the ice is experiencing and the ocean forcing applied to the model.

The slab ocean was first implemented in the Lagrangian version of neXtSIM
`(Rampal et al., 2016)`_ where the implementation benefitted from the monolithic
nature of the model core, allowing values to be calculated in whatever order
they needed to be. The modular design of neXtSIM_DG necessitated a close look
at the flow of data through the parts of the model that interact with the slab
ocean. These, along with the overall slab model, are described in this document.

.. (Rampal et al., 2016) Rampal, P., Bouillon, S., Ólason, E., and
   Morlighem, M., neXtSIM: a new Lagrangian sea ice model, *Cryosphere*,
   **10** 1055—1073 (2016)

Definitions
-----------

The model timestep begins at time *t₀* and integrates the model for a
period *δt* to the end of the timestep at *t₁*. The previous iteration of the
model, along with the forcing or coupling, provides the values of several
variables at *t* = *t₀*. The dynamics component of the model then advects these
values to time *t₁*, whereupon the physics part of the model attempts to update
the prognostic data fields to new values consistent with the motion due to the
dynamics and the effect of the passage of time on the model state.

.. The slab ocean functions by having two timescales over which the corrections apply, timeT and timeS