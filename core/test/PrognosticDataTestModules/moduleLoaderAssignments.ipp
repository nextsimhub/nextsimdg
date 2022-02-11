        if (module == "Nextsim::IFreezingPoint") {
            if (impl == "Nextsim::LinearFreezing") {
                p_IFreezingPoint = &i_LinearFreezing;
                pf_IFreezingPoint = &newLinearFreezing;
            } else if (impl == "Nextsim::UnescoFreezing") {
                p_IFreezingPoint = &i_UnescoFreezing;
                pf_IFreezingPoint = &newUnescoFreezing;
            } else {
                throwup(module, impl);
            }

        } else { }
