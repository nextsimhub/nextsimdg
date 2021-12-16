        if (module == "ITest") {
            if (impl == "Impl1") {
                p_ITest = &i_Impl1;
                pf_ITest = &newImpl1;
            } else if (impl == "Impl2") {
                p_ITest = &i_Impl2;
                pf_ITest = &newImpl2;
            } else {
                throwup(module, impl);
            }
        } else { }
