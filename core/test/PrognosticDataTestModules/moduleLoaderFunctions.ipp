static Nextsim::IFreezingPoint* p_IFreezingPoint;
template <> Nextsim::IFreezingPoint& ModuleLoader::getImplementation() { return *p_IFreezingPoint; }
std::unique_ptr<Nextsim::IFreezingPoint> (*pf_IFreezingPoint)();
template <> std::unique_ptr<Nextsim::IFreezingPoint> ModuleLoader::getInstance() const
{
    return (*pf_IFreezingPoint)();
}
static Nextsim::LinearFreezing i_LinearFreezing;
std::unique_ptr<Nextsim::IFreezingPoint> newLinearFreezing()
{
    return std::unique_ptr<Nextsim::LinearFreezing>(new Nextsim::LinearFreezing);
}
static Nextsim::UnescoFreezing i_UnescoFreezing;
std::unique_ptr<Nextsim::IFreezingPoint> newUnescoFreezing()
{
    return std::unique_ptr<Nextsim::UnescoFreezing>(new Nextsim::UnescoFreezing);
}
