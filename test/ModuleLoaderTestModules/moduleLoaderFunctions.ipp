static ITest* p_ITest;
template<>
ITest& ModuleLoader::getImplementation<ITest>()
{
    return *p_ITest;
}
std::unique_ptr<ITest> (*pf_ITest)();
template<>
std::unique_ptr<ITest> ModuleLoader::getInstance<ITest>() const
{
    return (*pf_ITest)();
}
static Impl1 i_Impl1;
std::unique_ptr<ITest> newImpl1()
{
    return std::unique_ptr<ITest>(new Impl1);
}
static Impl2 i_Impl2;
std::unique_ptr<ITest> newImpl2()
{
    return std::unique_ptr<ITest>(new Impl2);
}
