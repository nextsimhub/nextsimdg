/*!
 * @file IDynamics.hpp
 *
 * @date 7 Sep 2023
 * @author Tim Spain <timothy.spain@nersc.no>
 */

#ifndef IDYNAMICS_HPP
#define IDYNAMICS_HPP

#include "include/ModelComponent.hpp"
#include "include/Time.hpp"
#include "include/gridNames.hpp"

namespace Nextsim {
class IDynamics : public ModelComponent {
public:
    IDynamics()
        : IDynamics(false)
    {
    }
    IDynamics(bool usesDamageIn)
        : uice(ModelArray::Type::H)
        , vice(ModelArray::Type::H)
        , damage(ModelArray::Type::H)
        , hice(getStore())
        , cice(getStore())
        , hsnow(getStore())
        , damage0(getStore())
        , uwind(getStore())
        , vwind(getStore())
        , uocean(getStore())
        , vocean(getStore())
        , m_usesDamage(usesDamageIn)
    {
        getStore().registerArray(Shared::DAMAGE, &damage, RW);
    }
    virtual ~IDynamics() = default;

    ModelState getState() const override
    {
        return { {
                     { uName, mask(uice) },
                     { vName, mask(vice) },
                 },
            {} };
    }
    ModelState getState(const OutputLevel&) const override { return getState(); }
    ModelState getStateRecursive(const OutputSpec& os) const override
    {
        // Ensure the base class implementation of getState() is called
        return os ? IDynamics::getState() : ModelState();
    }

    std::string getName() const override { return "IDynamics"; }
    virtual void setData(const ModelState::DataMap& ms) override
    {
        uice.resize();
        vice.resize();
        damage.resize();
        if (!m_usesDamage) {
            damage = 0.;
        }
    }

    virtual void update(const TimestepTime& tst) = 0;

    /*!
     * Advects a ModelArray field of any type
     *
     * @param field The ModelArray data to be advected
     *
     * @param fieldName A name for the field, to allow storage of additional information
     */
    ModelArray& advectField(ModelArray& field, const std::string& fieldName)
    {
        if (field.nDimensions() > 2) {
            // Get the total number of higher dimensions. ModelArray::Type::H is assumed to be the
            // size of a 2D field.
            size_t size2D = ModelArray::size(ModelArray::Type::H);
            size_t n2DFields = field.size() / size2D;
            HField levelData;
            levelData.resize();
            for (size_t k = 0; k < n2DFields; ++k) {
                // Copy the data of one 2D field to levelData
                levelData.setData(field, k * size2D, 0);
                // generate a consistent, unique name for this level
                std::string field2DName = fieldName + generateSuffix(k, n2DFields);
                // Advection occurs!
                advectHField(levelData, field2DName);
                // Copy the data of levelData back to the correct 2D field
                field.setData(levelData, 0, k * size2D);
            }
            return field;
        } else {
            return advectHField(field, fieldName);
        }
    }

    /*!
     * Returns whether the dynamics implementation uses the damage field.
     */
    virtual bool usesDamage() const { return m_usesDamage; }

protected:
    // Shared ice velocity arrays
    HField uice;
    HField vice;
    // Updated damage array
    HField damage;
    // References to the DG0 finite volume data arrays
    ModelArrayRef<Shared::H_ICE, RW> hice;
    ModelArrayRef<Shared::C_ICE, RW> cice;
    ModelArrayRef<Shared::H_SNOW, RW> hsnow;
    ModelArrayRef<Protected::DAMAGE, RO> damage0;
    // ModelArrayRef<ModelComponent::SharedArray::D, MARBackingStore, RW> damage;

    // References to the forcing velocity arrays
    ModelArrayRef<Protected::WIND_U> uwind;
    ModelArrayRef<Protected::WIND_V> vwind;
    ModelArrayRef<Protected::OCEAN_U> uocean;
    ModelArrayRef<Protected::OCEAN_V> vocean;

    // Does this implementation of the dynamics use damage?
    bool m_usesDamage;

    /*
     * Checks and returns if the provided data map is spherical
     */
    static bool checkSpherical(const ModelState::DataMap& ms)
    {
        // Decide between Cartesian (x & y) and spherical (longitude & latitude)
        if (ms.count(longitudeName) > 0 && ms.count(latitudeName) > 0) {
            return true;
        } else if (ms.count(xName) > 0 && ms.count(yName) > 0) {
            return false;
        } else {
            // Throw a runtime_error exception which can either be handled or not
            throw std::runtime_error("Input data must contain either Cartesian (" + xName + ", "
                + yName + ") or spherical (" + longitudeName + ", " + latitudeName
                + ") coordinates.");
        }
    }

    /*!
     * Advects an HField, likely using the implementation's dynamics kernel
     */
    virtual ModelArray& advectHField(ModelArray& field, const std::string& fieldName) = 0;

private:
    // Generates a consistent suffix for multilevel fields
    std::string generateSuffix(size_t i, size_t n) {
        if (n < 2) return "";
        return "_" + std::to_string(i) + "OF" + std::to_string(n);
    }

};
}

#endif /* IDYNAMICS_HPP */
