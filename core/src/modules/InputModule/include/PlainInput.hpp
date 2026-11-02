/*!
 * @author  Einar Olason <einar.olason@nersc.no>
 */

#include "include/IInput.hpp"
#include "include/ParaGridIO.hpp"

#ifndef PLAININPUT_HPP
#define PLAININPUT_HPP

namespace Nextsim {

class PlainInput : public IInput {
public:
    PlainInput() = default;

    /*!
     *@brief Initialise object parameters and calculate interpolation weights. To be called by
     * parent's setData function.
     *
     * @param time Initial time
     * @param pathSpecIn Path of for the input file.
     * @param ncLonNameIn Not used
     * @param ncLatNameIn Not used
     * @param ncTimeNameIn Name of the time variable in the netCDF file
     * @param forcingsIn Set of forcing variables to read
     * @param vectorsIn Not used
     * @param modelLonsIn Not used
     * @param modelLatsIn Not used
     */
    void initialise(const TimePoint& time, const std::string& pathSpecIn,
        const std::string& ncLonNameIn, const std::string& ncLatNameIn,
        const std::string& ncTimeNameIn, const std::set<std::string>& forcingsIn,
        const std::set<std::pair<std::string, std::string>>& vectorsIn,
        const ModelArray& modelLonsIn, const ModelArray& modelLatsIn) override;

    /*!
     * @brief Update the forcing data for the current time step.
     *
     * @param time Current model time
     */
    void update(const TimePoint& time) override;

    /*!
     * @brief Get the forcing data for a specific field interpolated in time to currentTime.
     *
     * @param fieldName Name of the field to retrieve
     * @return Forcing data for the specified field
     */
    [[nodiscard]] ModelArray getField(const std::string& fieldName) const override;

    /*!
     * @brief Check if a field can be read in from the forcing data set.
     *
     * @param fieldName Name of the field to check
     * @return True if the field is valid, false otherwise
     */
    [[nodiscard]] bool isValid(const std::string& fieldName) const override;

    // ModelComponent overrides
    std::string getName() const override { return "ParaGridInputs"; };
    void setData(const ModelState::DataMap&) override { };

private:
    ModelState forcingState;
    std::set<std::string> forcings;
    std::string filePath;
    double deltaTime;
    std::map<std::string, bool> warn;
    TimePoint initTime;
};

}

#endif // PLAININPUT_HPP
