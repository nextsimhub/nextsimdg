/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @brief   XIOS interface header
 * @details
 *
 * Header file for XIOS interface
 *
 */
#ifndef SRC_INCLUDE_XIOS_HPP
#define SRC_INCLUDE_XIOS_HPP

#include "date.hpp"
#if USE_XIOS

#include "include/Configured.hpp"
#include "include/Iterator.hpp"
#include "include/Logged.hpp"
#include "include/ModelArray.hpp"
#include "include/ModelConfig.hpp"
#include "include/Time.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>
#include <include/xios_c_interface.hpp>
#include <mpi.h>

namespace Nextsim {

// Forward declarations to avoid circular dependencies
class ModelMetadata;
class ParaGridIO;

void enableXios();

class Xios : public Configured<Xios> {
private:
    //! Private constructor
    Xios(const std::string contextId = "nextSIM-DG", const std::string calendarType = "Gregorian");

    //! Performs some one-time initialization for the class. Returns true.
    static bool doOnce();

public:
    ~Xios();

    //! Prevent copying
    Xios(const Xios&) = delete;

    /*
     * Define Xios handler Singleton
     *
     * NOTE: The arguments will only be used the first time this is called.
     */
    inline static Xios& getInstance(
        const std::string contextId = "nextSIM-DG", const std::string calendarType = "Gregorian")
    {
        static Xios instance = Xios(contextId, calendarType);
        return instance;
    };

    /* Help config */
    ConfigMap getConfig() const;
    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    /* Initialization and finalization */
    void close_context_definition();
    void context_finalize();
    static void finalize();
    bool isInitialized();

    /* Configuration */
    void configure() override;
    void configureServer();

    /* MPI decomposition */
    int getClientMPISize();
    int getClientMPIRank();

    /* Calendar, date and duration */
    void setCalendarType(const std::string type);
    void setCalendarOrigin(const TimePoint origin);
    void setCalendarStart(const TimePoint start);
    void setCalendarTimestep(const Duration timestep);
    void setCalendarStep(const int stepNumber);
    void incrementCalendar();
    std::string getCalendarType();
    TimePoint getCalendarOrigin();
    TimePoint getCalendarStart();
    Duration getCalendarTimestep();
    int getCalendarStep();
    std::string getCurrentDate(const bool isoFormat = true);

    /* Axis */
    void createAxis(const std::string axisId);
    void setAxisSize(const std::string axisId, const size_t size);
    void setAxisValues(const std::string axisId, std::vector<double> values);
    size_t getAxisSize(const std::string axisId);
    std::vector<double> getAxisValues(const std::string axisId);

    /* Domain */
    void affixModelMetadata(ModelMetadata& metadata);

    /* Grid */
    void createGrid(const std::string gridId);
    void gridAddAxis(std::string axisId, const std::string gridId);
    std::vector<std::string> getGridAxisIds(const std::string gridId);

    /* Field */
    void createField(const std::string fieldId);
    void setFieldOperation(const std::string fieldId, const std::string operation);
    void setFieldGridRef(const std::string fieldId, const std::string gridRef);
    void setFieldFreqOffset(const std::string fieldId, const Duration freqOffset);
    std::string getFieldOperation(const std::string fieldId);
    std::string getFieldGridRef(const std::string fieldId);
    bool getFieldReadAccess(const std::string fieldId);
    Duration getFieldFreqOffset(const std::string fieldId);

    /* File */
    void createFile(const std::string fileId);
    void setFileType(const std::string fileId, const std::string fileType);
    void setFileOutputFreq(const std::string fileId, const Duration outputFreq);
    void setFileSplitFreq(const std::string fileId, const Duration splitFreq);
    void setFileParAccess(const std::string fileId, const std::string parAccess);
    std::string getFileType(const std::string fileId);
    Duration getFileOutputFreq(const std::string fileId);
    Duration getFileSplitFreq(const std::string fileId);
    std::string getFileMode(const std::string fileId);
    std::string getFileParAccess(const std::string fileId);
    void fileAddField(const std::string fileId, const std::string fieldId);
    std::vector<std::string> fileGetFieldIds(const std::string fileId);

    /* I/O */
    void read(const std::string fieldId, ModelArray& modelarray);

    enum {
        ENABLED_KEY,
        STARTTIME_KEY,
        STOPTIME_KEY,
        TIME_STEP_KEY,
        READ_MODE_KEY,
        OUTPUT_RESTARTPERIOD_KEY,
        OUTPUT_RESTARTFILE_KEY,
        OUTPUT_FIELD_NAMES_KEY,
        INPUT_RESTARTPERIOD_KEY,
        INPUT_RESTARTFILE_KEY,
        INPUT_FIELD_NAMES_KEY,
    };

    /* Length of C-strings passed to XIOS */
    int cStrLen { 20 };

protected:
    bool isConfigured;

private:
    inline static bool isEnabled = false;
    std::string clientId;
    std::string contextId;
    MPI_Comm clientComm;
    MPI_Fint clientComm_F;
    MPI_Fint nullComm_F;
    int mpi_rank { 0 };
    int mpi_size { 0 };

    /* Calendar, date and duration */
    std::string calendarType;
    Duration timestep;
    TimePoint startTime;
    TimePoint stopTime;
    xios::CCalendarWrapper* clientCalendar;
    std::string convertXiosDatetimeToString(const cxios_date datetime, const bool isoFormat = true);
    cxios_date convertStringToXiosDatetime(const std::string datetime, const bool isoFormat = true);
    std::string convertCStrToCppStr(const char* cStr, int cStrLen);
    Duration convertDurationFromXios(const cxios_duration duration);
    cxios_duration convertDurationToXios(const Duration duration);

    /* Axis */
    xios::CAxisGroup* getAxisGroup();
    xios::CAxis* getAxis(const std::string axisId);

    /* Domain */
    const std::string domainId = "xy_domain";
    xios::CDomainGroup* getDomainGroup();
    xios::CDomain* getDomain();

    /* Field */
    xios::CFieldGroup* getFieldGroup();
    xios::CField* getField(const std::string fieldId);
    void setFieldReadAccess(const std::string fieldId, const bool readAccess);
    std::set<std::string> configGetFieldNames(const bool reading);
    bool configCheckField(const std::string fieldId, const bool reading);

    /* Grid */
    xios::CGridGroup* getGridGroup();
    xios::CGrid* getGrid(const std::string gridId);

    /* File */
    xios::CFileGroup* getFileGroup();
    xios::CFile* getFile(const std::string fileId);
    void setFileMode(const std::string fileId, const std::string mode);
    std::string inputFilename;
    std::string inputFileId;
    std::string outputFilename;
    std::string outputFileId;

    /* I/O */
    void write(const std::string fieldId, ModelArray& modelarray);

    /* Declare any classes that need to access private members */
    friend ParaGridIO;
};

}

#endif // USE_XIOS
#endif
