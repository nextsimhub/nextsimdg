/*!
 * @file    Xios.hpp
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    02 Jan 2025
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

#include "Configured.hpp"
#include "Logged.hpp"
#include "ModelArray.hpp"
#include "Time.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>
#include <include/xios_c_interface.hpp>
#include <mpi.h>

namespace Nextsim {

void enableXios();

class Xios : public Configured<Xios> {
private:
    static Xios* instancePtr;

    // Private constructor
    Xios(const std::string dt = "P0-0T01:00:00", const std::string contextid = "nextSIM-DG",
        const std::string starttime = "1970-01-01T00:00:00Z",
        const std::string calendartype = "Gregorian");

public:
    Xios(const Xios& xiosHandler) = delete;
    ~Xios();

    /*
     * Define Xios handler Singleton
     *
     * NOTE: The arguments will only be used the first time this is called.
     */
    inline static Xios* getInstance(const std::string dt = "P0-0T01:00:00",
        const std::string contextId = "nextSIM-DG",
        const std::string startTime = "1970-01-01T00:00:00Z",
        const std::string calendarType = "Gregorian")
    {
        if (isNull()) {
            instancePtr = new Xios(dt, contextId, startTime, calendarType);
        }
        return instancePtr;
    };

    /* Utility for checking if the singleton has been created */
    inline static bool isNull() { return (instancePtr == nullptr); }

    /* Initialization and finalization */
    void close_context_definition();
    void context_finalize();
    void finalize();
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
    void createDomain(const std::string domainId);
    void setDomainType(const std::string domainId, const std::string domainType);
    void setDomainGlobalXSize(const std::string domainId, const size_t size);
    void setDomainGlobalYSize(const std::string domainId, const size_t size);
    void setDomainLocalXSize(const std::string domainId, const size_t size);
    void setDomainLocalYSize(const std::string domainId, const size_t size);
    void setDomainLocalXStart(const std::string domainId, const size_t start);
    void setDomainLocalYStart(const std::string domainId, const size_t start);
    void setDomainLocalXValues(const std::string domainId, std::vector<double> values);
    void setDomainLocalYValues(const std::string domainId, std::vector<double> values);
    std::string getDomainType(const std::string domainId);
    size_t getDomainGlobalXSize(const std::string domainId);
    size_t getDomainGlobalYSize(const std::string domainId);
    size_t getDomainLocalXSize(const std::string domainId);
    size_t getDomainLocalYSize(const std::string domainId);
    size_t getDomainLocalXStart(const std::string domainId);
    size_t getDomainLocalYStart(const std::string domainId);
    std::vector<double> getDomainLocalXValues(const std::string domainId);
    std::vector<double> getDomainLocalYValues(const std::string domainId);

    /* Grid */
    void createGrid(const std::string gridId);
    void setGridName(const std::string gridId, const std::string name);
    std::string getGridName(const std::string gridId);
    void gridAddAxis(std::string axisId, const std::string domainId);
    void gridAddDomain(const std::string gridId, const std::string domainId);
    std::vector<std::string> gridGetAxisIds(const std::string gridId);
    std::vector<std::string> gridGetDomainIds(const std::string gridId);

    /* Field */
    void createField(const std::string fieldId);
    void setFieldName(const std::string fieldId, const std::string name);
    void setFieldOperation(const std::string fieldId, const std::string operation);
    void setFieldGridRef(const std::string fieldId, const std::string gridRef);
    void setFieldReadAccess(const std::string fieldId, const bool readAccess);
    void setFieldFreqOffset(const std::string fieldId, const Duration freqOffset);
    std::string getFieldName(const std::string fieldId);
    std::string getFieldOperation(const std::string fieldId);
    std::string getFieldGridRef(const std::string fieldId);
    bool getFieldReadAccess(const std::string fieldId);
    Duration getFieldFreqOffset(const std::string fieldId);

    /* File */
    void createFile(const std::string fileId);
    void setFileName(const std::string fileId, const std::string fileName);
    void setFileType(const std::string fileId, const std::string fileType);
    void setFileOutputFreq(const std::string fileId, const Duration outputFreq);
    void setFileSplitFreq(const std::string fileId, const Duration splitFreq);
    void setFileMode(const std::string fileId, const std::string mode);
    void setFileParAccess(const std::string fileId, const std::string parAccess);
    std::string getFileName(const std::string fileId);
    std::string getFileType(const std::string fileId);
    Duration getFileOutputFreq(const std::string fileId);
    Duration getFileSplitFreq(const std::string fileId);
    std::string getFileMode(const std::string fileId);
    std::string getFileParAccess(const std::string fileId);
    void fileAddField(const std::string fileId, const std::string fieldId);
    std::vector<std::string> fileGetFieldIds(const std::string fileId);

    /* I/O */
    void write(const std::string fieldId, ModelArray& modelarray);
    void read(const std::string fieldId, ModelArray& modelarray);

    enum {
        ENABLED_KEY,
    };

    /* Length of C-strings passed to XIOS */
    int cStrLen { 20 };

protected:
    bool isConfigured;

private:
    bool isEnabled;

    std::string clientId;
    std::string calendarType;
    std::string contextId;
    Duration timestep;
    TimePoint startTime;
    MPI_Comm clientComm;
    MPI_Fint clientComm_F;
    MPI_Fint nullComm_F;
    int mpi_rank { 0 };
    int mpi_size { 0 };

    xios::CCalendarWrapper* clientCalendar;
    std::string convertXiosDatetimeToString(const cxios_date datetime, const bool isoFormat = true);
    cxios_date convertStringToXiosDatetime(const std::string datetime, const bool isoFormat = true);
    std::string convertCStrToCppStr(const char* cStr, int cStrLen);
    Duration convertDurationFromXios(const cxios_duration duration);
    cxios_duration convertDurationToXios(const Duration duration);

    xios::CAxisGroup* getAxisGroup();
    xios::CDomainGroup* getDomainGroup();
    xios::CFieldGroup* getFieldGroup();
    xios::CGridGroup* getGridGroup();
    xios::CFileGroup* getFileGroup();

    xios::CAxis* getAxis(const std::string axisId);
    xios::CDomain* getDomain(const std::string domainId);
    xios::CField* getField(const std::string fieldId);
    xios::CGrid* getGrid(const std::string gridId);
    xios::CFile* getFile(const std::string fileId);
};

}

#endif // USE_XIOS
#endif
