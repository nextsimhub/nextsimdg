/*!
 * @file    Xios.hpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    26 July 2024
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

class Xios : public Configured<Xios> {
public:
    Xios();
    ~Xios();

    /* Initialization and finalization */
    void close_context_definition();
    void context_finalize();
    void finalize();
    bool isInitialized();

    /* Configuration */
    void configure() override;
    void configureServer(const std::string calendarType = "Gregorian");

    /* MPI decomposition */
    int getClientMPISize();
    int getClientMPIRank();

    /* Calendar, date and duration */
    void setCalendarType(const std::string type);
    void setCalendarOrigin(const TimePoint origin);
    void setCalendarStart(const TimePoint start);
    void setCalendarTimestep(const Duration timestep);
    std::string getCalendarType();
    TimePoint getCalendarOrigin();
    TimePoint getCalendarStart();
    Duration getCalendarTimestep();
    int getCalendarStep();
    std::string getCurrentDate(const bool isoFormat = true);
    void updateCalendar(const int stepNumber);

    /* Axis */
    void createAxis(const std::string axisId);
    void setAxisSize(const std::string axisId, const size_t size);
    void setAxisValues(const std::string axisId, std::vector<double> values);
    size_t getAxisSize(const std::string axisId);
    std::vector<double> getAxisValues(const std::string axisId);

    /* Domain */
    void createDomain(const std::string domainId);
    void setDomainType(const std::string domainId, const std::string domainType);
    void setDomainGlobalLongitudeSize(const std::string domainId, const size_t size);
    void setDomainGlobalLatitudeSize(const std::string domainId, const size_t size);
    void setDomainLongitudeSize(const std::string domainId, const size_t size);
    void setDomainLatitudeSize(const std::string domainId, const size_t size);
    void setDomainLongitudeStart(const std::string domainId, const size_t start);
    void setDomainLatitudeStart(const std::string domainId, const size_t start);
    void setDomainLongitudeValues(const std::string domainId, std::vector<double> values);
    void setDomainLatitudeValues(const std::string domainId, std::vector<double> values);
    std::string getDomainType(const std::string domainId);
    size_t getDomainGlobalLongitudeSize(const std::string domainId);
    size_t getDomainGlobalLatitudeSize(const std::string domainId);
    size_t getDomainLongitudeSize(const std::string domainId);
    size_t getDomainLatitudeSize(const std::string domainId);
    size_t getDomainLongitudeStart(const std::string domainId);
    size_t getDomainLatitudeStart(const std::string domainId);
    std::vector<double> getDomainLongitudeValues(const std::string domainId);
    std::vector<double> getDomainLatitudeValues(const std::string domainId);

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
    std::string getFieldName(const std::string fieldId);
    std::string getFieldOperation(const std::string fieldId);
    std::string getFieldGridRef(const std::string fieldId);

    /* File */
    void createFile(const std::string fileId);
    void setFileName(const std::string fileId, const std::string fileName);
    void setFileType(const std::string fileId, const std::string fileType);
    void setFileOutputFreq(const std::string fileId, const std::string outputFreq);
    std::string getFileName(const std::string fileId);
    std::string getFileType(const std::string fileId);
    std::string getFileOutputFreq(const std::string fileId);
    bool validFileId(const std::string fileId);
    bool isDefinedFileName(const std::string fileId);
    bool isDefinedFileType(const std::string fileId);
    bool isDefinedFileOutputFreq(const std::string fileId);
    void fileAddField(const std::string fileId, const std::string fieldId);
    std::vector<std::string> fileGetFieldIds(const std::string fileId);

    /* I/O */
    void write(const std::string fieldId, ModelArray& modelarray);

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
    std::string contextId;
    MPI_Comm clientComm;
    MPI_Fint clientComm_F;
    MPI_Fint nullComm_F;
    int mpi_rank { 0 };
    int mpi_size { 0 };

    xios::CCalendarWrapper* clientCalendar;
    std::string convertXiosDatetimeToString(const cxios_date datetime, const bool isoFormat = true);
    cxios_date convertStringToXiosDatetime(const std::string datetime, const bool isoFormat = true);

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
