/*!
 * @file    Xios.hpp
 * @author  Joe Wallwork <jw2423@cam.ac.uk
 * @date    13 June 2024
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
    void configureServer(std::string calendarType = "Gregorian");

    /* MPI decomposition */
    int getClientMPISize();
    int getClientMPIRank();

    /* Date and duration */
    std::string convertXiosDatetimeToString(cxios_date datetime, bool isoFormat = true);
    void printCXiosDate(cxios_date date);
    void printCXiosDuration(cxios_duration duration);

    /* Calendar */
    void setCalendarType(std::string type);
    void setCalendarOrigin(cxios_date origin);
    void setCalendarStart(cxios_date start);
    void setCalendarTimestep(cxios_duration timestep);
    std::string getCalendarType();
    cxios_date getCalendarOrigin();
    cxios_date getCalendarStart();
    cxios_duration getCalendarTimestep();
    int getCalendarStep();
    std::string getCurrentDate(bool isoFormat = true);
    void updateCalendar(int stepNumber);

    /* Axis */
    void createAxis(std::string axisId);
    void setAxisSize(std::string axisId, size_t size);
    void setAxisValues(std::string axisId, std::vector<double> values);
    size_t getAxisSize(std::string axisId);
    std::vector<double> getAxisValues(std::string axisId);
    bool isDefinedAxisSize(std::string axisId);
    bool areDefinedAxisValues(std::string axisId);

    /* Domain */
    void createDomain(std::string domainId);
    void setDomainType(std::string domainId, std::string domainType);
    void setDomainGlobalLongitudeSize(std::string domainId, size_t size);
    void setDomainGlobalLatitudeSize(std::string domainId, size_t size);
    void setDomainLongitudeSize(std::string domainId, size_t size);
    void setDomainLatitudeSize(std::string domainId, size_t size);
    void setDomainLongitudeStart(std::string domainId, size_t start);
    void setDomainLatitudeStart(std::string domainId, size_t start);
    void setDomainLongitudeValues(std::string domainId, std::vector<double> values);
    void setDomainLatitudeValues(std::string domainId, std::vector<double> values);
    std::string getDomainType(std::string domainId);
    size_t getDomainGlobalLongitudeSize(std::string domainId);
    size_t getDomainGlobalLatitudeSize(std::string domainId);
    size_t getDomainLongitudeSize(std::string domainId);
    size_t getDomainLatitudeSize(std::string domainId);
    size_t getDomainLongitudeStart(std::string domainId);
    size_t getDomainLatitudeStart(std::string domainId);
    std::vector<double> getDomainLongitudeValues(std::string domainId);
    std::vector<double> getDomainLatitudeValues(std::string domainId);
    bool isDefinedDomainType(std::string domainId);
    bool isDefinedDomainGlobalLongitudeSize(std::string domainId);
    bool isDefinedDomainGlobalLatitudeSize(std::string domainId);
    bool isDefinedDomainLongitudeSize(std::string domainId);
    bool isDefinedDomainLatitudeSize(std::string domainId);
    bool isDefinedDomainLongitudeStart(std::string domainId);
    bool isDefinedDomainLatitudeStart(std::string domainId);
    bool areDefinedDomainLongitudeValues(std::string domainId);
    bool areDefinedDomainLatitudeValues(std::string domainId);

    /* Grid */
    void createGrid(std::string gridId);
    void setGridName(std::string gridId, std::string name);
    std::string getGridName(std::string gridId);
    bool isDefinedGridName(std::string GridId);
    void gridAddAxis(std::string axisId, std::string domainId);
    void gridAddDomain(std::string gridId, std::string domainId);

    /* Field */
    void createField(std::string fieldId);
    void setFieldName(std::string fieldId, std::string name);
    void setFieldOperation(std::string fieldId, std::string operation);
    void setFieldGridRef(std::string fieldId, std::string gridRef);
    std::string getFieldName(std::string fieldId);
    std::string getFieldOperation(std::string fieldId);
    std::string getFieldGridRef(std::string fieldId);
    bool isDefinedFieldName(std::string fieldId);
    bool isDefinedFieldOperation(std::string fieldId);
    bool isDefinedFieldGridRef(std::string fieldId);

    /* File */
    void createFile(std::string fileId);
    void setFileName(std::string fileId, std::string fileName);
    void setFileType(std::string fileId, std::string fileType);
    void setFileOutputFreq(std::string fileId, std::string outputFreq);
    std::string getFileName(std::string fileId);
    std::string getFileType(std::string fileId);
    std::string getFileOutputFreq(std::string fileId);
    bool validFileId(std::string fileId);
    bool isDefinedFileName(std::string fileId);
    bool isDefinedFileType(std::string fileId);
    bool isDefinedFileOutputFreq(std::string fileId);
    void fileAddField(std::string fileId, std::string fieldId);

    /* I/O */
    void write(const std::string fieldId, double* data, const size_t ni, const size_t nj);
    void write(
        const std::string fieldId, double* data, const size_t ni, const size_t nj, const size_t nk);

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

    xios::CAxisGroup* getAxisGroup();
    xios::CDomainGroup* getDomainGroup();
    xios::CFieldGroup* getFieldGroup();
    xios::CGridGroup* getGridGroup();
    xios::CFileGroup* getFileGroup();

    xios::CAxis* getAxis(std::string axisId);
    xios::CDomain* getDomain(std::string domainId);
    xios::CField* getField(std::string fieldId);
    xios::CGrid* getGrid(std::string gridId);
    xios::CFile* getFile(std::string fileId);
};

}

#endif // USE_XIOS
#endif
