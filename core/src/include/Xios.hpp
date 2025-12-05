/*!
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @author  Adeleke Bankole <ab3191@cam.ac.uk>
 * @brief   XIOS interface header
 * @details
 * Header file for XIOS interface
 */
#ifndef SRC_INCLUDE_XIOS_HPP
#define SRC_INCLUDE_XIOS_HPP

#include "date.hpp"
#if USE_XIOS

#include "include/Configured.hpp"
#include "include/Logged.hpp"
#include "include/ModelArray.hpp"
#include "include/Time.hpp"
#include <boost/date_time/posix_time/posix_time.hpp>
#include <boost/format.hpp>
#include <boost/format/group.hpp>
#include <include/xios_c_interface.hpp>
#include <mpi.h>

namespace Nextsim {

// Forward declarations to avoid circular dependencies
class ParaGridIO;

class Xios : public Configured<Xios> {
private:
    //! Private constructor
    Xios();

    //! Performs some one-time initialization for the class. Returns true.
    static bool doOnce();

public:
    ~Xios();

    //! Prevent copying
    Xios(const Xios&) = delete;

    //! Define Xios handler Singleton
    inline static Xios& getInstance()
    {
        static Xios instance = Xios();
        return instance;
    };

    /* Help config */
    static HelpMap& getHelpText(HelpMap& map, bool getAll);
    static HelpMap& getHelpRecursive(HelpMap& map, bool getAll);

    /* Initialization and finalization */
    void close_context_definition();
    void context_finalize();
    static void finalize();
    bool isInitialized();

    /* Configuration */
    void configure() override;

    /* Calendar, date and duration */
    void setCalendarStart(const TimePoint& start);
    void setCalendarStep(const int stepNumber);
    void incrementCalendar();
    TimePoint getCalendarStart();
    int getCalendarStep();
    TimePoint getCurrentDate();

    /* Axis */
    size_t getAxisSize(const std::string& axisId);

    /* Field */
    void createField(const std::string& fieldId);
    void setFieldGridRef(const std::string& fieldId, const std::string& gridRef);
    void setFieldFreqOffset(const std::string& fieldId, const Duration& freqOffset);
    std::string getFieldGridRef(const std::string& fieldId);
    bool getFieldReadAccess(const std::string& fieldId);
    Duration getFieldFreqOffset(const std::string& fieldId);
    std::set<std::string> configGetForcingFieldNames();
    ModelArray::Type getFieldType(const std::string& fieldId);
    void setFieldType(const std::string& fieldId, const ModelArray::Type& type);

    /* File */
    void createFile(const std::string& fileId, const int& fieldType);
    void fileAddField(const std::string& fileId, const std::string& fieldId);
    std::vector<std::string> fileGetFieldIds(const std::string& fileId);

    enum {
        ENABLED_KEY,
        OUTPUT_FIELD_NAMES_KEY,
        INPUT_FIELD_NAMES_KEY,
        DIAGNOSTIC_PERIOD_KEY,
        DIAGNOSTIC_FILE_KEY,
        DIAGNOSTIC_FIELD_NAMES_KEY,
        FORCING_PERIOD_KEY,
        FORCING_FILE_KEY,
        FORCING_FIELD_NAMES_KEY,
    };

    // TODO: Make the following enum private
    enum {
        OUTPUT_RESTART,
        INPUT_RESTART,
        DIAGNOSTIC,
        FORCING,
    };

protected:
    bool isConfigured;

private:
    inline static bool isEnabled = false;
    MPI_Comm clientComm;
    MPI_Fint clientComm_F;
    MPI_Fint nullComm_F;
    int mpi_rank { 0 };
    int mpi_size { 0 };
    int cStrLen { 20 }; // Length of C-strings passed to XIOS

    /* Client */
    const std::string clientId = "client";
    void setupClient();

    /* Context */
    const std::string contextId = "nextSIM-DG";
    void setupContext();

    /* Calendar, date and duration */
    xios::CCalendarWrapper* clientCalendar;
    std::string convertXiosDatetimeToString(
        const cxios_date& datetime, const bool isoFormat = true);
    cxios_date convertStringToXiosDatetime(std::string datetime, const bool isoFormat = true);
    std::string convertCStrToCppStr(const char* cStr, int cStrLen);
    Duration convertDurationFromXios(const cxios_duration& duration);
    cxios_duration convertDurationToXios(const Duration& duration);
    void setupCalendar();

    /* Axis */
    xios::CAxis* getAxis(const std::string& axisId);

    /* Domain */
    // NOTE: Dimension names get processed as <dim>_<domainId> by XIOS, so we define the domainIds
    //       so that these coincide with the altnames when applied to dimensions x and y.
    std::map<ModelArray::Type, std::string> domainIds = {
        // Standard cell-based x- and y-dimensions (alt. names x_dim and y_dim)
        { ModelArray::Type::H, "dim" },
        { ModelArray::Type::DG, "dim" },
        { ModelArray::Type::DGSTRESS, "dim" },
        // Vertex-based x- and y-dimensions (alt. names x_vertex and y_vertex)
        { ModelArray::Type::VERTEX, "vertex" },
        // CG-based x- and y-dimensions (alt. names x_cg and y_cg)
        { ModelArray::Type::CG, "cg" },
    };
    xios::CDomainGroup* getDomainGroup();
    xios::CDomain* getDomain(const std::string& domainId);
    void setupDomains();

    /* Grid */
    xios::CGrid* getGrid(const std::string& gridId);
    std::map<ModelArray::Type, std::string> gridIds = {
        { ModelArray::Type::H, "HGrid" },
        { ModelArray::Type::VERTEX, "VertexGrid" },
        { ModelArray::Type::DG, "DGGrid" },
        { ModelArray::Type::DGSTRESS, "DGSGrid" },
        { ModelArray::Type::CG, "CGGrid" },
    };
    void setupGrids();

    /* Field */
    xios::CFieldGroup* getFieldGroup();
    xios::CField* getField(const std::string& fieldId);
    std::set<std::string> configGetOutputRestartFieldNames();
    std::set<std::string> configGetInputRestartFieldNames();
    std::set<std::string> configGetDiagnosticFieldNames();
    bool configCheckField(const std::string& fieldId, const bool& readAccess);
    void setFieldReadAccess(const std::string& fieldId, const bool& readAccess);
    std::map<std::string, ModelArray::Type> fieldTypes;
    void setupFields();

    /* File */
    xios::CFileGroup* getFileGroup();
    xios::CFile* getFile(const std::string& fileId);
    std::string outputFileId;
    std::string inputFileId;
    std::string diagnosticFilename;
    std::string diagnosticFileId;
    std::string forcingFilename;
    std::string forcingFileId;
    const std::map<int, std::string&> fileMap = {
        { OUTPUT_RESTART, outputFileId },
        { INPUT_RESTART, inputFileId },
        { DIAGNOSTIC, diagnosticFileId },
        { FORCING, forcingFileId },
    };
    void setupFiles();

    /* I/O */
    void read(const std::string& fieldId, ModelArray& modelarray);
    void write(const std::string& fieldId, const ModelArray& modelarray);

    /* Declare any classes that need to access private members */
    friend ParaGridIO;
};
}

#endif // USE_XIOS
#endif
