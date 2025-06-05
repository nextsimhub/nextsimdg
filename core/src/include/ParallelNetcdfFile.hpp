#include <mpi.h>
#include <ncFile.h>
#include <ncVar.h>

namespace netCDF {

class NcFilePar : public netCDF::NcFile {
public:
    NcFilePar() = default;

    NcFilePar(const std::string& filePath, const FileMode fMode, MPI_Comm comm,
        MPI_Info info = MPI_INFO_NULL);

    void open_par(const std::string& path, const FileMode fMode, MPI_Comm comm, MPI_Info info);
};

void setVariableCollective(NcVar var, NcGroup group);

}
