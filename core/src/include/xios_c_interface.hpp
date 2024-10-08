/*!
 * @file    xios_c_interface.hpp
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @author  Joe Wallwork <jw2423@cam.ac.uk>
 * @date    12 August 2024
 * @brief   C interface for XIOS library
 * @details
 * This interface is based on an earlier version provided by Laurent as part of
 * the https://github.com/nextsimhub/xios_cpp_toy repo. This C interface is
 * designed to connect with the underlying Fortran interface of XIOS 2.
 *
 * This can be expanded as we add more XIOS functionality to the nextSIM-DG XIOS
 * C++ interface `Xios.cpp`.
 *
 */
#ifndef XIOS_C_INTERFACE
#define XIOS_C_INTERFACE

#include <axis.hpp>
#include <calendar_wrapper.hpp>
#include <domain.hpp>
#include <field.hpp>
#include <file.hpp>
#include <grid.hpp>
#include <icdate.hpp>

extern "C" {

// client methods
void cxios_init_client(
    const char* client_id, int len_client_id, MPI_Fint* f_local_comm, MPI_Fint* f_return_comm);
void cxios_finalize();

// context methods
void cxios_context_initialize(const char* context_id, int len_context_id, MPI_Fint* f_comm);
void cxios_context_is_initialized(const char* context_id, int len_context_id, bool* initialized);
void cxios_context_close_definition();
void cxios_context_finalize();

// conversions
void cxios_date_convert_to_string(cxios_date date_c, char* str, int str_size);
cxios_date cxios_date_convert_from_string(const char* str, int str_size);
void cxios_duration_convert_to_string(cxios_duration dur_c, char* str, int str_size);
cxios_duration cxios_duration_convert_from_string(const char* str, int str_size);

// calendar methods
void cxios_create_calendar(xios::CCalendarWrapper* calendar_wrapper_hdl);
void cxios_set_calendar_wrapper_date_time_origin(
    xios::CCalendarWrapper* calendarWrapper_hdl, cxios_date time_origin_c);
void cxios_set_calendar_wrapper_date_start_date(
    xios::CCalendarWrapper* calendarWrapper_hdl, cxios_date start_date_c);
void cxios_set_calendar_wrapper_type(
    xios::CCalendarWrapper* calendarWrapper_hdl, const char* type, int type_size);
void cxios_get_current_calendar_wrapper(xios::CCalendarWrapper** _ret);
void cxios_get_calendar_wrapper_date_start_date(
    xios::CCalendarWrapper* calendarWrapper_hdl, cxios_date* start_date_c);
void cxios_get_calendar_wrapper_date_time_origin(
    xios::CCalendarWrapper* calendarWrapper_hdl, cxios_date* time_origin_c);
void cxios_get_calendar_wrapper_type(
    xios::CCalendarWrapper* calendarWrapper_hdl, const char* type, int type_size);
void cxios_get_current_date(cxios_date* date);
void cxios_update_calendar(int step);

// timestep methods
void cxios_set_calendar_wrapper_timestep(
    xios::CCalendarWrapper* calendar_wrapper_hdl, cxios_duration timestep_c);
void cxios_get_calendar_wrapper_timestep(
    xios::CCalendarWrapper* calendar_wrapper_hdl, cxios_duration* timestep_c);
void cxios_update_calendar_timestep(xios::CCalendarWrapper* calendarWrapper_hdl);

// axis group methods
void cxios_axisgroup_handle_create(xios::CAxisGroup** _ret, const char* _id, int _id_len);
void cxios_xml_tree_add_axis(
    xios::CAxisGroup* axis_grp, xios::CAxis** axis, const char* _id, int _id_len);

// axis methods
void cxios_axis_handle_create(xios::CAxis** _ret, const char* _id, int _id_len);
void cxios_axis_valid_id(bool* _ret, const char* _id, int _id_len);
void cxios_set_axis_n_glo(xios::CAxis* axis_hdl, int n_glo);
void cxios_set_axis_value(xios::CAxis* axis_hdl, double* value, int* extent);
void cxios_get_axis_n_glo(xios::CAxis* axis_hdl, int* n_glo);
void cxios_get_axis_value(xios::CAxis* axis_hdl, double* value, int* extent);
bool cxios_is_defined_axis_n_glo(xios::CAxis* axis_hdl);
bool cxios_is_defined_axis_value(xios::CAxis* axis_hdl);

// domain group methods
void cxios_domaingroup_handle_create(xios::CDomainGroup** _ret, const char* _id, int _id_len);
void cxios_xml_tree_add_domain(
    xios::CDomainGroup* domain_grp, xios::CDomain** domain, const char* _id, int _id_len);

// domain methods
void cxios_domain_handle_create(xios::CDomain** _ret, const char* _id, int _id_len);
void cxios_domain_valid_id(bool* _ret, const char* _id, int _id_len);
void cxios_set_domain_type(xios::CDomain* domain_hdl, const char* type, int type_size);
void cxios_set_domain_ni_glo(xios::CDomain* domain_hdl, int ni_glo);
void cxios_set_domain_nj_glo(xios::CDomain* domain_hdl, int nj_glo);
void cxios_set_domain_ni(xios::CDomain* domain_hdl, int ni);
void cxios_set_domain_nj(xios::CDomain* domain_hdl, int nj);
void cxios_set_domain_ibegin(xios::CDomain* domain_hdl, int ibegin);
void cxios_set_domain_jbegin(xios::CDomain* domain_hdl, int jbegin);
void cxios_set_domain_lonvalue_1d(xios::CDomain* domain_hdl, double* lonvalue_1d, int* extent);
void cxios_set_domain_latvalue_1d(xios::CDomain* domain_hdl, double* latvalue_1d, int* extent);
void cxios_get_domain_type(xios::CDomain* domain_hdl, char* type, int type_size);
void cxios_get_domain_ni_glo(xios::CDomain* domain_hdl, int* ni_glo);
void cxios_get_domain_nj_glo(xios::CDomain* domain_hdl, int* nj_glo);
void cxios_get_domain_ni(xios::CDomain* domain_hdl, int* ni);
void cxios_get_domain_nj(xios::CDomain* domain_hdl, int* nj);
void cxios_get_domain_ibegin(xios::CDomain* domain_hdl, int* ibegin);
void cxios_get_domain_jbegin(xios::CDomain* domain_hdl, int* jbegin);
void cxios_get_domain_lonvalue_1d(xios::CDomain* domain_hdl, double* lonvalue_1d, int* extent);
void cxios_get_domain_latvalue_1d(xios::CDomain* domain_hdl, double* latvalue_1d, int* extent);
bool cxios_is_defined_domain_type(xios::CDomain* axis_hdl);
bool cxios_is_defined_domain_ni_glo(xios::CDomain* axis_hdl);
bool cxios_is_defined_domain_nj_glo(xios::CDomain* axis_hdl);
bool cxios_is_defined_domain_ni(xios::CDomain* axis_hdl);
bool cxios_is_defined_domain_nj(xios::CDomain* axis_hdl);
bool cxios_is_defined_domain_ibegin(xios::CDomain* axis_hdl);
bool cxios_is_defined_domain_jbegin(xios::CDomain* axis_hdl);
bool cxios_is_defined_domain_lonvalue_1d(xios::CDomain* axis_hdl);
bool cxios_is_defined_domain_latvalue_1d(xios::CDomain* axis_hdl);

// grid group methods
void cxios_gridgroup_handle_create(xios::CGridGroup** _ret, const char* _id, int _id_len);
void cxios_xml_tree_add_grid(
    xios::CGridGroup* grid_grp, xios::CGrid** grid, const char* _id, int _id_len);

// grid methods
void cxios_grid_handle_create(xios::CGrid** _ret, const char* _id, int _id_len);
void cxios_grid_valid_id(bool* _ret, const char* _id, int _id_len);
void cxios_set_grid_name(xios::CGrid* _ret, const char* name, int name_size);
void cxios_get_grid_name(xios::CGrid* _ret, char* name, int name_size);
bool cxios_is_defined_grid_name(xios::CGrid* file_hdl);
void cxios_xml_tree_add_axistogrid(
    xios::CGrid* grid, xios::CAxis** axis, const char* _id, int _id_len);
void cxios_xml_tree_add_domaintogrid(
    xios::CGrid* grid, xios::CDomain** domain, const char* _id, int _id_len);

// field group methods
void cxios_fieldgroup_handle_create(xios::CFieldGroup** _ret, const char* _id, int _id_len);
void cxios_xml_tree_add_field(
    xios::CFieldGroup* field_grp, xios::CField** field, const char* _id, int _id_len);

// field methods
void cxios_field_handle_create(xios::CField** _ret, const char* _id, int _id_len);
void cxios_field_valid_id(bool* _ret, const char* _id, int _id_len);
void cxios_set_field_name(xios::CField* _ret, const char* name, int name_size);
void cxios_set_field_operation(xios::CField* _ret, const char* operation, int operation_size);
void cxios_set_field_grid_ref(xios::CField* _ret, const char* grid_ref, int grid_ref_size);
void cxios_set_field_read_access(xios::CField* _ret, bool read_access);
void cxios_set_field_freq_offset(xios::CField* _ret, cxios_duration freq_offset);
void cxios_get_field_name(xios::CField* _ret, char* name, int name_size);
void cxios_get_field_operation(xios::CField* _ret, char* operation, int operation_size);
void cxios_get_field_grid_ref(xios::CField* _ret, char* grid_ref, int grid_ref_size);
void cxios_get_field_read_access(xios::CField* _ret, bool* read_access);
void cxios_get_field_freq_offset(xios::CField* _ret, cxios_duration* freq_offset);
bool cxios_is_defined_field_name(xios::CField* _ret);
bool cxios_is_defined_field_operation(xios::CField* _ret);
bool cxios_is_defined_field_grid_ref(xios::CField* _ret);
bool cxios_is_defined_field_read_access(xios::CField* _ret);
bool cxios_is_defined_field_freq_offset(xios::CField* _ret);

// file group methods
void cxios_filegroup_handle_create(xios::CFileGroup** _ret, const char* _id, int _id_len);
void cxios_xml_tree_add_file(
    xios::CFileGroup* file_grp, xios::CFile** file, const char* _id, int _id_len);

// file methods
void cxios_file_handle_create(xios::CFile** _ret, const char* _id, int _id_len);
void cxios_file_valid_id(bool* _ret, const char* _id, int _id_len);
void cxios_set_file_name(xios::CFile* file_hdl, const char* name, int name_size);
void cxios_set_file_type(xios::CFile* file_hdl, const char* type, int type_size);
void cxios_set_file_output_freq(xios::CFile* file_hdl, cxios_duration output_freq_c);
void cxios_set_file_split_freq(xios::CFile* file_hdl, cxios_duration split_freq_c);
void cxios_set_file_mode(xios::CFile* file_hdl, const char* mode, int mode_size);
void cxios_set_file_par_access(xios::CFile* file_hdl, const char* par_access, int par_access_size);
void cxios_get_file_name(xios::CFile* file_hdl, char* name, int name_size);
void cxios_get_file_type(xios::CFile* file_hdl, char* type, int type_size);
void cxios_get_file_output_freq(xios::CFile* file_hdl, cxios_duration* output_freq_c);
void cxios_get_file_split_freq(xios::CFile* file_hdl, cxios_duration* split_freq_c);
void cxios_get_file_mode(xios::CFile* file_hdl, char* mode, int mode_size);
void cxios_get_file_par_access(xios::CFile* file_hdl, char* par_access, int par_access_size);
bool cxios_is_defined_file_name(xios::CFile* file_hdl);
bool cxios_is_defined_file_type(xios::CFile* file_hdl);
bool cxios_is_defined_file_output_freq(xios::CFile* file_hdl);
bool cxios_is_defined_file_split_freq(xios::CFile* file_hdl);
bool cxios_is_defined_file_mode(xios::CFile* file_hdl);
bool cxios_is_defined_file_par_access(xios::CFile* file_hdl);
void cxios_xml_tree_add_fieldtofile(
    xios::CFile* file, xios::CField** field, const char* _id, int _id_len);

// writing methods
void cxios_write_data_k82(const char* fieldid, int fieldid_size, const double* data_k8,
    int data_size1, int data_size2, int tileid);
void cxios_write_data_k83(const char* fieldid, int fieldid_size, const double* data_k8,
    int data_size1, int data_size2, int data_size3, int tileid);
void cxios_write_data_k84(const char* fieldid, int fieldid_size, const double* data_k8,
    int data_size1, int data_size2, int data_size3, int date_size4, int tileid);

// reading methods
void cxios_read_data_k82(
    const char* fieldid, int fieldid_size, const double* data_k8, int data_size1, int data_size2);
void cxios_read_data_k83(const char* fieldid, int fieldid_size, const double* data_k8,
    int data_size1, int data_size2, int data_size3);
void cxios_read_data_k84(const char* fieldid, int fieldid_size, const double* data_k8,
    int data_size1, int data_size2, int data_size3, int data_size4);
};

#endif
