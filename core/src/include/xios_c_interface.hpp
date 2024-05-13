/*!
 * @file    xios_c_interface.hpp
 * @author  Tom Meltzer <tdm39@cam.ac.uk>
 * @date    Fri 23 Feb 13:43:16 GMT 2024
 * @brief   C interface for xios library
 * @details
 * This interface is based on an earlier version provided by Laurent as part of
 * the https://github.com/nextsimhub/xios_cpp_toy repo. This C interface is
 * designed to connect with the underlying Fortran interface of XIOS 2.
 *
 * This can be expanded as we add more XIOS functionality to the nextsim XIOS
 * C++ interface `Xios.cpp`.
 *
 */
#ifndef XIOS_C_INTERFACE
#define XIOS_C_INTERFACE

#include <axis.hpp>
#include <calendar_wrapper.hpp>
#include <domain.hpp>
#include <field.hpp>
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

void cxios_get_current_calendar_wrapper(xios::CCalendarWrapper** _ret);

// datetime conversions
void cxios_date_convert_to_string(cxios_date date_c, char* str, int str_size);
cxios_date cxios_date_convert_from_string(const char* str, int str_size);

void cxios_duration_convert_to_string(cxios_duration dur_c, char* str, int str_size);

// calendar methods
void cxios_set_calendar_wrapper_date_time_origin(
    xios::CCalendarWrapper* calendarWrapper_hdl, cxios_date time_origin_c);
void cxios_get_calendar_wrapper_date_time_origin(
    xios::CCalendarWrapper* calendarWrapper_hdl, cxios_date* time_origin_c);
void cxios_set_calendar_wrapper_date_start_date(
    xios::CCalendarWrapper* calendarWrapper_hdl, cxios_date start_date_c);
void cxios_get_calendar_wrapper_date_start_date(
    xios::CCalendarWrapper* calendarWrapper_hdl, cxios_date* start_date_c);
void cxios_get_current_date(cxios_date* date);

void cxios_update_calendar(int step);

// timestep methods
void cxios_set_calendar_wrapper_timestep(
    xios::CCalendarWrapper* calendar_wrapper_hdl, cxios_duration timestep_c);
void cxios_get_calendar_wrapper_timestep(
    xios::CCalendarWrapper* calendar_wrapper_hdl, cxios_duration* timestep_c);
void cxios_update_calendar_timestep(xios::CCalendarWrapper* calendarWrapper_hdl);

// grid methods
void cxios_axis_handle_create(xios::CAxis** _ret, const char* _id, int _id_len);
void cxios_get_axis_n_glo(xios::CAxis* axis_hdl, int* n_glo);
void cxios_get_axis_value(xios::CAxis* axis_hdl, double* value, int* extent);

// domain methods
void cxios_domain_handle_create(xios::CDomain** _ret, const char* _id, int _id_len);
void cxios_get_domain_type(xios::CDomain* domain_hdl, char* type, int type_size);
void cxios_get_domain_ni_glo(xios::CDomain* domain_hdl, int* ni_glo);
void cxios_get_domain_nj_glo(xios::CDomain* domain_hdl, int* nj_glo);
void cxios_set_domain_ni(xios::CDomain* domain_hdl, int ni);
void cxios_set_domain_nj(xios::CDomain* domain_hdl, int nj);
void cxios_set_domain_ibegin(xios::CDomain* domain_hdl, int ibegin);
void cxios_set_domain_jbegin(xios::CDomain* domain_hdl, int jbegin);
void cxios_set_domain_lonvalue_1d(xios::CDomain* domain_hdl, double* lonvalue_1d, int* extent);
void cxios_set_domain_latvalue_1d(xios::CDomain* domain_hdl, double* latvalue_1d, int* extent);

// writing methods
void cxios_write_data_k82(const char* fieldid, int fieldid_size, const double* data_k8,
    int data_size1, int data_size2, int tileid);
void cxios_write_data_k83(const char* fieldid, int fieldid_size, const double* data_k8,
    int data_size1, int data_size2, int data_size3, int tileid);
void cxios_write_data_k84(const char* fieldid, int fieldid_size, const double* data_k8,
    int data_size1, int data_size2, int data_size3, int data_size4, int tileid);
};

#endif
