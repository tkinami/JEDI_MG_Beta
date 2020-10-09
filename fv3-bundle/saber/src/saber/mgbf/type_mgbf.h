/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_MGBF_TYPE_MGBF_H_
#define SABER_MGBF_TYPE_MGBF_H_

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace saber {
extern "C" {
  void mgbf_create_f90(int &, const eckit::mpi::Comm *,
                       atlas::functionspace::FunctionSpaceImpl *, atlas::field::FieldSetImpl *,
                       const eckit::Configuration &,
                       const eckit::Configuration &);
  void mgbf_run_drivers_f90(const int &);
  void mgbf_apply_vbal_f90(const int &, atlas::field::FieldSetImpl *);
  void mgbf_apply_vbal_inv_f90(const int &, atlas::field::FieldSetImpl *);
  void mgbf_apply_vbal_ad_f90(const int &, atlas::field::FieldSetImpl *);
  void mgbf_apply_vbal_inv_ad_f90(const int &, atlas::field::FieldSetImpl *);
  void mgbf_apply_stddev_f90(const int &, atlas::field::FieldSetImpl *);
  void mgbf_apply_stddev_inv_f90(const int &, atlas::field::FieldSetImpl *);
  void mgbf_apply_filter_f90(const int &, atlas::field::FieldSetImpl *);
  void mgbf_get_cv_size_f90(const int &, int &);
  void mgbf_apply_filter_sqrt_f90(const int &, const double *, atlas::field::FieldSetImpl *);
  void mgbf_apply_filter_sqrt_ad_f90(const int &, atlas::field::FieldSetImpl *, const double *);
  void mgbf_randomize_f90(const int &, atlas::field::FieldSetImpl *);
  void mgbf_get_parameter_f90(const int &, const int &, const char *, atlas::field::FieldSetImpl *);
  void mgbf_set_parameter_f90(const int &, const int &, const char *, atlas::field::FieldSetImpl *);
  void mgbf_dealloc_f90(const int &);
}

}  // namespace saber

#endif  // SABER_MGBF_TYPE_MGBF_H_
