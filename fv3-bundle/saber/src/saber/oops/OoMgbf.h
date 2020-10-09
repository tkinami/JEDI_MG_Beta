/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_OOMGBF_H_
#define SABER_OOPS_OOMGBF_H_

#include <algorithm>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/assimilation/GMRESR.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Variables.h"
#if !ATLASIFIED
#include "oops/generic/UnstructuredGrid.h"
#endif
#include "oops/interface/Increment.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "saber/mgbf/type_mgbf.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
  class UnstructuredGrid;
}

namespace saber {

// -----------------------------------------------------------------------------
/// OoMgbf C++ interface

template<typename MODEL> class OoMgbf {
  typedef oops::Geometry<MODEL>    Geometry_;
  typedef oops::Increment<MODEL>   Increment_;

 public:
  OoMgbf(const Geometry_ &, const oops::Variables &, const util::DateTime &,
          const eckit::LocalConfiguration);
  explicit OoMgbf(OoMgbf &);
  ~OoMgbf();

  // C++ interfaces
  size_t getSize() {return keyOoMgbf_.size();}
  int getKey(int igrid) const {return keyOoMgbf_[igrid];}
  void clearKey() {keyOoMgbf_.clear();}

  // Fortran interfaces
  void runDrivers() const;
  void multiplyVbal(const Increment_ &, Increment_ &) const;
  void multiplyVbalInv(const Increment_ &, Increment_ &) const;
  void multiplyVbalAd(const Increment_ &, Increment_ &) const;
  void multiplyVbalInvAd(const Increment_ &, Increment_ &) const;
  void multiplyStdDev(const Increment_ &, Increment_ &) const;
  void multiplyStdDevInv(const Increment_ &, Increment_ &) const;
  void multiplyFilter(Increment_ &) const;
  void multiplyFilter(const Increment_ &, Increment_ &) const;
  void inverseMultiplyFilter(const Increment_ &, Increment_ &) const;
  void randomize(Increment_ &) const;
  void getParameter(const std::string &, Increment_ &) const;
  void setParameter(const std::string &, const Increment_ &) const;

  // Aliases for inversion with GMRESR
  void multiply(const Increment_ & dxi, Increment_ & dxo) const {multiplyFilter(dxi, dxo);}

 private:
  std::vector<int> keyOoMgbf_;
#if !ATLASIFIED
  std::unique_ptr<oops::UnstructuredGrid> ug_;
#endif
};

// -----------------------------------------------------------------------------
template<typename MODEL>
OoMgbf<MODEL>::OoMgbf(const Geometry_ & resol,
                        const oops::Variables & vars,
                        const util::DateTime & time,
                        const eckit::LocalConfiguration conf) : keyOoMgbf_() {
  // Grids
  std::vector<eckit::LocalConfiguration> grids;

  // Get global prefix
  std::string prefix;
  conf.get("prefix", prefix);

#if ATLASIFIED
  // Get the grids configuration from input configuration and complete it
  if (conf.has("grids")) {
    // Get grids from input configuration
    conf.get("grids", grids);
    ASSERT(grids.size() > 0);
  } else {
    // Create one empty configuration
    eckit::LocalConfiguration emptyConf;
    grids.push_back(emptyConf);
  }

  // Loop over grids
  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    // Add prefix
    if (!grids[jgrid].has("prefix")) {
      std::ostringstream ss;
      ss << std::setw(2) << std::setfill('0') << jgrid;
      grids[jgrid].set("prefix", prefix + "_" + ss.str());
    }

    // Add input variables to the grid configuration
    std::vector<std::string> vars_str;
    if (grids[jgrid].has("variables")) {
      grids[jgrid].get("variables", vars_str);
    } else {
      for (unsigned int jvar = 0; jvar < vars.size(); ++jvar) {
        vars_str.push_back(vars[jvar]);
      }
      grids[jgrid].set("variables", vars_str);
    }
    grids[jgrid].set("nv", vars_str.size());

    // Get the required number of levels add it to the grid configuration
    Increment_ dx(resol, vars, time);
    std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
    dx.setAtlas(atlasFieldSet.get());
    int nl = 0;
    for (unsigned int jvar = 0; jvar < vars_str.size(); ++jvar) {
      atlas::Field atlasField = atlasFieldSet->field(vars_str[jvar]);
      nl = std::max(nl, atlasField.levels());
    }
    grids[jgrid].set("nl", nl);

    // Add level index for 2D fields (first or last, first by default)
    if (!grids[jgrid].has("lev2d")) {
      grids[jgrid].set("lev2d", "first");
    }
  }
#else
  // Get the grids configuration from the unstructured grid configuration
  Increment_ dx(resol, vars, time);
  int colocated = 1;
  if (conf.has("colocated")) {
    colocated = conf.getInt("colocated");
  }
  ug_.reset(new oops::UnstructuredGrid(colocated, 1));
  dx.ug_coord(*ug_.get());
  ug_->defineGeometry();
  ug_->defineGrids(grids);

  // Modify grids
  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    int grid_index;
    grids[jgrid].get("grid_index", grid_index);
    std::ostringstream ss;
    ss << std::setw(2) << std::setfill('0') << grid_index;
    grids[jgrid].set("prefix", prefix + "_" + ss.str());
    std::vector<std::string> vars_str;
    grids[jgrid].get("variables", vars_str);
    grids[jgrid].set("nv", vars_str.size());
  }
#endif

  // Check grids number
  ASSERT(grids.size() > 0);

  // Print configuration
  oops::Log::info() << "Configuration: " << conf << std::endl;

  for (unsigned int jgrid = 0; jgrid < grids.size(); ++jgrid) {
    // Print configuration for this grid
    oops::Log::info() << "Grid " << jgrid << ": " << grids[jgrid] << std::endl;

    // Create OoMgbf instance
    int keyOoMgbf = 0;
#if ATLASIFIED
    mgbf_create_f90(keyOoMgbf, &resol.getComm(),
                    resol.atlasFunctionSpace()->get(),
                    resol.atlasFieldSet()->get(),
                    conf, grids[jgrid]);
#else
    mgbf_create_f90(keyOoMgbf, &resol.getComm(),
                    ug_->atlasFunctionSpace()->get(),
                    ug_->atlasFieldSet()->get(),
                    conf, grids[jgrid]);
#endif
    keyOoMgbf_.push_back(keyOoMgbf);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
OoMgbf<MODEL>::OoMgbf(OoMgbf & other) : keyOoMgbf_() {
  for (unsigned int jgrid = 0; jgrid < other.getSize(); ++jgrid) {
    keyOoMgbf_.push_back(other.getKey(jgrid));
  }
  other.clearKey();
#if !ATLASIFIED
  ug_ = std::move(other.ug_);
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
OoMgbf<MODEL>::~OoMgbf() {
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    if (keyOoMgbf_[jgrid] > 0) mgbf_dealloc_f90(keyOoMgbf_[jgrid]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::runDrivers() const {
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_run_drivers_f90(keyOoMgbf_[jgrid]);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::multiplyVbal(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_apply_vbal_f90(keyOoMgbf_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::multiplyVbalInv(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_apply_vbal_inv_f90(keyOoMgbf_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::multiplyVbalAd(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_apply_vbal_ad_f90(keyOoMgbf_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::multiplyVbalInvAd(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_apply_vbal_inv_ad_f90(keyOoMgbf_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::multiplyStdDev(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_apply_stddev_f90(keyOoMgbf_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::multiplyStdDevInv(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_apply_stddev_inv_f90(keyOoMgbf_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::multiplyFilter(Increment_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_apply_filter_f90(keyOoMgbf_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::multiplyFilter(const Increment_ & dxi, Increment_ & dxo) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dxi.toAtlas(atlasFieldSet.get());
#else
  dxi.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_apply_filter_f90(keyOoMgbf_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dxo.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dxo.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::inverseMultiplyFilter(const Increment_ & dxi, Increment_ & dxo) const {
  oops::IdentityMatrix<Increment_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::randomize(Increment_ & dx) const {
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_randomize_f90(keyOoMgbf_[jgrid], atlasFieldSet->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::getParameter(const std::string & param, Increment_ & dx) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.setAtlas(atlasFieldSet.get());
#else
  ug_->setAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_get_parameter_f90(keyOoMgbf_[jgrid], nstr, cstr, atlasFieldSet.get()->get());
  }
#if ATLASIFIED
  dx.fromAtlas(atlasFieldSet.get());
#else
  ug_->fromAtlas(atlasFieldSet.get());
  dx.field_from_ug(*ug_.get());
#endif
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void OoMgbf<MODEL>::setParameter(const std::string & param, const Increment_ & dx) const {
  const int nstr = param.size();
  const char *cstr = param.c_str();
  std::unique_ptr<atlas::FieldSet> atlasFieldSet(new atlas::FieldSet());
#if ATLASIFIED
  dx.toAtlas(atlasFieldSet.get());
#else
  dx.field_to_ug(*ug_.get());
  ug_->toAtlas(atlasFieldSet.get());
#endif
  for (unsigned int jgrid = 0; jgrid < keyOoMgbf_.size(); ++jgrid) {
    mgbf_set_parameter_f90(keyOoMgbf_[jgrid], nstr, cstr, atlasFieldSet->get());
  }
}
// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_OOMGBF_H_
