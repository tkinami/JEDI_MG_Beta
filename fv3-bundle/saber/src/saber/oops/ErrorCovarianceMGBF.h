/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_ERRORCOVARIANCEMGBF_H_
#define SABER_OOPS_ERRORCOVARIANCEMGBF_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

#include "saber/oops/OoMgbf.h"
#include "saber/oops/ParametersMGBF.h"

namespace eckit {
  class LocalConfiguration;
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------

/// Model space error covariance with MGBF

template <typename MODEL>
class ErrorCovarianceMGBF : public oops::ModelSpaceCovarianceBase<MODEL>,
                             public util::Printable,
                             private util::ObjectCounter<ErrorCovarianceMGBF<MODEL>> {
  typedef oops::Geometry<MODEL>    Geometry_;
  typedef oops::Increment<MODEL>   Increment_;
  typedef OoMgbf<MODEL>           OoMgbf_;
  typedef oops::State<MODEL>       State_;
  typedef ParametersMGBF<MODEL>   Parameters_;

 public:
  static const std::string classname() {return "saber::ErrorCovarianceMGBF";}

  ErrorCovarianceMGBF(const Geometry_ &, const oops::Variables &,
                       const eckit::Configuration &, const State_ &, const State_ &);
  virtual ~ErrorCovarianceMGBF();

 private:
  ErrorCovarianceMGBF(const ErrorCovarianceMGBF&);
  ErrorCovarianceMGBF& operator=(const ErrorCovarianceMGBF&);

 private:
  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;

  void print(std::ostream &) const override;

  std::unique_ptr<OoMgbf_> ooMgbf_;
};

// =============================================================================

template<typename MODEL>
ErrorCovarianceMGBF<MODEL>::ErrorCovarianceMGBF(const Geometry_ & resol,
                                                  const oops::Variables & vars,
                                                  const eckit::Configuration & conf,
                                                  const State_ & xb, const State_ & fg)
  : oops::ModelSpaceCovarianceBase<MODEL>(xb, fg, resol, conf), ooMgbf_()
{
  oops::Log::trace() << "ErrorCovarianceMGBF::ErrorCovarianceMGBF starting" << std::endl;

// Setup parameters
  Parameters_ param(resol, vars, xb.validTime(), conf);

// Transfer OoMgbf pointer
  ooMgbf_.reset(new OoMgbf_(param.getOoMgbf()));

  oops::Log::trace() << "ErrorCovarianceMGBF::ErrorCovarianceMGBF done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ErrorCovarianceMGBF<MODEL>::~ErrorCovarianceMGBF() {
  oops::Log::trace() << "ErrorCovarianceMGBF<MODEL>::~ErrorCovarianceMGBF starting" << std::endl;
  util::Timer timer(classname(), "~ErrorCovarianceMGBF");
  oops::Log::trace() << "ErrorCovarianceMGBF<MODEL>::~ErrorCovarianceMGBF done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceMGBF<MODEL>::doRandomize(Increment_ & dx) const {
  oops::Log::trace() << "ErrorCovarianceMGBF<MODEL>::doRandomize starting" << std::endl;
  util::Timer timer(classname(), "doRandomize");
  ooMgbf_->randomize(dx);
  oops::Log::trace() << "ErrorCovarianceMGBF<MODEL>::doRandomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceMGBF<MODEL>::doMultiply(const Increment_ & dxi,
                                             Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovarianceMGBF<MODEL>::doMultiply starting" << std::endl;
  util::Timer timer(classname(), "doMultiply");
  ooMgbf_->multiplyFilter(dxi, dxo);
  oops::Log::trace() << "ErrorCovarianceMGBF<MODEL>::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceMGBF<MODEL>::doInverseMultiply(const Increment_ & dxi,
                                                    Increment_ & dxo) const {
  oops::Log::trace() << "ErrorCovarianceMGBF<MODEL>::doInverseMultiply starting" << std::endl;
  util::Timer timer(classname(), "doInverseMultiply");
  ooMgbf_->inverseMultiplyFilter(dxi, dxo);
  oops::Log::trace() << "ErrorCovarianceMGBF<MODEL>::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ErrorCovarianceMGBF<MODEL>::print(std::ostream & os) const {
  oops::Log::trace() << "ErrorCovarianceMGBF<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << "ErrorCovarianceMGBF<MODEL>::print not implemented";
  oops::Log::trace() << "ErrorCovarianceMGBF<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_ERRORCOVARIANCEMGBF_H_
