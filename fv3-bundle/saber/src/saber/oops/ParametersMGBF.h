/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef SABER_OOPS_PARAMETERSMGBF_H_
#define SABER_OOPS_PARAMETERSMGBF_H_

#include <fstream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

#include "saber/oops/OoMgbf.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace saber {

// -----------------------------------------------------------------------------
/// MGBF parameters

template<typename MODEL>
class ParametersMGBF {
  typedef oops::Geometry<MODEL>                           Geometry_;
  typedef oops::Increment<MODEL>                          Increment_;
  typedef OoMgbf<MODEL>                                  OoMgbf_;
  typedef oops::State<MODEL>                              State_;

 public:
  static const std::string classname() {return "oops::ParametersMGBF";}
  ParametersMGBF(const Geometry_ &,
                  const oops::Variables &,
                  const util::DateTime &,
                  const eckit::Configuration &);
  ~ParametersMGBF();

  OoMgbf_ & getOoMgbf() {return *ooMgbf_;}
  void write() const;

 private:
  const Geometry_ resol_;
  const oops::Variables vars_;
  util::DateTime time_;
  const eckit::LocalConfiguration conf_;
  std::unique_ptr<OoMgbf_> ooMgbf_;
};

// =============================================================================

template<typename MODEL>
ParametersMGBF<MODEL>::ParametersMGBF(const Geometry_ & resol,
                                        const oops::Variables & vars,
                                        const util::DateTime & time,
                                        const eckit::Configuration & conf)
  : resol_(resol), vars_(vars), time_(time), conf_(conf), ooMgbf_()
{
  oops::Log::trace() << "ParametersMGBF<MODEL>::ParametersMGBF construction starting" << std::endl;
  util::Timer timer(classname(), "ParametersMGBF");

  // Setup MGBF configuration
  eckit::LocalConfiguration MGBFConf(conf_, "mgbf");

  // Get missing value
  const double msvalr = util::missingValue(msvalr);
  MGBFConf.set("msvalr", msvalr);

  // Create MGBF
  oops::Log::info() << "Create MGBF" << std::endl;
  ooMgbf_.reset(new OoMgbf_(resol, vars, time_, MGBFConf));

// Read data from files
  oops::Log::info() << "Read data from files" << std::endl;
  if (conf_.has("input")) {
  // Set MGBF input parameters
    std::vector<eckit::LocalConfiguration> inputConfs;
    conf_.get("input", inputConfs);

    for (const auto & inputConf : inputConfs) {
    // Read parameter for the specified timeslot
      const util::DateTime date(inputConf.getString("date"));

    // Setup increment
      Increment_ dx(resol_, vars_, time_);
      dx.read(inputConf);

    // Set parameter to MGBF
      std::string param = inputConf.getString("parameter");
      ooMgbf_->setParameter(param, dx);
    }
  }

// Estimate parameters
//  ooMgbf_->runDrivers();

  oops::Log::trace() << "ParametersMGBF:ParametersMGBF constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
ParametersMGBF<MODEL>::~ParametersMGBF() {
  oops::Log::trace() << "ParametersMGBF<MODEL>::~ParametersMGBF destruction starting" << std::endl;
  util::Timer timer(classname(), "~ParametersMGBF");
  oops::Log::trace() << "ParametersMGBF:~ParametersMGBF destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ParametersMGBF<MODEL>::write() const {
  oops::Log::trace() << "ParametersMGBF::write starting" << std::endl;
  util::Timer timer(classname(), "write");

// Write parameters
  std::vector<eckit::LocalConfiguration> outputConfs;
  conf_.get("output", outputConfs);
  for (const auto & outputConf : outputConfs) {
  // Setup dummy increment
    Increment_ dx(resol_, vars_, time_);
    dx.zero();

  // Get parameter from MGBF
    std::string param = outputConf.getString("parameter");
    ooMgbf_->getParameter(param, dx);

  // Write parameter for the specified timeslot
    const util::DateTime date(outputConf.getString("date"));
    dx.write(outputConf);
    oops::Log::test() << "Norm of " << param << " at " << date << ": " << std::scientific
                      << std::setprecision(3) << dx.norm() << std::endl;
  }
  oops::Log::trace() << "ParametersMGBF::write done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace saber

#endif  // SABER_OOPS_PARAMETERSMGBF_H_
