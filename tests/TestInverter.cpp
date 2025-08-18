#include "converter_models/Topology.h"
#include "json.hpp"
#include <cmath>
#include <UnitTest++.h>

using namespace OpenMagnetics;
using json = nlohmann::json;

// Verify that the ripple estimation logic returns a sensible value for a
// representative operating point.  The numbers used below mirror the Python
// example supplied by the user.  A fairly loose tolerance is used because the
// ripple depends on many parameters (time resolution, FFT length, etc.).
TEST(Test_TwoLevelInverter_CurrentRipple){
    json j;
    // Inverter-level parameters
    json dcBusVoltage; dcBusVoltage["nominal"] = 700.0; j["dcBusVoltage"] = dcBusVoltage;
    json vdcRipple; vdcRipple["nominal"] = 0.0; j["vdcRipple"] = vdcRipple;
    j["inverterLegPowerRating"] = 10000.0/3.0;
    json lineRms; lineRms["nominal"] = 14.7; j["lineRmsCurrent"] = lineRms;
    j["switchingFrequency"] = 10000.0;
    j["riseTime"] = 0.0;
    j["deadtime"] = 0.0;
    j["pwmType"] = "triangular";
    j["modulationStrategy"] = "SPWM";
    j["thirdHarmonicInjectionCoefficient"] = 0.0;
    j["modulationDepth"] = 0.94;

    // Single operating point representing a grid-connected leg delivering power
    json op;
    op["fundamentalFrequency"] = 50.0;
    op["powerFactor"] = 0.98;
    op["outputPower"] = 10000.0/3.0;
    json load;
    load["type"] = "grid";
    json phaseVoltage; phaseVoltage["nominal"] = 400.0/std::sqrt(3.0); load["phaseVoltage"] = phaseVoltage;
    load["gridFrequency"] = 50.0;
    json gridRes; gridRes["nominal"] = 10.3e-3; load["gridResistance"] = gridRes;
    json gridInd; gridInd["nominal"] = 66e-6; load["gridInductance"] = gridInd;
    op["load"] = load;
    j["operatingPoints"] = json::array({op});

    TwoLevelInverter inverter(j);
    L1Params l1; l1.R1_ohm = 0.03; l1.L1_h = 1e-3; // filter inductor parameters

    double ripple = inverter.compute_current_ripple(l1, 0, 1, 20, 1000.0);
    CHECK_CLOSE(3.4, ripple, 0.8); // allow tolerance
}

