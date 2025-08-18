#pragma once
#include "Constants.h"
#include <MAS.hpp>
#include "processors/Inputs.h"
#include "constructive_models/Magnetic.h"
#include "json.hpp"
#include "Definitions.h"
#include <vector>
#include <string>
#include <optional>
#include <complex>

using namespace MAS;

namespace OpenMagnetics {

class Topology {
private:
protected:
public:
    Topology() = default;
    ~Topology() = default;
};

class FlybackOperatingPoint : public MAS::FlybackOperatingPoint{
public:
    double resolve_switching_frequency(double inputVoltage, double diodeVoltageDrop, std::optional<double> inductance = std::nullopt, std::optional<std::vector<double>> turnsRatios = std::nullopt, double efficiency = 0.85);
    FlybackModes resolve_mode(std::optional<double> currentRippleRatio = std::nullopt);
};


class Flyback : public MAS::Flyback {
private:
    std::optional<double> maximumDrainSourceVoltage = 600;
    std::vector<OpenMagnetics::FlybackOperatingPoint> operatingPoints;
    std::optional<double> maximumDutyCycle = 0.5;
    double efficiency = 1;

public:
    const std::vector<OpenMagnetics::FlybackOperatingPoint> & get_operating_points() const { return operatingPoints; }
    std::vector<OpenMagnetics::FlybackOperatingPoint> & get_mutable_operating_points() { return operatingPoints; }
    void set_operating_points(const std::vector<OpenMagnetics::FlybackOperatingPoint> & value) { this->operatingPoints = value; }


protected:
public:
    bool _assertErrors = false;

    Flyback(const json& j);
    Flyback() {
        maximumDrainSourceVoltage = 600;
        maximumDutyCycle = 0.5;
        efficiency = 1;
    };

    bool run_checks(bool assert = false);

    // According to Worked Example (7), pages 135-144 — Designing the Flyback Transformer of Switching Power Supplies A - Z (Second Edition) by Sanjaya Maniktala
    Inputs process();
    DesignRequirements process_design_requirements();
    std::vector<OperatingPoint> process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance);
    std::vector<OperatingPoint> process_operating_points(Magnetic magnetic);

    OperatingPoint processOperatingPointsForInputVoltage(double inputVoltage, OpenMagnetics::FlybackOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, std::optional<FlybackModes> customMode=std::nullopt, std::optional<double> customDutyCycle=std::nullopt, std::optional<double> customDeadTime=std::nullopt);
    static double get_total_input_power(std::vector<double> outputCurrents, std::vector<double> outputVoltages, double efficiency, double diodeVoltageDrop);
    static double get_total_input_power(double outputCurrent, double outputVoltage, double efficiency, double diodeVoltageDrop);
    static double get_minimum_output_reflected_voltage(double maximumDrainSourceVoltage, double maximumInputVoltage, double safetyMargin=0.85);

};

class AdvancedFlyback : public Flyback {
private:
    std::vector<double> desiredTurnsRatios;
    double desiredInductance;
    std::vector<std::vector<double>> desiredDutyCycle;
    std::optional<std::vector<double>> desiredDeadTime;

protected:
public:
    bool _assertErrors = false;

    AdvancedFlyback() = default;
    ~AdvancedFlyback() = default;

    AdvancedFlyback(const json& j);

    Inputs process();

    const double & get_desired_inductance() const { return desiredInductance; }
    double & get_mutable_desired_inductance() { return desiredInductance; }
    void set_desired_inductance(const double & value) { this->desiredInductance = value; }

    const std::vector<std::vector<double>> & get_desired_duty_cycle() const { return desiredDutyCycle; }
    std::vector<std::vector<double>> & get_mutable_desired_duty_cycle() { return desiredDutyCycle; }
    void set_desired_duty_cycle(const std::vector<std::vector<double>> & value) { this->desiredDutyCycle = value; }

    std::optional<std::vector<double>> get_desired_dead_time() const { return desiredDeadTime; }
    void set_desired_dead_time(std::optional<std::vector<double>> value) { this->desiredDeadTime = value; }

    const std::vector<double> & get_desired_turns_ratios() const { return desiredTurnsRatios; }
    std::vector<double> & get_mutable_desired_turns_ratios() { return desiredTurnsRatios; }
    void set_desired_turns_ratios(const std::vector<double> & value) { this->desiredTurnsRatios = value; }

};


void from_json(const json & j, FlybackOperatingPoint & x);
void to_json(json & j, const FlybackOperatingPoint & x);

inline void from_json(const json & j, FlybackOperatingPoint& x) {
    x.set_output_voltages(j.at("outputVoltages").get<std::vector<double>>());
    x.set_output_currents(j.at("outputCurrents").get<std::vector<double>>());
    x.set_switching_frequency(get_stack_optional<double>(j, "switchingFrequency"));
    x.set_mode(get_stack_optional<FlybackModes>(j, "mode"));
    x.set_ambient_temperature(j.at("ambientTemperature").get<double>());
}

inline void to_json(json & j, const FlybackOperatingPoint & x) {
    j = json::object();
    j["outputVoltages"] = x.get_output_voltages();
    j["outputCurrents"] = x.get_output_currents();
    j["switchingFrequency"] = x.get_switching_frequency();
    j["mode"] = x.get_mode();
    j["ambientTemperature"] = x.get_ambient_temperature();
}

void from_json(const json & j, Flyback & x);
void to_json(json & j, const Flyback & x);

inline void from_json(const json & j, Flyback& x) {
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_maximum_drain_source_voltage(get_stack_optional<double>(j, "maximumDrainSourceVoltage"));
    x.set_maximum_duty_cycle(get_stack_optional<double>(j, "maximumDutyCycle"));
    x.set_current_ripple_ratio(j.at("currentRippleRatio").get<double>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<FlybackOperatingPoint>>());
    x.set_efficiency(j.at("efficiency").get<double>());
}

inline void to_json(json & j, const Flyback & x) {
    j = json::object();
    j["inputVoltage"] = x.get_input_voltage();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["maximumDrainSourceVoltage"] = x.get_maximum_drain_source_voltage();
    j["maximumDutyCycle"] = x.get_maximum_duty_cycle();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
    j["operatingPoints"] = x.get_operating_points();
    j["efficiency"] = x.get_efficiency();
}

void from_json(const json & j, AdvancedFlyback & x);
void to_json(json & j, const AdvancedFlyback & x);

inline void from_json(const json & j, AdvancedFlyback& x) {
    x.set_input_voltage(j.at("inputVoltage").get<DimensionWithTolerance>());
    x.set_diode_voltage_drop(j.at("diodeVoltageDrop").get<double>());
    x.set_desired_inductance(j.at("desiredInductance").get<double>());
    x.set_desired_dead_time(get_stack_optional<std::vector<double>>(j, "desiredDeadTime"));
    x.set_desired_duty_cycle(j.at("desiredDutyCycle").get<std::vector<std::vector<double>>>());
    x.set_desired_turns_ratios(j.at("desiredTurnsRatios").get<std::vector<double>>());
    x.set_operating_points(j.at("operatingPoints").get<std::vector<FlybackOperatingPoint>>());
    x.set_efficiency(j.at("efficiency").get<double>());
    x.set_current_ripple_ratio(std::numeric_limits<double>::quiet_NaN());
}

inline void to_json(json & j, const AdvancedFlyback & x) {
    j = json::object();
    j["inputVoltage"] = x.get_input_voltage();
    j["diodeVoltageDrop"] = x.get_diode_voltage_drop();
    j["desiredInductance"] = x.get_desired_inductance();
    j["desiredDutyCycle"] = x.get_desired_duty_cycle();
    j["desiredDeadTime"] = x.get_desired_dead_time();
    j["desiredTurnsRatios"] = x.get_desired_turns_ratios();
    j["operatingPoints"] = x.get_operating_points();
    j["efficiency"] = x.get_efficiency();
    j["currentRippleRatio"] = x.get_current_ripple_ratio();
}
} // namespace OpenMagnetics

namespace OpenMagnetics {
// Simple description of the external grid seen by the converter
// VLL_rms: line-to-line RMS voltage of the grid
// f_hz:    fundamental frequency in Hertz
// Rgrid_ohm / Lgrid_h: series impedance modelling the grid
struct GridParams {
    double VLL_rms;    // [V] line-to-line RMS voltage
    double f_hz;       // [Hz] fundamental frequency
    double Rgrid_ohm;  // [Ohm] grid resistance
    double Lgrid_h;    // [H] grid inductance
};

// Parameters of the first inductor of an L or LCL output filter
// Only the resistance and inductance are required for the ripple
// estimation carried out by the TwoLevelInverter class
struct L1Params {
    double R1_ohm = 0.0; // [Ohm] series resistance
    double L1_h = 0.0;   // [H] inductance value
};

//------------------------------------------------------------------------------
// TwoLevelInverter
//------------------------------------------------------------------------------
// Convenience C++ port of the reference Python script supplied by the user.
// The class stores all the information needed to reproduce the switching
// waveform of a single inverter leg and to estimate the high–frequency ripple
// of the series inductor current.  Only the parameters requested in the JSON
// schema are stored so that the class mirrors the user's input structure.
class TwoLevelInverter {
public:
    // Description of the electrical load connected to the inverter.  Depending
    // on the "type" field the relevant members are used:
    //  - "grid": phaseVoltage, gridFrequency, gridResistance, gridInductance
    //  - "R-L":  resistance, inductance
    struct Load {
        std::string type;      // "grid" or "R-L"
        double phaseVoltage = 0.0;   // [V] phase RMS voltage when grid connected
        double gridFrequency = 0.0;  // [Hz] grid fundamental frequency
        double gridResistance = 0.0; // [Ohm] equivalent grid resistance
        double gridInductance = 0.0; // [H] equivalent grid inductance
        double resistance = 0.0;     // [Ohm] load resistance for R‑L loads
        double inductance = 0.0;     // [H] load inductance for R‑L loads
    };

    // Electrical operating point for which ripple will be calculated.
    // The power factor may be specified either directly or via the current
    // phase angle.  When connected to a grid the output power is required.
    struct OperatingPoint {
        double fundamentalFrequency = 0.0;      // [Hz]
        std::optional<double> powerFactor;      // cos(phi)
        std::optional<double> currentPhaseAngle; // [deg]
        std::optional<double> outputPower;      // [W]
        Load load;                              // Load description
    };

private:
    // Basic inverter configuration parameters taken from the JSON input
    double dcBusVoltage = 0.0;                  // [V]
    double vdcRipple = 0.0;                     // [V] DC bus ripple (unused)
    double inverterLegPowerRating = 0.0;        // [W]
    double lineRmsCurrent = 0.0;                // [A]
    double switchingFrequency = 0.0;            // [Hz]
    double riseTime = 0.0;                      // [s] transistor rise time
    double deadtime = 0.0;                      // [s]
    std::string pwmType;                        // carrier waveform type
    std::string modulationStrategy;             // SPWM or SVPWM
    double thirdHarmonicInjectionCoefficient = 0.0;
    double modulationDepth = 0.0;               // modulation index
    std::vector<OperatingPoint> operatingPoints; // list of operating points

public:
    TwoLevelInverter() = default;
    explicit TwoLevelInverter(const nlohmann::json& j);

    // Estimate high‑frequency ripple current of the first filter inductor.
    //  l1               – series inductor parameters
    //  opIndex          – index into the operatingPoints vector
    //  fund_cycles      – number of fundamental cycles to simulate
    //  samples_per_carrier – time resolution of the carrier waveform
    //  f_cut            – frequency above which RMS ripple is accumulated
    double compute_current_ripple(const L1Params& l1,
                                  size_t opIndex,
                                  int fund_cycles = 1,
                                  int samples_per_carrier = 20,
                                  double f_cut = 1000.0) const;
};
} // namespace OpenMagnetics
