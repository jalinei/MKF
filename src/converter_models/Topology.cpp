#include "converter_models/Topology.h"
#include "physical_models/MagnetizingInductance.h"
#include "support/Utils.h"
#include <cfloat>
#include <cmath>
#include <algorithm>
#include <complex>

namespace OpenMagnetics {

    Flyback::Flyback(const json& j) {
        from_json(j, *this);
    }

    AdvancedFlyback::AdvancedFlyback(const json& j) {
        from_json(j, *this);
    }

    double calculate_BMO_duty_cycle(double outputVoltage, double inputVoltage, double turnsRatio){
        return (turnsRatio * outputVoltage) / (inputVoltage + turnsRatio * outputVoltage);
    }

    double calculate_BMO_primary_current_peak(double outputCurrent, double efficiency, double dutyCycle, double turnsRatio){
        return (2 * outputCurrent) / (efficiency * (1 - dutyCycle) * turnsRatio);
    }


    double calculate_QRM_frequency(double magnetizingInductance, double outputPower, double outputVoltage, double minimumInputVoltage, double turnsRatio, 
                                   double diodeVoltageDrop, double efficiency, double drainSourceCapacitance = 100e-12) {

        double dt = std::numbers::pi * sqrt(magnetizingInductance * drainSourceCapacitance);
        double a = pow((outputVoltage + diodeVoltageDrop + 1.0 / turnsRatio * minimumInputVoltage),2);
        double b = efficiency * minimumInputVoltage * minimumInputVoltage * pow((outputVoltage + diodeVoltageDrop),2);
        double c = outputVoltage + diodeVoltageDrop + 1.0 / turnsRatio * minimumInputVoltage;
        double d = sqrt(outputPower/(efficiency * magnetizingInductance));
        double e = minimumInputVoltage * (outputVoltage + diodeVoltageDrop);
        double f = sqrt(4 * dt + (2 * magnetizingInductance * outputPower* a) / b);
        double g = (1.414 * magnetizingInductance * c * d) / e;
        double h = 4.0 / pow((f + g),2);

        return h;
    }


    FlybackModes FlybackOperatingPoint::resolve_mode(std::optional<double> currentRippleRatio) {
        if (get_mode()) {
            return get_mode().value();
        }
        else {
            if (!currentRippleRatio) {
                throw std::runtime_error("Either current ripple ratio or mode is needed for the Flyback OperatingPoint Mode");
            }
            auto mode = currentRippleRatio.value() < 1? FlybackModes::CONTINUOUS_CONDUCTION_MODE : FlybackModes::DISCONTINUOUS_CONDUCTION_MODE;
            return mode;
        }
    }
    double FlybackOperatingPoint::resolve_switching_frequency(double inputVoltage, double diodeVoltageDrop, std::optional<double> inductance, std::optional<std::vector<double>> turnsRatios, double efficiency) {
        if (get_switching_frequency()) {
            return get_switching_frequency().value();
        }
        else {
            if (!get_mode()) {
                throw std::runtime_error("Either switching frequency or mode is needed for the Flyback OperatingPoint");
            }
            auto mode = get_mode().value();
            switch (mode) {
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    throw std::runtime_error("Switching Frequency is needed for CCM");
                }
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    throw std::runtime_error("Switching Frequency is needed for DCM");
                }
                case FlybackModes::QUASI_RESONANT_MODE: {
                    if (!inductance) {
                        throw std::runtime_error("Inductance in missing for switching frequency calculation");
                    }
                    if (!turnsRatios) {
                        throw std::runtime_error("TurnsRatios in missing for switching frequency calculation");
                    }
                    double totalOutputVoltageReflectedPrimaryMinusDiode = 0;

                    for (size_t secondaryIndex = 0; secondaryIndex < get_output_voltages().size(); ++secondaryIndex) {
                        auto outputVoltage = get_output_voltages()[secondaryIndex];
                        auto turnsRatio = turnsRatios.value()[secondaryIndex];
                        totalOutputVoltageReflectedPrimaryMinusDiode += outputVoltage * turnsRatio;
                    }
                    
                    double totalOutputPower = Flyback::get_total_input_power(get_output_currents(), get_output_voltages(), 1, 0);

                    double switchingFrequency = calculate_QRM_frequency(inductance.value(), totalOutputPower, totalOutputVoltageReflectedPrimaryMinusDiode / turnsRatios.value()[0], inputVoltage, turnsRatios.value()[0], diodeVoltageDrop, efficiency);
                    return switchingFrequency;
                }
                case FlybackModes::BOUNDARY_MODE_OPERATION: {
                    if (!inductance) {
                        throw std::runtime_error("Inductance in missing for switching frequency calculation");
                    }
                    if (!turnsRatios) {
                        throw std::runtime_error("TurnsRatios in missing for switching frequency calculation");
                    }
                    double currentPeak = 0;
                    double switchingFrequency = 0;
                    for (size_t secondaryIndex = 0; secondaryIndex < get_output_voltages().size(); ++secondaryIndex) {
                        auto outputCurrent = get_output_currents()[secondaryIndex];
                        auto outputVoltage = get_output_voltages()[secondaryIndex];
                        auto turnsRatio = turnsRatios.value()[secondaryIndex];
                        auto dutyCycleMaximum = calculate_BMO_duty_cycle((outputVoltage + diodeVoltageDrop), outputVoltage, turnsRatio);
                        currentPeak = std::max(currentPeak, calculate_BMO_primary_current_peak(outputCurrent, efficiency, dutyCycleMaximum, turnsRatio)); // hardcoded
                        double tOn = (currentPeak * inductance.value()) / inputVoltage;
                        double tOff = (currentPeak * inductance.value()) / (turnsRatio * outputVoltage);

                        switchingFrequency = std::max(switchingFrequency, 1.0 / (tOn + tOff));
                    }
                    return switchingFrequency;
                }
                default:
                  throw std::runtime_error("Unknown mode in Flyback");
            }
        }
    }

    OperatingPoint Flyback::processOperatingPointsForInputVoltage(double inputVoltage, FlybackOperatingPoint outputOperatingPoint, std::vector<double> turnsRatios, double inductance, std::optional<FlybackModes> customMode, std::optional<double> customDutyCycle, std::optional<double> customDeadTime) {

        OperatingPoint operatingPoint;
        double switchingFrequency = outputOperatingPoint.resolve_switching_frequency(inputVoltage, get_diode_voltage_drop(), inductance, turnsRatios, get_efficiency());

        double deadTime = 0;
        double maximumReflectedOutputVoltage = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
            double outputVoltage = outputOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop();
            maximumReflectedOutputVoltage = std::max(maximumReflectedOutputVoltage, outputVoltage * turnsRatios[secondaryIndex]);
        }

        double primaryVoltavePeaktoPeak = inputVoltage + maximumReflectedOutputVoltage;

        double totalOutputPower = get_total_input_power(outputOperatingPoint.get_output_currents(), outputOperatingPoint.get_output_voltages(), 1, 0);
        double maximumEffectiveLoadCurrent = totalOutputPower / outputOperatingPoint.get_output_voltages()[0];
        double maximumEffectiveLoadCurrentReflected = maximumEffectiveLoadCurrent / turnsRatios[0];
        double totalInputPower = get_total_input_power(outputOperatingPoint.get_output_currents(), outputOperatingPoint.get_output_voltages(), get_efficiency(), 0);
        double averageInputCurrent = totalInputPower / inputVoltage;

        double dutyCycle;
        if (customDutyCycle) {
            dutyCycle = customDutyCycle.value();
        }
        else {
            dutyCycle = averageInputCurrent / (averageInputCurrent + maximumEffectiveLoadCurrentReflected);
        }

        double centerSecondaryCurrentRampLumped = maximumEffectiveLoadCurrent / (1 - dutyCycle);
        double centerPrimaryCurrentRamp = centerSecondaryCurrentRampLumped / turnsRatios[0];


        if (customDeadTime) {
            deadTime = customDeadTime.value();
        }

        if (dutyCycle > 1 ) {
            throw std::runtime_error("dutyCycle cannot be larger than one: " + std::to_string(dutyCycle));
        }

        double primaryCurrentAverage = centerPrimaryCurrentRamp;
        double currentRippleRatio;
        if (std::isnan(get_current_ripple_ratio())) {
            double primaryCurrentPeakToPeak = inputVoltage * dutyCycle / switchingFrequency / inductance;
            currentRippleRatio = primaryCurrentPeakToPeak / centerPrimaryCurrentRamp;
        }
        else {
            currentRippleRatio = get_current_ripple_ratio();
        }
        double primaryCurrentPeakToPeak = centerPrimaryCurrentRamp * currentRippleRatio * 2;
        double primaryCurrentOffset = primaryCurrentAverage - primaryCurrentPeakToPeak / 2;
        primaryCurrentOffset = std::max(0.0, primaryCurrentOffset);

        FlybackModes mode;
        if (customMode) {
            mode = customMode.value();
        }
        else {
            if (primaryCurrentOffset > 0) {
                mode = FlybackModes::CONTINUOUS_CONDUCTION_MODE;
            }
            else {
                mode = FlybackModes::DISCONTINUOUS_CONDUCTION_MODE;
            }
        }

        // Primary
        {
            OperatingPointExcitation excitation;
            Waveform currentWaveform;
            Waveform voltageWaveform;
            Processed currentProcessed;
            Processed voltageProcessed;

            currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_PRIMARY, primaryCurrentPeakToPeak, switchingFrequency, dutyCycle, primaryCurrentOffset, deadTime);
            currentProcessed.set_label(WaveformLabel::FLYBACK_PRIMARY);
            currentProcessed.set_peak_to_peak(primaryCurrentPeakToPeak);
            currentProcessed.set_peak(primaryCurrentOffset + primaryCurrentPeakToPeak / 2);
            currentProcessed.set_duty_cycle(dutyCycle);
            currentProcessed.set_offset(primaryCurrentOffset);
            currentProcessed.set_dead_time(deadTime);

            voltageProcessed.set_peak_to_peak(primaryVoltavePeaktoPeak);
            voltageProcessed.set_peak(inputVoltage);
            voltageProcessed.set_duty_cycle(dutyCycle);
            voltageProcessed.set_offset(0);
            voltageProcessed.set_dead_time(deadTime);
            switch (mode) {
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    voltageProcessed.set_label(WaveformLabel::RECTANGULAR);
                    break;
                }
                case FlybackModes::QUASI_RESONANT_MODE:
                case FlybackModes::BOUNDARY_MODE_OPERATION:
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, primaryVoltavePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    voltageProcessed.set_label(WaveformLabel::RECTANGULAR_WITH_DEADTIME);
                    break;
                }
            }

            excitation.set_frequency(switchingFrequency);
            SignalDescriptor current;
            current.set_waveform(currentWaveform);
            currentProcessed = Inputs::calculate_processed_data(currentWaveform, switchingFrequency, true, currentProcessed);
            auto sampledCurrentWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(currentWaveform, switchingFrequency);
            auto currentHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency);
            current.set_processed(currentProcessed);
            current.set_harmonics(currentHarmonics);
            excitation.set_current(current);
            SignalDescriptor voltage;
            voltage.set_waveform(voltageWaveform);
            voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, switchingFrequency, true, voltageProcessed);
            auto sampledVoltageWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
            auto voltageHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency);
            voltage.set_processed(voltageProcessed);
            voltage.set_harmonics(voltageHarmonics);
            excitation.set_voltage(voltage);
            json isolationSideJson;
            to_json(isolationSideJson, get_isolation_side_from_index(0));
            excitation.set_name(isolationSideJson);
            excitation = Inputs::prune_harmonics(excitation, Defaults().harmonicAmplitudeThreshold, 1);

            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        // Secondaries
        for (size_t secondaryIndex = 0; secondaryIndex < turnsRatios.size(); ++secondaryIndex) {

            OperatingPointExcitation excitation;
            Waveform currentWaveform;
            Waveform voltageWaveform;
            Processed currentProcessed;
            Processed voltageProcessed;

            double secondaryPower = get_total_input_power(outputOperatingPoint.get_output_currents()[secondaryIndex], outputOperatingPoint.get_output_voltages()[secondaryIndex], 1, 0);
            double powerDivider = secondaryPower / totalOutputPower;

            double secondaryVoltagePeaktoPeak = inputVoltage / turnsRatios[secondaryIndex] + get_diode_voltage_drop() + outputOperatingPoint.get_output_voltages()[secondaryIndex];
            double secondaryCurrentAverage = centerPrimaryCurrentRamp * turnsRatios[secondaryIndex] * powerDivider;
            double secondaryCurrentPeaktoPeak = secondaryCurrentAverage * currentRippleRatio * 2;
            double secondaryCurrentOffset = std::max(0.0, secondaryCurrentAverage - secondaryCurrentPeaktoPeak / 2);


            currentProcessed.set_peak_to_peak(secondaryCurrentPeaktoPeak);
            currentProcessed.set_peak(secondaryCurrentOffset + secondaryCurrentPeaktoPeak / 2);
            currentProcessed.set_duty_cycle(dutyCycle);
            currentProcessed.set_offset(secondaryCurrentOffset);
            currentProcessed.set_dead_time(deadTime);

            voltageProcessed.set_peak_to_peak(secondaryVoltagePeaktoPeak);
            voltageProcessed.set_peak(outputOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop());
            voltageProcessed.set_duty_cycle(dutyCycle);
            voltageProcessed.set_offset(0);
            voltageProcessed.set_dead_time(deadTime);

            switch (mode) {
                case FlybackModes::CONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR, secondaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_SECONDARY, secondaryCurrentPeaktoPeak, switchingFrequency, dutyCycle, secondaryCurrentOffset, deadTime);
                    voltageProcessed.set_label(WaveformLabel::SECONDARY_RECTANGULAR);
                    currentProcessed.set_label(WaveformLabel::FLYBACK_SECONDARY);
                    break;
                }
                case FlybackModes::QUASI_RESONANT_MODE:
                case FlybackModes::BOUNDARY_MODE_OPERATION:
                case FlybackModes::DISCONTINUOUS_CONDUCTION_MODE: {
                    voltageWaveform = Inputs::create_waveform(WaveformLabel::RECTANGULAR_WITH_DEADTIME, secondaryVoltagePeaktoPeak, switchingFrequency, dutyCycle, 0, deadTime);
                    currentWaveform = Inputs::create_waveform(WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME, secondaryCurrentPeaktoPeak, switchingFrequency, dutyCycle, secondaryCurrentOffset, deadTime);
                    voltageProcessed.set_label(WaveformLabel::SECONDARY_RECTANGULAR_WITH_DEADTIME);
                    currentProcessed.set_label(WaveformLabel::FLYBACK_SECONDARY_WITH_DEADTIME);
                    break;
                }
            }

            excitation.set_frequency(switchingFrequency);
            SignalDescriptor current;
            current.set_waveform(currentWaveform);
            currentProcessed = Inputs::calculate_processed_data(currentWaveform, switchingFrequency, true, currentProcessed);
            auto sampledCurrentWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(currentWaveform, switchingFrequency);
            auto currentHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledCurrentWaveform, switchingFrequency);
            current.set_processed(currentProcessed);
            current.set_harmonics(currentHarmonics);
            excitation.set_current(current);
            SignalDescriptor voltage;
            voltage.set_waveform(voltageWaveform);
            voltageProcessed = Inputs::calculate_processed_data(voltageWaveform, switchingFrequency, true, voltageProcessed);
            auto sampledVoltageWaveform = OpenMagnetics::Inputs::calculate_sampled_waveform(voltageWaveform, switchingFrequency);
            auto voltageHarmonics = OpenMagnetics::Inputs::calculate_harmonics_data(sampledVoltageWaveform, switchingFrequency);
            voltage.set_processed(voltageProcessed);
            voltage.set_harmonics(voltageHarmonics);
            excitation.set_voltage(voltage);
            json isolationSideJson;
            to_json(isolationSideJson, get_isolation_side_from_index(secondaryIndex + 1));
            excitation.set_name(isolationSideJson);
            excitation = Inputs::prune_harmonics(excitation, Defaults().harmonicAmplitudeThreshold, 1);
            operatingPoint.get_mutable_excitations_per_winding().push_back(excitation);
        }

        OperatingConditions conditions;
        conditions.set_ambient_temperature(outputOperatingPoint.get_ambient_temperature());
        conditions.set_cooling(std::nullopt);
        operatingPoint.set_conditions(conditions);

        return operatingPoint;
    }

    double Flyback::get_total_input_power(std::vector<double> outputCurrents, std::vector<double> outputVoltages, double efficiency, double diodeVoltageDrop) {
        double totalPower = 0;
        for (size_t secondaryIndex = 0; secondaryIndex < outputCurrents.size(); ++secondaryIndex) {
            totalPower += outputCurrents[secondaryIndex] * (outputVoltages[secondaryIndex] + diodeVoltageDrop);
        }

        return totalPower / efficiency;
    }


    double Flyback::get_total_input_power(double outputCurrent, double outputVoltage, double efficiency, double diodeVoltageDrop) {
        double totalPower = outputCurrent * (outputVoltage + diodeVoltageDrop);

        return totalPower / efficiency;
    }

    double Flyback::get_minimum_output_reflected_voltage(double maximumDrainSourceVoltage, double maximumInputVoltage, double safetyMargin) {
        return maximumDrainSourceVoltage * safetyMargin - maximumInputVoltage;
    }

    bool Flyback::run_checks(bool assert) {
        if (get_operating_points().size() == 0) {
            if (!assert) {
                return false;
            }
            throw std::runtime_error("At least one operating point is needed");
        }
        for (size_t flybackOperatingPointIndex = 1; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
            if (get_operating_points()[flybackOperatingPointIndex].get_output_voltages().size() != get_operating_points()[0].get_output_voltages().size()) {
                if (!assert) {
                    return false;
                }
                throw std::runtime_error("Different operating points cannot have different number of output voltages");
            }
            if (get_operating_points()[flybackOperatingPointIndex].get_output_currents().size() != get_operating_points()[0].get_output_currents().size()) {
                if (!assert) {
                    return false;
                }
                throw std::runtime_error("Different operating points cannot have different number of output currents");
            }
        }
        if (!get_input_voltage().get_nominal() && !get_input_voltage().get_maximum() && !get_input_voltage().get_minimum()) {
            if (!assert) {
                return false;
            }
            throw std::runtime_error("No input voltage introduced");
        }

        return true;
    }

    DesignRequirements Flyback::process_design_requirements() {
        double minimumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MINIMUM);
        double maximumInputVoltage = resolve_dimensional_values(get_input_voltage(), DimensionalValues::MAXIMUM);

        if (!get_maximum_drain_source_voltage() && !get_maximum_duty_cycle()) {
            throw std::invalid_argument("Missing both maximum duty cycle and maximum drain source voltage");
        }
        double maximumNeededInductance = 0;
        std::vector<double> turnsRatios(get_operating_points()[0].get_output_voltages().size(), 0);

        if (get_maximum_duty_cycle()) {
            double maximumDutyCycle = get_maximum_duty_cycle().value();
            if (maximumDutyCycle > 1 || maximumDutyCycle < 0) {
                throw std::invalid_argument("maximumDutyCycle must be between 0 and 1");
            }
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];

                double totalOutputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), 1, 0);
                double totalInputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), get_efficiency(), 0);
                double maximumEffectiveLoadCurrent = totalOutputPower / flybackOperatingPoint.get_output_voltages()[0];
                double averageInputCurrent = totalInputPower / minimumInputVoltage;
                double maximumEffectiveLoadCurrentReflected = averageInputCurrent * (1 - maximumDutyCycle) / maximumDutyCycle;

                auto turnsRatioFirstOutput = maximumEffectiveLoadCurrent / maximumEffectiveLoadCurrentReflected;
                turnsRatios[0] = std::max(turnsRatios[0], turnsRatioFirstOutput);

                for (size_t secondaryIndex = 1; secondaryIndex < flybackOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                    auto turnsRatio = turnsRatioFirstOutput * (flybackOperatingPoint.get_output_voltages()[0] + get_diode_voltage_drop()) / (flybackOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop());
                    turnsRatios[secondaryIndex] = std::max(turnsRatios[secondaryIndex], turnsRatio);
                }
            }
        }

        if (get_maximum_drain_source_voltage()) {
            std::vector<double> turnsRatiosFromMaximumDrainSourceVoltage(get_operating_points()[0].get_output_voltages().size(), 0);
            double maximumDrainSourceVoltage = get_maximum_drain_source_voltage().value();
            auto minimumOutputReflectedVoltage = get_minimum_output_reflected_voltage(maximumDrainSourceVoltage, maximumInputVoltage);
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];
                for (size_t secondaryIndex = 0; secondaryIndex < flybackOperatingPoint.get_output_voltages().size(); ++secondaryIndex) {
                    auto turnsRatio = minimumOutputReflectedVoltage / (flybackOperatingPoint.get_output_voltages()[secondaryIndex] + get_diode_voltage_drop());
                    turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex] = std::max(turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex], turnsRatio);
                }
            }

            for (size_t secondaryIndex = 0; secondaryIndex < get_operating_points()[0].get_output_voltages().size(); ++secondaryIndex) {
                if (turnsRatios[secondaryIndex] > 1) {
                    turnsRatios[secondaryIndex] = std::min(turnsRatios[secondaryIndex], turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex]);
                }
                else {
                    turnsRatios[secondaryIndex] = std::max(turnsRatios[secondaryIndex], turnsRatiosFromMaximumDrainSourceVoltage[secondaryIndex]);
                }
            }

        }

        for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
            auto flybackOperatingPoint = get_operating_points()[flybackOperatingPointIndex];
            double switchingFrequency = flybackOperatingPoint.resolve_switching_frequency(minimumInputVoltage, get_diode_voltage_drop());
            double totalOutputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), 1, 0);
            double maximumEffectiveLoadCurrent = totalOutputPower / flybackOperatingPoint.get_output_voltages()[0];
            double dutyCycle = 0;
            if (get_maximum_duty_cycle()) {
                dutyCycle = get_maximum_duty_cycle().value();
            }
            else {
                double maximumEffectiveLoadCurrentReflected = maximumEffectiveLoadCurrent / turnsRatios[0];
                double totalInputPower = get_total_input_power(flybackOperatingPoint.get_output_currents(), flybackOperatingPoint.get_output_voltages(), get_efficiency(), 0);
                double averageInputCurrent = totalInputPower / minimumInputVoltage;
                dutyCycle = averageInputCurrent / (averageInputCurrent + maximumEffectiveLoadCurrentReflected);
            }

            double centerSecondaryCurrentRampLumped = maximumEffectiveLoadCurrent / (1 - dutyCycle);
            double centerPrimaryCurrentRamp = centerSecondaryCurrentRampLumped / turnsRatios[0];
            double tOn = dutyCycle / switchingFrequency;
            double voltsSeconds = minimumInputVoltage * tOn;
            double neededInductance = voltsSeconds / get_current_ripple_ratio() / centerPrimaryCurrentRamp;
            maximumNeededInductance = std::max(maximumNeededInductance, neededInductance);
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        for (auto turnsRatio : turnsRatios) {
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(roundFloat(turnsRatio, 2));
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_minimum(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        for (size_t windingIndex = 0; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::FLYBACK_CONVERTER);
        return designRequirements;
    }

    std::vector<OperatingPoint> Flyback::process_operating_points(std::vector<double> turnsRatios, double magnetizingInductance) {
        std::vector<OperatingPoint> operatingPoints;
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;

        if (get_input_voltage().get_nominal()) {
            inputVoltages.push_back(get_input_voltage().get_nominal().value());
            inputVoltagesNames.push_back("Nom.");
        }
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
            inputVoltagesNames.push_back("Min.");
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
            inputVoltagesNames.push_back("Max.");
        }

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                auto mode = get_mutable_operating_points()[flybackOperatingPointIndex].resolve_mode(get_current_ripple_ratio());
                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, magnetizingInductance, mode);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(flybackOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }
        return operatingPoints;
    }

    Inputs Flyback::process() {
        Flyback::run_checks(_assertErrors);

        Inputs inputs;
        auto designRequirements = process_design_requirements();
        std::vector<double> turnsRatios;
        for (auto turnsRatio : designRequirements.get_turns_ratios()) {
            turnsRatios.push_back(resolve_dimensional_values(turnsRatio));
        }
        auto desiredMagnetizingInductance = resolve_dimensional_values(designRequirements.get_magnetizing_inductance());
        auto operatingPoints = process_operating_points(turnsRatios, desiredMagnetizingInductance);

        inputs.set_design_requirements(designRequirements);
        inputs.set_operating_points(operatingPoints);

        return inputs;
    }

    std::vector<OperatingPoint> Flyback::process_operating_points(Magnetic magnetic) {
        Flyback::run_checks(_assertErrors);

        std::vector<OperatingPoint> operatingPoints;

        if (!get_maximum_drain_source_voltage() && !get_maximum_duty_cycle()) {
            throw std::invalid_argument("Missing both maximum duty cycle and maximum drain source voltage");
        }
        OpenMagnetics::MagnetizingInductance magnetizingInductanceModel("ZHANG");  // hardcoded
        double magnetizingInductance = magnetizingInductanceModel.calculate_inductance_from_number_turns_and_gapping(magnetic.get_mutable_core(), magnetic.get_mutable_coil()).get_magnetizing_inductance().get_nominal().value();
        std::vector<double> turnsRatios = magnetic.get_turns_ratios();

        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;


        if (get_input_voltage().get_nominal()) {
            inputVoltages.push_back(get_input_voltage().get_nominal().value());
            inputVoltagesNames.push_back("Nom.");
        }
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
            inputVoltagesNames.push_back("Min.");
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
            inputVoltagesNames.push_back("Max.");
        }

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                auto mode = get_mutable_operating_points()[flybackOperatingPointIndex].resolve_mode(get_current_ripple_ratio());
                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, magnetizingInductance, mode);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(flybackOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                operatingPoints.push_back(operatingPoint);
            }
        }

        return operatingPoints;
    }

    Inputs AdvancedFlyback::process() {
        Flyback::run_checks(_assertErrors);

        Inputs inputs;

        double maximumNeededInductance = get_desired_inductance();
        std::vector<double> turnsRatios = get_desired_turns_ratios();

        inputs.get_mutable_operating_points().clear();
        std::vector<double> inputVoltages;
        std::vector<std::string> inputVoltagesNames;


        if (get_input_voltage().get_nominal()) {
            inputVoltages.push_back(get_input_voltage().get_nominal().value());
            inputVoltagesNames.push_back("Nom.");
        }
        if (get_input_voltage().get_maximum()) {
            inputVoltages.push_back(get_input_voltage().get_maximum().value());
            inputVoltagesNames.push_back("Max.");
        }
        if (get_input_voltage().get_minimum()) {
            inputVoltages.push_back(get_input_voltage().get_minimum().value());
            inputVoltagesNames.push_back("Min.");
        }

        DesignRequirements designRequirements;
        designRequirements.get_mutable_turns_ratios().clear();
        for (auto turnsRatio : turnsRatios) {
            DimensionWithTolerance turnsRatioWithTolerance;
            turnsRatioWithTolerance.set_nominal(roundFloat(turnsRatio, 2));
            designRequirements.get_mutable_turns_ratios().push_back(turnsRatioWithTolerance);
        }
        DimensionWithTolerance inductanceWithTolerance;
        inductanceWithTolerance.set_nominal(roundFloat(maximumNeededInductance, 10));
        designRequirements.set_magnetizing_inductance(inductanceWithTolerance);
        std::vector<IsolationSide> isolationSides;
        for (size_t windingIndex = 0; windingIndex < turnsRatios.size() + 1; ++windingIndex) {
            isolationSides.push_back(get_isolation_side_from_index(windingIndex));
        }
        designRequirements.set_isolation_sides(isolationSides);
        designRequirements.set_topology(Topologies::FLYBACK_CONVERTER);

        inputs.set_design_requirements(designRequirements);

        for (size_t inputVoltageIndex = 0; inputVoltageIndex < inputVoltages.size(); ++inputVoltageIndex) {
            auto inputVoltage = inputVoltages[inputVoltageIndex];
            for (size_t flybackOperatingPointIndex = 0; flybackOperatingPointIndex < get_operating_points().size(); ++flybackOperatingPointIndex) {
                std::optional<double> customDeadTime = std::nullopt;
                if (get_desired_duty_cycle().size() <= flybackOperatingPointIndex) {
                    throw std::runtime_error("Missing duty cycle for flybackOperatingPointIndex: " + std::to_string(flybackOperatingPointIndex));
                }
                double customDutyCycle = get_desired_duty_cycle()[flybackOperatingPointIndex][inputVoltageIndex];

                if (get_desired_dead_time()) {
                    if (get_desired_dead_time()->size() <= flybackOperatingPointIndex) {
                        throw std::runtime_error("Missing dead time for flybackOperatingPointIndex: " + std::to_string(flybackOperatingPointIndex));
                    }
                    customDeadTime = get_desired_dead_time().value()[flybackOperatingPointIndex];
                }

                auto operatingPoint = processOperatingPointsForInputVoltage(inputVoltage, get_operating_points()[flybackOperatingPointIndex], turnsRatios, maximumNeededInductance, std::nullopt, customDutyCycle, customDeadTime);

                std::string name = inputVoltagesNames[inputVoltageIndex] + " input volt.";
                if (get_operating_points().size() > 1) {
                    name += " with op. point " + std::to_string(flybackOperatingPointIndex);
                }
                operatingPoint.set_name(name);
                inputs.get_mutable_operating_points().push_back(operatingPoint);
            }
        }

        return inputs;
    }

} // namespace OpenMagnetics
namespace OpenMagnetics {

using json = nlohmann::json;

struct LegReference { double Vg_phase_rms; double I_leg_rms; double I_leg_angle_deg; double P_leg_w; double Q_leg_var; std::complex<double> Zg_ohm; std::complex<double> V_pcc_ref; std::complex<double> V_inv_ref; double V_pcc_rms; double V_pcc_angle_deg; double V_inv_rms; double V_inv_angle_deg; };

static std::vector<int> delay_rising_edges(const std::vector<int>& bit, int Ndt) {
    if (Ndt <= 0) return bit;
    std::vector<int> out(bit);
    for (size_t i = 1; i < bit.size(); ++i) {
        if (bit[i] == 1 && bit[i-1] == 0) {
            size_t k2 = std::min(static_cast<size_t>(i + Ndt), out.size());
            for (size_t k = i; k < k2; ++k) out[k] = 0;
        }
    }
    return out;
}

static void apply_deadtime_to_complementary(const std::vector<int>& s_top,
                                            double t_dead_s,
                                            double fs,
                                            std::vector<int>& g_top,
                                            std::vector<int>& g_bot,
                                            int& Ndt) {
    Ndt = static_cast<int>(std::round(t_dead_s * fs));
    g_top = delay_rising_edges(s_top, Ndt);
    std::vector<int> s_bot_ideal(s_top.size());
    for (size_t i = 0; i < s_top.size(); ++i) s_bot_ideal[i] = 1 - s_top[i];
    g_bot = delay_rising_edges(s_bot_ideal, Ndt);
    for (size_t i = 0; i < g_top.size(); ++i) {
        if (g_top[i] == 1 && g_bot[i] == 1) {
            g_top[i] = 0; g_bot[i] = 0;
        }
    }
}

static std::vector<double> pole_voltage_with_deadtime(const std::vector<int>& g_top,
                                                      const std::vector<int>& g_bot,
                                                      double Vdc,
                                                      const std::vector<double>& i_leg_inst) {
    std::vector<double> v(i_leg_inst.size());
    for (size_t i = 0; i < v.size(); ++i) {
        bool pos = g_top[i] == 1;
        bool neg = g_bot[i] == 1;
        bool off = !(pos || neg);
        if (pos) v[i] = +0.5 * Vdc;
        else if (neg) v[i] = -0.5 * Vdc;
        else if (off && i_leg_inst[i] >= 0) v[i] = -0.5 * Vdc;
        else v[i] = +0.5 * Vdc;
    }
    return v;
}

static std::vector<double> fundamental_leg_current_wave(double I_leg_rms,
                                                        double I_leg_angle_deg,
                                                        double f1,
                                                        const std::vector<double>& t) {
    std::vector<double> i(t.size());
    double Ipk = std::sqrt(2.0) * I_leg_rms;
    double angle = I_leg_angle_deg * M_PI / 180.0;
    for (size_t k = 0; k < t.size(); ++k) {
        i[k] = Ipk * std::sin(2.0 * M_PI * f1 * t[k] + angle);
    }
    return i;
}

static LegReference compute_leg_reference(double P_leg_w,
                                          double PF,
                                          bool lagging,
                                          const GridParams& grid,
                                          const L1Params& l1) {
    LegReference res;
    double V_phase_rms = grid.VLL_rms / std::sqrt(3.0);
    double omega = 2.0 * M_PI * grid.f_hz;
    std::complex<double> Zg(grid.Rgrid_ohm, omega * grid.Lgrid_h);
    double I_leg_rms = P_leg_w / (V_phase_rms * PF);
    double phi = std::acos(PF);
    double I_angle = lagging ? -phi : +phi;
    std::complex<double> I_ph = std::polar(I_leg_rms, I_angle);
    double Q_leg = V_phase_rms * I_leg_rms * std::sin(phi);
    if (!lagging) Q_leg = -Q_leg;
    std::complex<double> Vg_ph(V_phase_rms, 0);
    std::complex<double> V_pcc = Vg_ph - I_ph * Zg;
    std::complex<double> Z1(l1.R1_ohm, omega * l1.L1_h);
    std::complex<double> V_inv = V_pcc + I_ph * Z1;

    res.Vg_phase_rms = V_phase_rms;
    res.I_leg_rms = I_leg_rms;
    res.I_leg_angle_deg = I_angle * 180.0 / M_PI;
    res.P_leg_w = P_leg_w;
    res.Q_leg_var = Q_leg;
    res.Zg_ohm = Zg;
    res.V_pcc_ref = V_pcc;
    res.V_inv_ref = V_inv;
    res.V_pcc_rms = std::abs(V_pcc);
    res.V_pcc_angle_deg = std::arg(V_pcc) * 180.0 / M_PI;
    res.V_inv_rms = std::abs(V_inv);
    res.V_inv_angle_deg = std::arg(V_inv) * 180.0 / M_PI;
    return res;
}

static std::vector<double> tri_wave(const std::vector<double>& t, double f_carrier) {
    std::vector<double> tri(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        double x = std::fmod(t[i] * f_carrier, 1.0);
        tri[i] = 4.0 * std::abs(x - 0.5) - 1.0;
    }
    return tri;
}

static void three_phase_references(double V_phase_rms,
                                   double f1,
                                   double phase_deg,
                                   const std::vector<double>& t,
                                   double third_harm_coeff,
                                   std::vector<double>& va,
                                   std::vector<double>& vb,
                                   std::vector<double>& vc) {
    double A = std::sqrt(2.0) * V_phase_rms;
    double w1 = 2.0 * M_PI * f1;
    double phase = phase_deg * M_PI / 180.0;
    va.resize(t.size()); vb.resize(t.size()); vc.resize(t.size());
    for (size_t i = 0; i < t.size(); ++i) {
        double th = w1 * t[i] + phase;
        double baseA = std::sin(th);
        double baseB = std::sin(th - 2.0 * M_PI / 3.0);
        double baseC = std::sin(th + 2.0 * M_PI / 3.0);
        double v3 = (third_harm_coeff != 0.0) ? std::sin(3.0 * th) * third_harm_coeff : 0.0;
        va[i] = A * (baseA + v3);
        vb[i] = A * (baseB + v3);
        vc[i] = A * (baseC + v3);
    }
}

static void svpwm_zero_sequence(std::vector<double>& va,
                                std::vector<double>& vb,
                                std::vector<double>& vc) {
    for (size_t i = 0; i < va.size(); ++i) {
        double vmax = std::max({va[i], vb[i], vc[i]});
        double vmin = std::min({va[i], vb[i], vc[i]});
        double e0 = -0.5 * (vmax + vmin);
        va[i] += e0; vb[i] += e0; vc[i] += e0;
    }
}

static std::vector<int> carrier_compare(const std::vector<double>& v_norm,
                                        const std::vector<double>& tri) {
    std::vector<int> s(v_norm.size());
    for (size_t i = 0; i < v_norm.size(); ++i) {
        double v = std::clamp(v_norm[i], -1.0, 1.0);
        s[i] = (v >= tri[i]) ? 1 : 0;
    }
    return s;
}

static std::pair<std::vector<double>, std::vector<double>> rms_spectrum(const std::vector<double>& x, double fs) {
    size_t N = x.size();
    size_t K = N / 2 + 1;
    std::vector<std::complex<double>> X(K, {0.0, 0.0});
    for (size_t k = 0; k < K; ++k) {
        std::complex<double> sum(0.0, 0.0);
        for (size_t n = 0; n < N; ++n) {
            double ang = -2.0 * M_PI * k * n / static_cast<double>(N);
            sum += std::complex<double>(x[n] * std::cos(ang), x[n] * std::sin(ang));
        }
        X[k] = sum / static_cast<double>(N);
    }
    std::vector<double> mag(K), Vrms(K, 0.0), freqs(K);
    for (size_t k = 0; k < K; ++k) {
        mag[k] = std::abs(X[k]);
        freqs[k] = k * fs / static_cast<double>(N);
    }
    Vrms[0] = mag[0];
    if (N % 2 == 0) {
        for (size_t k = 1; k < K - 1; ++k) Vrms[k] = (2.0 * mag[k]) / std::sqrt(2.0);
        Vrms[K - 1] = mag[K - 1] / std::sqrt(2.0);
    } else {
        for (size_t k = 1; k < K; ++k) Vrms[k] = (2.0 * mag[k]) / std::sqrt(2.0);
    }
    return {freqs, Vrms};
}

static std::vector<double> denom_mag(const std::vector<double>& freqs,
                                     const L1Params& l1,
                                     const GridParams& grid) {
    std::vector<double> denom(freqs.size());
    for (size_t i = 0; i < freqs.size(); ++i) {
        double w = 2.0 * M_PI * freqs[i];
        std::complex<double> Z1(l1.R1_ohm, w * l1.L1_h);
        std::complex<double> Zg(grid.Rgrid_ohm, w * grid.Lgrid_h);
        denom[i] = std::abs(Z1 + Zg);
    }
    return denom;
}

TwoLevelInverter::TwoLevelInverter(const json& j) {
    dcBusVoltage = resolve_dimensional_values(j.at("dcBusVoltage"));
    vdcRipple = resolve_dimensional_values(j.at("vdcRipple"));
    inverterLegPowerRating = j.at("inverterLegPowerRating").get<double>();
    lineRmsCurrent = resolve_dimensional_values(j.at("lineRmsCurrent"));
    switchingFrequency = j.at("switchingFrequency").get<double>();
    riseTime = j.at("riseTime").get<double>();
    deadtime = j.at("deadtime").get<double>();
    pwmType = j.at("pwmType").get<std::string>();
    modulationStrategy = j.at("modulationStrategy").get<std::string>();
    thirdHarmonicInjectionCoefficient = j.value("thirdHarmonicInjectionCoefficient", 0.0);
    modulationDepth = j.at("modulationDepth").get<double>();

    for (const auto& opj : j.at("operatingPoints")) {
        OperatingPoint op;
        op.fundamentalFrequency = opj.at("fundamentalFrequency").get<double>();
        if (opj.contains("powerFactor")) op.powerFactor = opj.at("powerFactor").get<double>();
        if (opj.contains("currentPhaseAngle")) op.currentPhaseAngle = opj.at("currentPhaseAngle").get<double>();
        if (opj.contains("outputPower")) op.outputPower = opj.at("outputPower").get<double>();
        auto loadj = opj.at("load");
        op.load.type = loadj.at("type").get<std::string>();
        if (op.load.type == "grid") {
            op.load.phaseVoltage = resolve_dimensional_values(loadj.at("phaseVoltage"));
            op.load.gridFrequency = loadj.at("gridFrequency").get<double>();
            op.load.gridResistance = resolve_dimensional_values(loadj.at("gridResistance"));
            op.load.gridInductance = resolve_dimensional_values(loadj.at("gridInductance"));
        } else {
            op.load.resistance = resolve_dimensional_values(loadj.at("resistance"));
            op.load.inductance = resolve_dimensional_values(loadj.at("inductance"));
        }
        operatingPoints.push_back(op);
    }
}

double TwoLevelInverter::compute_current_ripple(const L1Params& l1,
                                                size_t opIndex,
                                                int fund_cycles,
                                                int samples_per_carrier,
                                                double f_cut) const {
    const auto& op = operatingPoints.at(opIndex);
    GridParams grid{op.load.phaseVoltage * std::sqrt(3.0), op.load.gridFrequency,
                    op.load.gridResistance, op.load.gridInductance};

    double PF;
    bool lagging = true;
    if (op.powerFactor.has_value()) {
        PF = op.powerFactor.value();
    } else if (op.currentPhaseAngle.has_value()) {
        double ang = op.currentPhaseAngle.value() * M_PI / 180.0;
        PF = std::cos(ang);
        lagging = ang < 0;
    } else {
        PF = 1.0;
    }
    double P_leg = op.outputPower.value_or(inverterLegPowerRating);
    auto ref = compute_leg_reference(P_leg, PF, lagging, grid, l1);

    double f1 = op.fundamentalFrequency;
    double f_carrier = switchingFrequency;
    double fsamp = f_carrier * samples_per_carrier;
    double T = static_cast<double>(fund_cycles) / f1;
    size_t N = static_cast<size_t>(std::round(T * fsamp));
    T = N / fsamp;
    std::vector<double> t(N);
    for (size_t i = 0; i < N; ++i) t[i] = i / fsamp;

    auto i_leg = fundamental_leg_current_wave(ref.I_leg_rms, ref.I_leg_angle_deg, f1, t);
    auto tri = tri_wave(t, f_carrier);
    std::vector<double> va, vb, vc;
    three_phase_references(ref.V_inv_rms, f1, ref.V_inv_angle_deg, t,
                           thirdHarmonicInjectionCoefficient, va, vb, vc);
    if (modulationStrategy == "SVPWM") {
        svpwm_zero_sequence(va, vb, vc);
    }
    std::vector<int> s_a = carrier_compare(va, tri);
    std::vector<int> g_top, g_bot; int Ndt;
    apply_deadtime_to_complementary(s_a, deadtime, fsamp, g_top, g_bot, Ndt);
    auto v_pole = pole_voltage_with_deadtime(g_top, g_bot, dcBusVoltage, i_leg);
    auto [freqs, Vrms] = rms_spectrum(v_pole, fsamp);
    auto denom = denom_mag(freqs, l1, grid);
    std::vector<double> Irms(freqs.size(), 0.0);
    for (size_t i = 0; i < freqs.size(); ++i) {
        if (denom[i] > 0) Irms[i] = Vrms[i] / denom[i];
    }
    double sum = 0.0;
    for (size_t i = 0; i < freqs.size(); ++i) {
        if (freqs[i] >= f_cut) sum += Irms[i] * Irms[i];
    }
    return std::sqrt(sum);
}

} // namespace OpenMagnetics
