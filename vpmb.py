#!/usr/bin/python

# Copyright 2010, Bryan Waite, Erik C. Baker. All rights reserved.
# Redistribution and use in source and binary forms, with or without modification, are
# permitted provided that the following conditions are met:

#    1. Redistributions of source code must retain the above copyright notice, this list of
#       conditions and the following disclaimer.

#    2. Redistributions in binary form must reproduce the above copyright notice, this list
#       of conditions and the following disclaimer in the documentation and/or other materials
#       provided with the distribution.

# THIS SOFTWARE IS PROVIDED BY Bryan Waite, Erik C. Baker ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
# FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Bryan Waite, or Erik C. Baker OR
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
# ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

# The views and conclusions contained in the software and documentation are those of the
# authors and should not be interpreted as representing official policies, either expressed
# or implied, of Bryan Waite, or Erik C. Baker.

# import pycallgraph

import json
from math import log, trunc, exp, sqrt
import datetime
from optparse import OptionParser


class AltitudeException(Exception):
    """Thrown when altitude is invalid, or diver acclimatized is invalid."""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class MaxIterationException(Exception):
    """Thrown when root finding fails."""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class InputFileException(Exception):
    """Thrown when there are errors with the input file values."""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class DecompressionStepException(Exception):
    """Thrown when the decompression step is too large."""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class RootException(Exception):
    """Thrown when root calculated is not within brackets"""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class OffGassingException(Exception):
    """Thrown when Off Gassing gradient is too small"""
    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)


class DiveState(object):
    """Contains the program state so that this isn't a huge mess"""

    def __init__(self, input_file_name=None, json_input=None):
        """Take the input_file_name  parsed by parse_settings or raw json data
        and use it to initialize the program state"""

        if input_file_name is None and json_input is None:
            raise ValueError("""DiveState must be given a file name (to load the input from), or
            raw json data""")

        data = None
        if json_input:
            data = json_input
        else:
            # load the file data
            input_file = open(input_file_name)
            data = json.loads(input_file.read())
            input_file.close()

        self.input_values = data["input"]
        self.settings_values = data["settings"]
        self.altitude_values = data["altitude"]
        self.output_object = Output(self)

        # init the instance variables
        # strings
        # self.Word = ""
        self.Units = ""
        self.Units_Word1 = ""
        self.Units_Word2 = ""

        # integers
        # self.Number_of_Mixes = 0
        self.Number_of_Changes = 0
        self.Segment_Number_Start_of_Ascent = 0
        self.Repetitive_Dive_Flag = 0

        # bools
        self.Schedule_Converged = False
        self.Critical_Volume_Algorithm_Off = False
        self.Altitude_Dive_Algorithm_Off = False

        # floats
        self.Ascent_Ceiling_Depth = 0.0
        self.Deco_Stop_Depth = 0.0
        self.Step_Size = 0.0
        self.Sum_Check = 0.0
        self.Depth = 0.0
        self.Ending_Depth = 0.0
        self.Starting_Depth = 0.0
        self.Rate = 0.0
        self.Run_Time_End_of_Segment = 0.0
        self.Last_Run_Time = 0.0
        self.Stop_Time = 0.0
        self.Depth_Start_of_Deco_Zone = 0.0
        self.Deepest_Possible_Stop_Depth = 0.0
        self.First_Stop_Depth = 0.0
        self.Critical_Volume_Comparison = 0.0
        self.Next_Stop = 0.0
        self.Run_Time_Start_of_Deco_Zone = 0.0
        self.Critical_Radius_N2_Microns = 0.0
        self.Critical_Radius_He_Microns = 0.0
        self.Run_Time_Start_of_Ascent = 0.0
        self.Altitude_of_Dive = 0.0
        self.Deco_Phase_Volume_Time = 0.0
        self.Surface_Interval_Time = 0.0
        self.Regenerated_Radius_He = [0.0 for i in range(16)]
        self.Regenerated_Radius_N2 = [0.0 for i in range(16)]

        # Global Arrays
        self.Mix_Change = []
        self.Depth_Change = []
        self.Rate_Change = []
        self.Step_Size_Change = []
        self.He_Pressure_Start_of_Ascent = [0.0 for i in range(16)]
        self.N2_Pressure_Start_of_Ascent = [0.0 for i in range(16)]
        self.He_Pressure_Start_of_Deco_Zone = [0.0 for i in range(16)]
        self.N2_Pressure_Start_of_Deco_Zone = [0.0 for i in range(16)]
        self.Phase_Volume_Time = [0.0 for i in range(16)]
        self.Last_Phase_Volume_Time = [0.0 for i in range(16)]

        self.Allowable_Gradient_He = [0.0 for i in range(16)]
        self.Allowable_Gradient_N2 = [0.0 for i in range(16)]

        self.Adjusted_Crushing_Pressure_He = [0.0 for i in range(16)]
        self.Adjusted_Crushing_Pressure_N2 = [0.0 for i in range(16)]

        self.Initial_Allowable_Gradient_N2 = [0.0 for i in range(16)]
        self.Initial_Allowable_Gradient_He = [0.0 for i in range(16)]

        self.Deco_Gradient_He = [0.0 for i in range(16)]
        self.Deco_Gradient_N2 = [0.0 for i in range(16)]

        # GLOBAL CONSTANTS

        self.Water_Vapor_Pressure = 0.0
        self.Surface_Tension_Gamma = 0.0
        self.Skin_Compression_GammaC = 0.0
        self.Crit_Volume_Parameter_Lambda = 0.0
        self.Minimum_Deco_Stop_Time = 0.0
        self.Regeneration_Time_Constant = 0.0
        self.Constant_Pressure_Other_Gases = 0.0
        self.Gradient_Onset_of_Imperm_Atm = 0.0

        self.ATM = 101325.0  # 1 atm of pressure
        self.fraction_inert_gas = 0.79  # oxygen = 21% of air so this is what's left over

        # GLOBAL VARIABLES

        self.Segment_Number = 0
        self.Run_Time = 0.0
        self.Segment_Time = 0.0
        self.Ending_Ambient_Pressure = 0.0
        self.Mix_Number = 0
        self.Barometric_Pressure = 0.0

        self.units_fsw = False
        self.Units_Factor = 0.0

        # GLOBAL ARRAYS

        # Float
        self.Helium_Time_Constant = [0.0 for i in range(16)]
        self.Nitrogen_Time_Constant = [0.0 for i in range(16)]

        self.Helium_Pressure = [0.0 for i in range(16)]
        self.Nitrogen_Pressure = [0.0 for i in range(16)]

        self.Initial_Helium_Pressure = [0.0 for i in range(16)]
        self.Initial_Nitrogen_Pressure = [0.0 for i in range(16)]

        # self.Fraction_Oxygen = []
        self.Fraction_Helium = []
        self.Fraction_Nitrogen = []

        self.Initial_Critical_Radius_He = [0.0 for i in range(16)]
        self.Initial_Critical_Radius_N2 = [0.0 for i in range(16)]
        self.Adjusted_Critical_Radius_He = [0.0 for i in range(16)]
        self.Adjusted_Critical_Radius_N2 = [0.0 for i in range(16)]

        self.Max_Crushing_Pressure_He = [0.0 for i in range(16)]
        self.Max_Crushing_Pressure_N2 = [0.0 for i in range(16)]

        self.Surface_Phase_Volume_Time = [0.0 for i in range(16)]

        self.Max_Actual_Gradient = [0.0 for i in range(16)]

        self.Amb_Pressure_Onset_of_Imperm = [0.0 for i in range(16)]
        self.Gas_Tension_Onset_of_Imperm = [0.0 for i in range(16)]

        self.Diver_Acclimatized = None

        # ASSIGN HALF-TIME VALUES TO BUHLMANN COMPARTMENT ARRAYS
        self.Helium_Half_Time = [1.88, 3.02, 4.72, 6.99, 10.21, 14.48, 20.53, 29.11, 41.20,
                                 55.19, 70.69, 90.34, 115.29, 147.42, 188.24, 240.03]

        self.Nitrogen_Half_Time = [5.0, 8.0, 12.5, 18.5, 27.0, 38.3, 54.3, 77.0, 109.0,
                                   146.0, 187.0, 239.0, 305.0, 390.0, 498.0, 635.0]

    def get_json(self):
        return self.output_object.get_json()

    def GAS_LOADINGS_SURFACE_INTERVAL(self, Surface_Interval_Time):
        """
        Purpose: This subprogram calculates the gas loading (off-gassing) during
        a surface interval.

        Side Effects: Sets
        `self.Helium_Pressure`,
        `self.Nitrogen_Pressure`

        Returns: None
        """

        Inspired_Helium_Pressure = 0.0
        Inspired_Nitrogen_Pressure = (self.Barometric_Pressure - self.Water_Vapor_Pressure) * self.fraction_inert_gas

        for i in range(16):
            Temp_Helium_Pressure = self.Helium_Pressure[i]
            Temp_Nitrogen_Pressure = self.Nitrogen_Pressure[i]

            self.Helium_Pressure[i] = HALDANE_EQUATION(Temp_Helium_Pressure, Inspired_Helium_Pressure, self.Helium_Time_Constant[i], Surface_Interval_Time)

            self.Nitrogen_Pressure[i] = HALDANE_EQUATION(Temp_Nitrogen_Pressure, Inspired_Nitrogen_Pressure, self.Nitrogen_Time_Constant[i], Surface_Interval_Time)

    def _new_critical_radius(self, Max_Actual_Gradient_Pascals, Adj_Crush_Pressure_Pascals):
        """Calculates the new radius for the `VPM_REPETITIVE_ALGORITHM`

        Side Effects: None
        Returns: A floating point value
        """
        return ((2.0 * self.Surface_Tension_Gamma * (self.Skin_Compression_GammaC - self.Surface_Tension_Gamma))) / (Max_Actual_Gradient_Pascals * self.Skin_Compression_GammaC - self.Surface_Tension_Gamma * Adj_Crush_Pressure_Pascals)

    def VPM_REPETITIVE_ALGORITHM(self, Surface_Interval_Time):
        """
        Purpose: This subprogram implements the VPM Repetitive Algorithm that was
        envisioned by Professor David E. Yount only months before his passing.

        Side Effects: Sets
        `self.Adjusted_Critical_Radius_He`,
        `self.Adjusted_Critical_Radius_N2`

        Returns: None
        """

        for i in range(16):
            Max_Actual_Gradient_Pascals = (self.Max_Actual_Gradient[i] / self.Units_Factor) * self.ATM

            Adj_Crush_Pressure_He_Pascals = (self.Adjusted_Crushing_Pressure_He[i] / self.Units_Factor) * self.ATM
            Adj_Crush_Pressure_N2_Pascals = (self.Adjusted_Crushing_Pressure_N2[i] / self.Units_Factor) * self.ATM

            if self.Max_Actual_Gradient[i] > self.Initial_Allowable_Gradient_N2[i]:
                New_Critical_Radius_N2 = self._new_critical_radius(Max_Actual_Gradient_Pascals, Adj_Crush_Pressure_N2_Pascals)

                self.Adjusted_Critical_Radius_N2[i] = self.Initial_Critical_Radius_N2[i] + (self.Initial_Critical_Radius_N2[i] - New_Critical_Radius_N2) * exp(-Surface_Interval_Time / self.Regeneration_Time_Constant)

            else:
                self.Adjusted_Critical_Radius_N2[i] = self.Initial_Critical_Radius_N2[i]

            if self.Max_Actual_Gradient[i] > self.Initial_Allowable_Gradient_He[i]:
                New_Critical_Radius_He = self._new_critical_radius(Max_Actual_Gradient_Pascals, Adj_Crush_Pressure_He_Pascals)

                self.Adjusted_Critical_Radius_He[i] = self.Initial_Critical_Radius_He[i] + (self.Initial_Critical_Radius_He[i] - New_Critical_Radius_He) * exp(-Surface_Interval_Time / self.Regeneration_Time_Constant)
            else:
                self.Adjusted_Critical_Radius_He[i] = self.Initial_Critical_Radius_He[i]

    def CALC_MAX_ACTUAL_GRADIENT(self, Deco_Stop_Depth):
        """
        Purpose: This subprogram calculates the actual supersaturation gradient
        obtained in each compartment as a result of the ascent profile during
        decompression.  Similar to the concept with crushing pressure, the
        supersaturation gradients are not cumulative over a multi-level, staged
        ascent.  Rather, it will be the maximum value obtained in any one discrete
        step of the overall ascent.  Thus, the program must compute and store the
        maximum actual gradient for each compartment that was obtained across all
        steps of the ascent profile.  This subroutine is invoked on the last pass
        through the deco stop loop block when the final deco schedule is being
        generated.

        The max actual gradients are later used by the VPM Repetitive Algorithm to
        determine if adjustments to the critical radii are required.  If the max
        actual gradient did not exceed the initial alllowable gradient, then no
        adjustment will be made.  However, if the max actual gradient did exceed
        the intitial allowable gradient, such as permitted by the Critical Volume
        Algorithm, then the critical radius will be adjusted (made larger) on the
        repetitive dive to compensate for the bubbling that was allowed on the
        previous dive.  The use of the max actual gradients is intended to prevent
        the repetitive algorithm from being overly conservative.

        Side Effects: Sets
        `self.Max_Actual_Gradient`

        Returns: None
        """

        # Note: negative supersaturation gradients are meaningless for this
        # application, so the values must be equal to or greater than zero.

        for i in range(16):
            Compartment_Gradient = (self.Helium_Pressure[i] + self.Nitrogen_Pressure[i] + self.Constant_Pressure_Other_Gases) - (Deco_Stop_Depth + self.Barometric_Pressure)
            if Compartment_Gradient <= 0.0:
                Compartment_Gradient = 0.0

            self.Max_Actual_Gradient[i] = max(self.Max_Actual_Gradient[i], Compartment_Gradient)

    def VPM_ALTITUDE_DIVE_ALGORITHM(self, altitude_settings):
        """
        Purpose:  This subprogram updates gas loadings and adjusts critical radii
        (as required) based on whether or not diver is acclimatized at altitude or
        makes an ascent to altitude before the dive.

        Side Effects: Sets
        `self.Adjusted_Critical_Radius_He`,
        `self.Adjusted_Critical_Radius_N2`,
        `self.Barometric_Pressure`,
        `self.Helium_Pressure`,
        `self.Initial_Critical_Radius_He`,
        `self.Initial_Critical_Radius_N2`
        `self.Nitrogen_Pressure`,

        or

        Raises an AltitudeException

        Returns: None
        """

        Ascent_to_Altitude_Time = altitude_settings['Ascent_to_Altitude_Hours'] * 60.0
        Time_at_Altitude_Before_Dive = altitude_settings['Hours_at_Altitude_Before_Dive'] * 60.0

        if self.Diver_Acclimatized:
            self.Barometric_Pressure = CALC_BAROMETRIC_PRESSURE(altitude_settings['Altitude_of_Dive'], self.units_fsw)

            for i in range(16):
                self.Adjusted_Critical_Radius_N2[i] = self.Initial_Critical_Radius_N2[i]
                self.Adjusted_Critical_Radius_He[i] = self.Initial_Critical_Radius_He[i]
                self.Helium_Pressure[i] = 0.0
                self.Nitrogen_Pressure[i] = (self.Barometric_Pressure - self.Water_Vapor_Pressure) * self.fraction_inert_gas
        else:
            if (altitude_settings['Starting_Acclimatized_Altitude'] >= altitude_settings['Altitude_of_Dive']) or (altitude_settings['Starting_Acclimatized_Altitude'] < 0.0):
                raise AltitudeException("ERROR! STARTING ACCLIMATIZED ALTITUDE MUST BE LESS THAN ALTITUDE OF DIVE AND GREATER THAN OR EQUAL TO ZERO")

            self.Barometric_Pressure = CALC_BAROMETRIC_PRESSURE(altitude_settings['Starting_Acclimatized_Altitude'], self.units_fsw)

            Starting_Ambient_Pressure = self.Barometric_Pressure

            for i in range(16):
                self.Helium_Pressure[i] = 0.0
                self.Nitrogen_Pressure[i] = (self.Barometric_Pressure - self.Water_Vapor_Pressure) * self.fraction_inert_gas

            self.Barometric_Pressure = CALC_BAROMETRIC_PRESSURE(altitude_settings['Altitude_of_Dive'], self.units_fsw)
            Ending_Ambient_Pressure = self.Barometric_Pressure
            Initial_Inspired_N2_Pressure = (Starting_Ambient_Pressure - self.Water_Vapor_Pressure) * self.fraction_inert_gas
            Rate = (Ending_Ambient_Pressure - Starting_Ambient_Pressure) / Ascent_to_Altitude_Time
            Nitrogen_Rate = Rate * self.fraction_inert_gas

            for i in range(16):
                Initial_Nitrogen_Pressure = self.Nitrogen_Pressure[i]

                self.Nitrogen_Pressure[i] = SCHREINER_EQUATION(Initial_Inspired_N2_Pressure, Nitrogen_Rate, Ascent_to_Altitude_Time, self.Nitrogen_Time_Constant[i], Initial_Nitrogen_Pressure)

                Compartment_Gradient = (self.Nitrogen_Pressure[i] + self.Constant_Pressure_Other_Gases) - Ending_Ambient_Pressure

                Compartment_Gradient_Pascals = (Compartment_Gradient / self.Units_Factor) * self.ATM

                Gradient_He_Bubble_Formation = ((2.0 * self.Surface_Tension_Gamma * (self.Skin_Compression_GammaC - self.Surface_Tension_Gamma)) / (self.Initial_Critical_Radius_He[i] * self.Skin_Compression_GammaC))

                if Compartment_Gradient_Pascals > Gradient_He_Bubble_Formation:

                    New_Critical_Radius_He = ((2.0 * self.Surface_Tension_Gamma * (self.Skin_Compression_GammaC - self.Surface_Tension_Gamma))) / (Compartment_Gradient_Pascals * self.Skin_Compression_GammaC)

                    self.Adjusted_Critical_Radius_He[i] = self.Initial_Critical_Radius_He[i] + (self.Initial_Critical_Radius_He[i] - New_Critical_Radius_He) * exp(-Time_at_Altitude_Before_Dive / self.Regeneration_Time_Constant)

                    self.Initial_Critical_Radius_He[i] = self.Adjusted_Critical_Radius_He[i]

                else:
                    Ending_Radius_He = 1.0 / (Compartment_Gradient_Pascals / (2.0 * (self.Surface_Tension_Gamma - self.Skin_Compression_GammaC)) + 1.0 / self.Initial_Critical_Radius_He[i])

                    Regenerated_Radius_He = self.Initial_Critical_Radius_He[i] + (Ending_Radius_He - self.Initial_Critical_Radius_He[i]) * exp(-Time_at_Altitude_Before_Dive / self.Regeneration_Time_Constant)

                    self.Initial_Critical_Radius_He[i] = Regenerated_Radius_He

                    self.Adjusted_Critical_Radius_He[i] = self.Initial_Critical_Radius_He[i]

                Gradient_N2_Bubble_Formation = ((2.0 * self.Surface_Tension_Gamma * (self.Skin_Compression_GammaC - self.Surface_Tension_Gamma)) / (self.Initial_Critical_Radius_N2[i] * self.Skin_Compression_GammaC))

                if Compartment_Gradient_Pascals > Gradient_N2_Bubble_Formation:

                    New_Critical_Radius_N2 = ((2.0 * self.Surface_Tension_Gamma * (self.Skin_Compression_GammaC - self.Surface_Tension_Gamma))) / (Compartment_Gradient_Pascals * self.Skin_Compression_GammaC)

                    self.Adjusted_Critical_Radius_N2[i] = self.Initial_Critical_Radius_N2[i] + (self.Initial_Critical_Radius_N2[i] - New_Critical_Radius_N2) * exp(-Time_at_Altitude_Before_Dive / self.Regeneration_Time_Constant)

                    self.Initial_Critical_Radius_N2[i] = self.Adjusted_Critical_Radius_N2[i]

                else:
                    Ending_Radius_N2 = 1.0 / (Compartment_Gradient_Pascals / (2.0 * (self.Surface_Tension_Gamma - self.Skin_Compression_GammaC)) + 1.0 / self.Initial_Critical_Radius_N2[i])

                    Regenerated_Radius_N2 = self.Initial_Critical_Radius_N2[i] + (Ending_Radius_N2 - self.Initial_Critical_Radius_N2[i]) * exp(-Time_at_Altitude_Before_Dive / self.Regeneration_Time_Constant)

                    self.Initial_Critical_Radius_N2[i] = Regenerated_Radius_N2

                    self.Adjusted_Critical_Radius_N2[i] = self.Initial_Critical_Radius_N2[i]

            Inspired_Nitrogen_Pressure = (self.Barometric_Pressure - self.Water_Vapor_Pressure) * self.fraction_inert_gas

            for i in range(16):
                Initial_Nitrogen_Pressure = self.Nitrogen_Pressure[i]

                self.Nitrogen_Pressure[i] = HALDANE_EQUATION(Initial_Nitrogen_Pressure, Inspired_Nitrogen_Pressure, self.Nitrogen_Time_Constant[i], Time_at_Altitude_Before_Dive)

    def GAS_LOADINGS_ASCENT_DESCENT(self, Starting_Depth, Ending_Depth, Rate):
        """
         Purpose: This subprogram applies the Schreiner equation to update the
         gas loadings (partial pressures of helium and nitrogen) in the half-time
         compartments due to a linear ascent or descent segment at a constant rate.

         Side Effects: Sets `self.Segment_Time`,
         `self.Ending_Ambient_Pressure`,
         `self.Helium_Pressure`,
         `self.Initial_Helium_Pressure`,
         `self.Initial_Nitrogen_Pressure`,
         `self.Nitrogen_Pressure`
         `self.Run_Time`,
         `self.Segment_Number`,

         Returns: None
         """
        self.Segment_Time = float(Ending_Depth - Starting_Depth) / Rate

        self.Run_Time = self.Run_Time + self.Segment_Time
        self.Segment_Number = self.Segment_Number + 1
        self.Ending_Ambient_Pressure = Ending_Depth + self.Barometric_Pressure

        Starting_Ambient_Pressure = Starting_Depth + self.Barometric_Pressure
        Initial_Inspired_He_Pressure = (Starting_Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Helium[self.Mix_Number - 1]
        Initial_Inspired_N2_Pressure = (Starting_Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Nitrogen[self.Mix_Number - 1]
        Helium_Rate = Rate * self.Fraction_Helium[self.Mix_Number - 1]
        Nitrogen_Rate = Rate * self.Fraction_Nitrogen[self.Mix_Number - 1]

        for i in range(16):
            self.Initial_Helium_Pressure[i] = self.Helium_Pressure[i]
            self.Initial_Nitrogen_Pressure[i] = self.Nitrogen_Pressure[i]

            self.Helium_Pressure[i] = SCHREINER_EQUATION(Initial_Inspired_He_Pressure, Helium_Rate, self.Segment_Time, self.Helium_Time_Constant[i], self.Initial_Helium_Pressure[i])

            self.Nitrogen_Pressure[i] = SCHREINER_EQUATION(Initial_Inspired_N2_Pressure, Nitrogen_Rate, self.Segment_Time, self.Nitrogen_Time_Constant[i], self.Initial_Nitrogen_Pressure[i])

    def _crushing_pressure_helper(self, Radius_Onset_of_Imperm_Molecule, Ending_Ambient_Pressure_Pa, Amb_Press_Onset_of_Imperm_Pa, Gas_Tension_Onset_of_Imperm_Pa, Gradient_Onset_of_Imperm_Pa):
        """Calculate the crushing pressure for a molecule(He or N2) (a helper for CALC_CRUSHING_PRESSURE)

        Side Effects: None

        Returns: A floating point value
        """

        A = Ending_Ambient_Pressure_Pa - Amb_Press_Onset_of_Imperm_Pa + Gas_Tension_Onset_of_Imperm_Pa + (2.0 * (self.Skin_Compression_GammaC - self.Surface_Tension_Gamma)) / Radius_Onset_of_Imperm_Molecule
        B = 2.0 * (self.Skin_Compression_GammaC - self.Surface_Tension_Gamma)
        C = Gas_Tension_Onset_of_Imperm_Pa * Radius_Onset_of_Imperm_Molecule ** 3

        High_Bound = Radius_Onset_of_Imperm_Molecule
        Low_Bound = B / A

        Ending_Radius = RADIUS_ROOT_FINDER(A, B, C, Low_Bound, High_Bound)
        Crushing_Pressure_Pascals = Gradient_Onset_of_Imperm_Pa + Ending_Ambient_Pressure_Pa - Amb_Press_Onset_of_Imperm_Pa + Gas_Tension_Onset_of_Imperm_Pa * (1.0 - Radius_Onset_of_Imperm_Molecule ** 3 / Ending_Radius ** 3)

        return (Crushing_Pressure_Pascals / self.ATM) * self.Units_Factor

    def CALC_CRUSHING_PRESSURE(self, Starting_Depth, Ending_Depth, Rate):
        """
         Purpose: Compute the effective "crushing pressure" in each compartment as
         a result of descent segment(s).  The crushing pressure is the gradient
         (difference in pressure) between the outside ambient pressure and the
         gas tension inside a VPM nucleus (bubble seed).  This gradient acts to
         reduce (shrink) the radius smaller than its initial value at the surface.
         This phenomenon has important ramifications because the smaller the radius
         of a VPM nucleus, the greater the allowable supersaturation gradient upon
         ascent.  Gas loading (uptake) during descent, especially in the fast
         compartments, will reduce the magnitude of the crushing pressure.  The
         crushing pressure is not cumulative over a multi-level descent.  It will
         be the maximum value obtained in any one discrete segment of the overall
         descent.  Thus, the program must compute and store the maximum crushing
         pressure for each compartment that was obtained across all segments of
         the descent profile.

         The calculation of crushing pressure will be different depending on
         whether or not the gradient is in the VPM permeable range (gas can diffuse
         across skin of VPM nucleus) or the VPM impermeable range (molecules in
         skin of nucleus are squeezed together so tight that gas can no longer
         diffuse in or out of nucleus; the gas becomes trapped and further resists
         the crushing pressure).  The solution for crushing pressure in the VPM
         permeable range is a simple linear equation.  In the VPM impermeable
         range, a cubic equation must be solved using a numerical method.

         Separate crushing pressures are tracked for helium and nitrogen because
         they can have different critical radii.  The crushing pressures will be
         the same for helium and nitrogen in the permeable range of the model, but
         they will start to diverge in the impermeable range.  This is due to
         the differences between starting radius, radius at the onset of
         impermeability, and radial compression in the impermeable range.

         Side Effects: Sets
         `self.Max_Crushing_Pressure_He`,
         `self.Max_Crushing_Pressure_N2`

         Returns: None
         """

        # First, convert the Gradient for Onset of Impermeability from units of
        # atmospheres to diving pressure units (either fsw or msw) and to Pascals
        # (SI units).  The reason that the Gradient for Onset of Impermeability is
        # given in the program settings in units of atmospheres is because that is
        # how it was reported in the original research papers by Yount and
        # colleauges.

        Gradient_Onset_of_Imperm = self.settings_values['Gradient_Onset_of_Imperm_Atm'] * self.Units_Factor  # convert to diving units
        Gradient_Onset_of_Imperm_Pa = self.settings_values['Gradient_Onset_of_Imperm_Atm'] * self.ATM     # convert to Pascals

        # Assign values of starting and ending ambient pressures for descent segment

        Starting_Ambient_Pressure = Starting_Depth + self.Barometric_Pressure
        Ending_Ambient_Pressure = Ending_Depth + self.Barometric_Pressure

        # MAIN LOOP WITH NESTED DECISION TREE
        # For each compartment, the program computes the starting and ending
        # gas tensions and gradients.  The VPM is different than some dissolved gas
        # algorithms, Buhlmann for example, in that it considers the pressure due to
        # oxygen, carbon dioxide, and water vapor in each compartment in addition to
        # the inert gases helium and nitrogen.  These "other gases" are included in
        # the calculation of gas tensions and gradients.

        for i in range(16):
            Starting_Gas_Tension = self.Initial_Helium_Pressure[i] + self.Initial_Nitrogen_Pressure[i] + self.Constant_Pressure_Other_Gases

            Starting_Gradient = Starting_Ambient_Pressure - Starting_Gas_Tension

            Ending_Gas_Tension = self.Helium_Pressure[i] + self.Nitrogen_Pressure[i] + self.Constant_Pressure_Other_Gases

            Ending_Gradient = Ending_Ambient_Pressure - Ending_Gas_Tension

            # Compute radius at onset of impermeability for helium and nitrogen
            # critical radii

            Radius_Onset_of_Imperm_He = 1.0 / (Gradient_Onset_of_Imperm_Pa / (2.0 * (self.settings_values['Skin_Compression_GammaC'] - self.settings_values['Surface_Tension_Gamma'])) + 1.0 / self.Adjusted_Critical_Radius_He[i])

            Radius_Onset_of_Imperm_N2 = 1.0 / (Gradient_Onset_of_Imperm_Pa / (2.0 * (self.settings_values['Skin_Compression_GammaC'] - self.settings_values['Surface_Tension_Gamma'])) + 1.0 / self.Adjusted_Critical_Radius_N2[i])

            # FIRST BRANCH OF DECISION TREE - PERMEABLE RANGE
            # Crushing pressures will be the same for helium and nitrogen
            if Ending_Gradient <= Gradient_Onset_of_Imperm:
                Crushing_Pressure_He = Ending_Ambient_Pressure - Ending_Gas_Tension
                Crushing_Pressure_N2 = Ending_Ambient_Pressure - Ending_Gas_Tension

            # SECOND BRANCH OF DECISION TREE - IMPERMEABLE RANGE
            # Both the ambient pressure and the gas tension at the onset of
            # impermeability must be computed in order to properly solve for the ending
            # radius and resultant crushing pressure.  The first decision block
            # addresses the special case when the starting gradient just happens to be
            # equal to the gradient for onset of impermeability (not very likely!).

            # if Ending_Gradient > Gradient_Onset_of_Imperm:
            else:

                if Starting_Gradient == Gradient_Onset_of_Imperm:
                    self.Amb_Pressure_Onset_of_Imperm[i] = Starting_Ambient_Pressure
                    self.Gas_Tension_Onset_of_Imperm[i] = Starting_Gas_Tension
                # In most cases, a subroutine will be called to find these values using a
                # numerical method.
                if Starting_Gradient < Gradient_Onset_of_Imperm:
                    self.ONSET_OF_IMPERMEABILITY(Starting_Ambient_Pressure, Ending_Ambient_Pressure, Rate, i)

                # Next, using the values for ambient pressure and gas tension at the onset
                # of impermeability, the equations are set up to process the calculations
                # through the radius root finder subroutine.  This subprogram will find the
                # root (solution) to the cubic equation using a numerical method.  In order
                # to do this efficiently, the equations are placed in the form
                # Ar^3 - Br^2 - C = 0, where r is the ending radius after impermeable
                # compression.  The coefficients A, B, and C for helium and nitrogen are
                # computed and passed to the subroutine as arguments.  The high and low
                # bounds to be used by the numerical method of the subroutine are also
                # computed (see separate page posted on Deco List ftp site entitled
                # "VPM: Solving for radius in the impermeable regime").  The subprogram
                # will return the value of the ending radius and then the crushing
                # pressures for helium and nitrogen can be calculated.
                Ending_Ambient_Pressure_Pa = (Ending_Ambient_Pressure / self.Units_Factor) * self.ATM

                Amb_Press_Onset_of_Imperm_Pa = (self.Amb_Pressure_Onset_of_Imperm[i] / self.Units_Factor) * self.ATM

                Gas_Tension_Onset_of_Imperm_Pa = (self.Gas_Tension_Onset_of_Imperm[i] / self.Units_Factor) * self.ATM

                Crushing_Pressure_He = self._crushing_pressure_helper(Radius_Onset_of_Imperm_He, Ending_Ambient_Pressure_Pa, Amb_Press_Onset_of_Imperm_Pa, Gas_Tension_Onset_of_Imperm_Pa, Gradient_Onset_of_Imperm_Pa)

                Crushing_Pressure_N2 = self._crushing_pressure_helper(Radius_Onset_of_Imperm_N2, Ending_Ambient_Pressure_Pa, Amb_Press_Onset_of_Imperm_Pa, Gas_Tension_Onset_of_Imperm_Pa, Gradient_Onset_of_Imperm_Pa)

            # UPDATE VALUES OF MAX CRUSHING PRESSURE IN Object ARRAYS
            self.Max_Crushing_Pressure_He[i] = max(self.Max_Crushing_Pressure_He[i], Crushing_Pressure_He)
            self.Max_Crushing_Pressure_N2[i] = max(self.Max_Crushing_Pressure_N2[i], Crushing_Pressure_N2)

    def ONSET_OF_IMPERMEABILITY(self, Starting_Ambient_Pressure, Ending_Ambient_Pressure, Rate, i):
        """
        Purpose:  This subroutine uses the Bisection Method to find the ambient
        pressure and gas tension at the onset of impermeability for a given
        compartment.  Source:  "Numerical Recipes in Fortran 77",
        Cambridge University Press, 1992.

        Side Effects: Sets
        `self.Amb_Pressure_Onset_of_Imperm`,
        `self.Gas_Tension_Onset_of_Imperm`

        or

        Raises a RootException

        Returns: None
        """
        # First convert the Gradient for Onset of Impermeability to the diving
        # pressure units that are being used

        Gradient_Onset_of_Imperm = self.Gradient_Onset_of_Imperm_Atm * self.Units_Factor

        # ESTABLISH THE BOUNDS FOR THE ROOT SEARCH USING THE BISECTION METHOD
        # In this case, we are solving for time - the time when the ambient pressure
        # minus the gas tension will be equal to the Gradient for Onset of
        # Impermeabliity.  The low bound for time is set at zero and the high
        # bound is set at the elapsed time (segment time) it took to go from the
        # starting ambient pressure to the ending ambient pressure.  The desired
        # ambient pressure and gas tension at the onset of impermeability will
        # be found somewhere between these endpoints.  The algorithm checks to
        # make sure that the solution lies in between these bounds by first
        # computing the low bound and high bound function values.

        Initial_Inspired_He_Pressure = (Starting_Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Helium[self.Mix_Number - 1]

        Initial_Inspired_N2_Pressure = (Starting_Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Nitrogen[self.Mix_Number - 1]

        Helium_Rate = Rate * self.Fraction_Helium[self.Mix_Number - 1]
        Nitrogen_Rate = Rate * self.Fraction_Nitrogen[self.Mix_Number - 1]
        Low_Bound = 0.0

        High_Bound = (Ending_Ambient_Pressure - Starting_Ambient_Pressure) / Rate

        Starting_Gas_Tension = self.Initial_Helium_Pressure[i] + self.Initial_Nitrogen_Pressure[i] + self.Constant_Pressure_Other_Gases

        Function_at_Low_Bound = Starting_Ambient_Pressure - Starting_Gas_Tension - Gradient_Onset_of_Imperm

        High_Bound_Helium_Pressure = SCHREINER_EQUATION(Initial_Inspired_He_Pressure, Helium_Rate, High_Bound, self.Helium_Time_Constant[i], self.Initial_Helium_Pressure[i])

        High_Bound_Nitrogen_Pressure = SCHREINER_EQUATION(Initial_Inspired_N2_Pressure, Nitrogen_Rate, High_Bound, self.Nitrogen_Time_Constant[i], self.Initial_Nitrogen_Pressure[i])

        Ending_Gas_Tension = High_Bound_Helium_Pressure + High_Bound_Nitrogen_Pressure + self.Constant_Pressure_Other_Gases

        Function_at_High_Bound = Ending_Ambient_Pressure - Ending_Gas_Tension - Gradient_Onset_of_Imperm

        if(Function_at_High_Bound * Function_at_Low_Bound) >= 0.0:
            raise RootException("ERROR! ROOT IS NOT WITHIN BRACKETS")

        # APPLY THE BISECTION METHOD IN SEVERAL ITERATIONS UNTIL A SOLUTION WITH
        # THE DESIRED ACCURACY IS FOUND
        # Note: the program allows for up to 100 iterations.  Normally an exit will
        # be made from the loop well before that number.  If, for some reason, the
        # program exceeds 100 iterations, there will be a pause to alert the user.

        if Function_at_Low_Bound < 0.0:
            Time = Low_Bound
            Differential_Change = High_Bound - Low_Bound
        else:
            Time = High_Bound
            Differential_Change = Low_Bound - High_Bound

        for j in range(100):
            Last_Diff_Change = Differential_Change
            Differential_Change = Last_Diff_Change * 0.5
            Mid_Range_Time = Time + Differential_Change

            Mid_Range_Ambient_Pressure = (Starting_Ambient_Pressure + Rate * Mid_Range_Time)

            Mid_Range_Helium_Pressure = SCHREINER_EQUATION(Initial_Inspired_He_Pressure, Helium_Rate, Mid_Range_Time, self.Helium_Time_Constant[i], self.Initial_Helium_Pressure[i])

            Mid_Range_Nitrogen_Pressure = SCHREINER_EQUATION(Initial_Inspired_N2_Pressure, Nitrogen_Rate, Mid_Range_Time, self.Nitrogen_Time_Constant[i], self.Initial_Nitrogen_Pressure[i])

            Gas_Tension_at_Mid_Range = Mid_Range_Helium_Pressure + Mid_Range_Nitrogen_Pressure + self.Constant_Pressure_Other_Gases

            Function_at_Mid_Range = Mid_Range_Ambient_Pressure - Gas_Tension_at_Mid_Range - Gradient_Onset_of_Imperm

            if Function_at_Mid_Range <= 0.0:
                Time = Mid_Range_Time

            # When a solution with the desired accuracy is found, the program breaks
            if (abs(Differential_Change) < 1.0E-3) or (Function_at_Mid_Range == 0.0):
                break

        self.Amb_Pressure_Onset_of_Imperm[i] = Mid_Range_Ambient_Pressure
        self.Gas_Tension_Onset_of_Imperm[i] = Gas_Tension_at_Mid_Range

    def GAS_LOADINGS_CONSTANT_DEPTH(self, Depth, Run_Time_End_of_Segment):
        """
        Purpose: This subprogram applies the Haldane equation to update the
        gas loadings (partial pressures of helium and nitrogen) in the half-time
        compartments for a segment at constant depth.

        Side Effects: Sets
        `self.Ending_Ambient_Pressure`,
        `self.Helium_Pressure`,
        `self.Nitrogen_Pressure`
        `self.Run_Time`,
        `self.Segment_Number`,
        `self.Segment_Time`,

        Returns: None
        """

        self.Segment_Time = Run_Time_End_of_Segment - self.Run_Time
        self.Run_Time = Run_Time_End_of_Segment
        self.Segment_Number = self.Segment_Number + 1
        Ambient_Pressure = Depth + self.Barometric_Pressure

        Inspired_Helium_Pressure = (Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Helium[self.Mix_Number - 1]

        Inspired_Nitrogen_Pressure = (Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Nitrogen[self.Mix_Number - 1]

        self.Ending_Ambient_Pressure = Ambient_Pressure

        Temp_Helium_Pressure = 0.0
        Temp_Nitrogen_Pressure = 0.0

        for i in range(16):
            Temp_Helium_Pressure = self.Helium_Pressure[i]
            Temp_Nitrogen_Pressure = self.Nitrogen_Pressure[i]

            self.Helium_Pressure[i] = HALDANE_EQUATION(Temp_Helium_Pressure, Inspired_Helium_Pressure, self.Helium_Time_Constant[i], self.Segment_Time)

            self.Nitrogen_Pressure[i] = HALDANE_EQUATION(Temp_Nitrogen_Pressure, Inspired_Nitrogen_Pressure, self.Nitrogen_Time_Constant[i], self.Segment_Time)

    def NUCLEAR_REGENERATION(self, Dive_Time):
        """
        Purpose: This subprogram calculates the regeneration of VPM critical
        radii that takes place over the dive time.  The regeneration time constant
        has a time scale of weeks so this will have very little impact on dives of
        normal length, but will have a major impact for saturation dives.

        Side Effects: Sets
        `self.Adjusted_Crushing_Pressure_He`,
        `self.Adjusted_Crushing_Pressure_N2`
        `self.Regenerated_Radius_He`,
        `self.Regenerated_Radius_N2`,

        Returns: None
        """

        # First convert the maximum crushing pressure obtained for each compartment
        # to Pascals.  Next, compute the ending radius for helium and nitrogen
        # critical nuclei in each compartment.
        for i in range(16):
            Crushing_Pressure_Pascals_He = (self.Max_Crushing_Pressure_He[i] / self.Units_Factor) * self.ATM

            Crushing_Pressure_Pascals_N2 = (self.Max_Crushing_Pressure_N2[i] / self.Units_Factor) * self.ATM

            Ending_Radius_He = 1.0 / (Crushing_Pressure_Pascals_He / (2.0 * (self.settings_values['Skin_Compression_GammaC'] - self.settings_values['Surface_Tension_Gamma'])) + 1.0 / self.Adjusted_Critical_Radius_He[i])

            Ending_Radius_N2 = 1.0 / (Crushing_Pressure_Pascals_N2 / (2.0 * (self.settings_values['Skin_Compression_GammaC'] - self.settings_values['Surface_Tension_Gamma'])) + 1.0 / self.Adjusted_Critical_Radius_N2[i])
            # A "regenerated" radius for each nucleus is now calculated based on the
            # regeneration time constant.  This means that after application of
            # crushing pressure and reduction in radius, a nucleus will slowly grow
            # back to its original initial radius over a period of time.  This
            # phenomenon is probabilistic in nature and depends on absolute temperature.
            # It is independent of crushing pressure.
            self.Regenerated_Radius_He[i] = self.Adjusted_Critical_Radius_He[i] + (Ending_Radius_He - self.Adjusted_Critical_Radius_He[i]) * exp(-Dive_Time / self.settings_values['Regeneration_Time_Constant'])

            self.Regenerated_Radius_N2[i] = self.Adjusted_Critical_Radius_N2[i] + (Ending_Radius_N2 - self.Adjusted_Critical_Radius_N2[i]) * exp(-Dive_Time / self.settings_values['Regeneration_Time_Constant'])

            # In order to preserve reference back to the initial critical radii after
            # regeneration, an "adjusted crushing pressure" for the nuclei in each
            # compartment must be computed.  In other words, this is the value of
            # crushing pressure that would have reduced the original nucleus to the
            # to the present radius had regeneration not taken place.  The ratio
            # for adjusting crushing pressure is obtained from algebraic manipulation
            # of the standard VPM equations.  The adjusted crushing pressure, in lieu
            # of the original crushing pressure, is then applied in the VPM Critical
            # Volume Algorithm and the VPM Repetitive Algorithm.

            Crush_Pressure_Adjust_Ratio_He = (Ending_Radius_He * (self.Adjusted_Critical_Radius_He[i] - self.Regenerated_Radius_He[i])) / (self.Regenerated_Radius_He[i] * (self.Adjusted_Critical_Radius_He[i] - Ending_Radius_He))

            Crush_Pressure_Adjust_Ratio_N2 = (Ending_Radius_N2 * (self.Adjusted_Critical_Radius_N2[i] - self.Regenerated_Radius_N2[i])) / (self.Regenerated_Radius_N2[i] * (self.Adjusted_Critical_Radius_N2[i] - Ending_Radius_N2))

            Adj_Crush_Pressure_He_Pascals = Crushing_Pressure_Pascals_He * Crush_Pressure_Adjust_Ratio_He
            Adj_Crush_Pressure_N2_Pascals = Crushing_Pressure_Pascals_N2 * Crush_Pressure_Adjust_Ratio_N2

            self.Adjusted_Crushing_Pressure_He[i] = (Adj_Crush_Pressure_He_Pascals / self.ATM) * self.Units_Factor
            self.Adjusted_Crushing_Pressure_N2[i] = (Adj_Crush_Pressure_N2_Pascals / self.ATM) * self.Units_Factor

    def CALC_INITIAL_ALLOWABLE_GRADIENT(self):
        """
        Purpose: This subprogram calculates the initial allowable gradients for
        helium and nitrogren in each compartment.  These are the gradients that
        will be used to set the deco ceiling on the first pass through the deco
        loop.  If the Critical Volume Algorithm is set to "off", then these
        gradients will determine the final deco schedule.  Otherwise, if the
        Critical Volume Algorithm is set to "on", these gradients will be further
        "relaxed" by the Critical Volume Algorithm subroutine.  The initial
        allowable gradients are referred to as "PssMin" in the papers by Yount
        and colleauges, i.e., the minimum supersaturation pressure gradients
        that will probe bubble formation in the VPM nuclei that started with the
        designated minimum initial radius (critical radius).

        The initial allowable gradients are computed directly from the
        "regenerated" radii after the Nuclear Regeneration subroutine.  These
        gradients are tracked separately for helium and nitrogen.

        Side Effects: Sets
        `self.Allowable_Gradient_He`,
        `self.Allowable_Gradient_N2`
        `self.Initial_Allowable_Gradient_He`,
        `self.Initial_Allowable_Gradient_N2`,

        Returns: None
        """

        # The initial allowable gradients are computed in Pascals and then converted
        # to the diving pressure units.  Two different sets of arrays are used to
        # save the calculations - Initial Allowable Gradients and Allowable
        # Gradients.  The Allowable Gradients are assigned the values from Initial
        # Allowable Gradients however the Allowable Gradients can be changed later
        # by the Critical Volume subroutine.  The values for the Initial Allowable
        # Gradients are saved in a global array for later use by both the Critical
        # Volume subroutine and the VPM Repetitive Algorithm subroutine.

        for i in range(16):
            Initial_Allowable_Grad_N2_Pa = ((2.0 * self.settings_values['Surface_Tension_Gamma'] * (self.settings_values['Skin_Compression_GammaC'] - self.settings_values['Surface_Tension_Gamma'])) / (self.Regenerated_Radius_N2[i] * self.settings_values['Skin_Compression_GammaC']))

            Initial_Allowable_Grad_He_Pa = ((2.0 * self.settings_values['Surface_Tension_Gamma'] * (self.settings_values['Skin_Compression_GammaC'] - self.settings_values['Surface_Tension_Gamma'])) / (self.Regenerated_Radius_He[i] * self.settings_values['Skin_Compression_GammaC']))

            self.Initial_Allowable_Gradient_N2[i] = (Initial_Allowable_Grad_N2_Pa / self.ATM) * self.Units_Factor
            self.Initial_Allowable_Gradient_He[i] = (Initial_Allowable_Grad_He_Pa / self.ATM) * self.Units_Factor

            self.Allowable_Gradient_He[i] = self.Initial_Allowable_Gradient_He[i]
            self.Allowable_Gradient_N2[i] = self.Initial_Allowable_Gradient_N2[i]

    def CALC_START_OF_DECO_ZONE(self, Starting_Depth, Rate):
        """
        Purpose: This subroutine uses the Bisection Method to find the depth at
        which the leading compartment just enters the decompression zone.
        Source:  "Numerical Recipes in Fortran 77", Cambridge University Press,
        1992.

        Side Effects: Sets
        `self.Depth_Start_of_Deco_Zone`

        or

        Raises a RootException

        Returns: None
        """

        # First initialize some variables
        self.Depth_Start_of_Deco_Zone = 0.0
        Starting_Ambient_Pressure = Starting_Depth + self.Barometric_Pressure

        Initial_Inspired_He_Pressure = (Starting_Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Helium[self.Mix_Number - 1]

        Initial_Inspired_N2_Pressure = (Starting_Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Nitrogen[self.Mix_Number - 1]

        Helium_Rate = Rate * self.Fraction_Helium[self.Mix_Number - 1]
        Nitrogen_Rate = Rate * self.Fraction_Nitrogen[self.Mix_Number - 1]

        # ESTABLISH THE BOUNDS FOR THE ROOT SEARCH USING THE BISECTION METHOD
        # AND CHECK TO MAKE SURE THAT THE ROOT WILL BE WITHIN BOUNDS.  PROCESS
        # EACH COMPARTMENT INDIVIDUALLY AND FIND THE MAXIMUM DEPTH ACROSS ALL
        # COMPARTMENTS (LEADING COMPARTMENT)
        # In this case, we are solving for time - the time when the gas tension in
        # the compartment will be equal to ambient pressure.  The low bound for time
        # is set at zero and the high bound is set at the time it would take to
        # ascend to zero ambient pressure (absolute).  Since the ascent rate is
        # negative, a multiplier of -1.0 is used to make the time positive.  The
        # desired point when gas tension equals ambient pressure is found at a time
        # somewhere between these endpoints.  The algorithm checks to make sure that
        # the solution lies in between these bounds by first computing the low bound
        # and high bound function values.

        Low_Bound = 0.0
        High_Bound = -1.0 * (Starting_Ambient_Pressure / Rate)

        for i in range(16):
            Initial_Helium_Pressure = self.Helium_Pressure[i]
            Initial_Nitrogen_Pressure = self.Nitrogen_Pressure[i]

            Function_at_Low_Bound = Initial_Helium_Pressure + Initial_Nitrogen_Pressure + self.Constant_Pressure_Other_Gases - Starting_Ambient_Pressure

            High_Bound_Helium_Pressure = SCHREINER_EQUATION(Initial_Inspired_He_Pressure, Helium_Rate, High_Bound, self.Helium_Time_Constant[i], Initial_Helium_Pressure)

            High_Bound_Nitrogen_Pressure = SCHREINER_EQUATION(Initial_Inspired_N2_Pressure, Nitrogen_Rate, High_Bound, self.Nitrogen_Time_Constant[i], Initial_Nitrogen_Pressure)

            Function_at_High_Bound = High_Bound_Helium_Pressure + High_Bound_Nitrogen_Pressure + self.Constant_Pressure_Other_Gases

            if (Function_at_High_Bound * Function_at_Low_Bound) >= 0.0:
                raise RootException("ERROR! ROOT IS NOT WITHIN BRACKETS")

            # APPLY THE BISECTION METHOD IN SEVERAL ITERATIONS UNTIL A SOLUTION WITH
            # THE DESIRED ACCURACY IS FOUND
            # Note: the program allows for up to 100 iterations.  Normally an exit will
            # be made from the loop well before that number.  If, for some reason, the
            # program exceeds 100 iterations, there will be a pause to alert the user.

            if Function_at_Low_Bound < 0.0:
                Time_to_Start_of_Deco_Zone = Low_Bound
                Differential_Change = High_Bound - Low_Bound
            else:
                Time_to_Start_of_Deco_Zone = High_Bound
                Differential_Change = Low_Bound - High_Bound

            for j in range(100):
                Last_Diff_Change = Differential_Change
                Differential_Change = Last_Diff_Change * 0.5

                Mid_Range_Time = Time_to_Start_of_Deco_Zone + Differential_Change

                Mid_Range_Helium_Pressure = SCHREINER_EQUATION(Initial_Inspired_He_Pressure, Helium_Rate, Mid_Range_Time, self.Helium_Time_Constant[i], Initial_Helium_Pressure)

                Mid_Range_Nitrogen_Pressure = SCHREINER_EQUATION(Initial_Inspired_N2_Pressure, Nitrogen_Rate, Mid_Range_Time, self.Nitrogen_Time_Constant[i], Initial_Nitrogen_Pressure)

                Function_at_Mid_Range = Mid_Range_Helium_Pressure + Mid_Range_Nitrogen_Pressure + self.Constant_Pressure_Other_Gases - (Starting_Ambient_Pressure + Rate * Mid_Range_Time)

                if Function_at_Mid_Range <= 0.0:
                    Time_to_Start_of_Deco_Zone = Mid_Range_Time

                if (abs(Differential_Change) < 1.0E-3) or (Function_at_Mid_Range == 0.0):
                    break

                if j == 100:
                    raise MaxIterationException('ERROR! ROOT SEARCH EXCEEDED MAXIMUM ITERATIONS')

            # When a solution with the desired accuracy is found, the program jumps out
            # of the loop to Line 170 and assigns the solution value for the individual
            # compartment.

            Cpt_Depth_Start_of_Deco_Zone = (Starting_Ambient_Pressure + Rate * Time_to_Start_of_Deco_Zone) - self.Barometric_Pressure
            # The overall solution will be the compartment with the maximum depth where
            # gas tension equals ambient pressure (leading compartment).

            self.Depth_Start_of_Deco_Zone = max(self.Depth_Start_of_Deco_Zone, Cpt_Depth_Start_of_Deco_Zone)

    def CALC_ASCENT_CEILING(self):
        """
        Purpose: This subprogram calculates the ascent ceiling (the safe ascent
        depth) in each compartment, based on the allowable gradients, and then
        finds the deepest ascent ceiling across all compartments.

        Side Effects: Sets
        `self.Ascent_Ceiling_Depth`

        Returns: None
        """

        Gas_Loading = 0.0
        Weighted_Allowable_Gradient = 0.0
        Tolerated_Ambient_Pressure = 0.0
        Compartment_Ascent_Ceiling = [0.0 for i in range(16)]

        # Since there are two sets of allowable gradients being tracked, one for
        # helium and one for nitrogen, a "weighted allowable gradient" must be
        # computed each time based on the proportions of helium and nitrogen in
        # each compartment.  This proportioning follows the methodology of
        # Buhlmann/Keller.  If there is no helium and nitrogen in the compartment,
        # such as after extended periods of oxygen breathing, then the minimum value
        # across both gases will be used.  It is important to note that if a
        # compartment is empty of helium and nitrogen, then the weighted allowable
        # gradient formula cannot be used since it will result in division by zero.

        for i in range(16):
            Gas_Loading = self.Helium_Pressure[i] + self.Nitrogen_Pressure[i]

            Weighted_Allowable_Gradient = None
            Tolerated_Ambient_Pressure = None

            if Gas_Loading > 0.0:
                Weighted_Allowable_Gradient = (self.Allowable_Gradient_He[i] * self.Helium_Pressure[i] + self.Allowable_Gradient_N2[i] * self.Nitrogen_Pressure[i]) / (self.Helium_Pressure[i] + self.Nitrogen_Pressure[i])
                Tolerated_Ambient_Pressure = (Gas_Loading + self.Constant_Pressure_Other_Gases) - Weighted_Allowable_Gradient

            else:
                Weighted_Allowable_Gradient = min(self.Allowable_Gradient_He[i], self.Allowable_Gradient_N2[i])
                Tolerated_Ambient_Pressure = self.Constant_Pressure_Other_Gases - Weighted_Allowable_Gradient

            #     The tolerated ambient pressure cannot be less than zero absolute, i.e.,
            #     the vacuum of outer space!
            if Tolerated_Ambient_Pressure < 0.0:
                Tolerated_Ambient_Pressure = 0.0

            #     The Ascent Ceiling Depth is computed in a loop after all of the individual
            #     compartment ascent ceilings have been calculated.  It is important that
            #     the Ascent Ceiling Depth (max ascent ceiling across all compartments) only
            #     be extracted from the compartment values and not be compared against some
            #     initialization value.  For example, if MAX(Ascent_Ceiling_Depth . .) was
            #     compared against zero, this could cause a program lockup because sometimes
            #     the Ascent Ceiling Depth needs to be negative (but not less than zero
            #     absolute ambient pressure) in order to decompress to the last stop at zero
            #     depth.

            Compartment_Ascent_Ceiling[i] = Tolerated_Ambient_Pressure - self.Barometric_Pressure

        self.Ascent_Ceiling_Depth = Compartment_Ascent_Ceiling[0]

        for i in range(1, 16):
            self.Ascent_Ceiling_Depth = max(self.Ascent_Ceiling_Depth, Compartment_Ascent_Ceiling[i])

    def PROJECTED_ASCENT(self, Starting_Depth, Rate, Step_Size):
        """
        Purpose: This subprogram performs a simulated ascent outside of the main
        program to ensure that a deco ceiling will not be violated due to unusual
        gas loading during ascent (on-gassing).  If the deco ceiling is violated,
        the stop depth will be adjusted deeper by the step size until a safe
        ascent can be made.

        Side Effects: Sets
        `self.Deco_Stop_Depth`

        Returns: None
        """

        New_Ambient_Pressure = self.Deco_Stop_Depth + self.Barometric_Pressure
        Starting_Ambient_Pressure = Starting_Depth + self.Barometric_Pressure

        Initial_Inspired_He_Pressure = (Starting_Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Helium[self.Mix_Number - 1]

        Initial_Inspired_N2_Pressure = (Starting_Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Nitrogen[self.Mix_Number - 1]

        Helium_Rate = Rate * self.Fraction_Helium[self.Mix_Number - 1]
        Nitrogen_Rate = Rate * self.Fraction_Nitrogen[self.Mix_Number - 1]

        Temp_Gas_Loading = [0.0 for i in range(16)]
        Allowable_Gas_Loading = [0.0 for i in range(16)]
        Initial_Helium_Pressure = [0.0 for i in range(16)]
        Initial_Nitrogen_Pressure = [0.0 for i in range(16)]

        for i in range(16):
            Initial_Helium_Pressure[i] = self.Helium_Pressure[i]
            Initial_Nitrogen_Pressure[i] = self.Nitrogen_Pressure[i]

        while(True):
            Ending_Ambient_Pressure = New_Ambient_Pressure

            Segment_Time = (Ending_Ambient_Pressure - Starting_Ambient_Pressure) / Rate

            for i in range(16):
                Temp_Helium_Pressure = SCHREINER_EQUATION(Initial_Inspired_He_Pressure, Helium_Rate, Segment_Time, self.Helium_Time_Constant[i], Initial_Helium_Pressure[i])

                Temp_Nitrogen_Pressure = SCHREINER_EQUATION(Initial_Inspired_N2_Pressure, Nitrogen_Rate, Segment_Time, self.Nitrogen_Time_Constant[i], Initial_Nitrogen_Pressure[i])

                Temp_Gas_Loading[i] = Temp_Helium_Pressure + Temp_Nitrogen_Pressure

                if Temp_Gas_Loading[i] > 0.0:
                    Weighted_Allowable_Gradient = (self.Allowable_Gradient_He[i] * Temp_Helium_Pressure + self.Allowable_Gradient_N2[i] * Temp_Nitrogen_Pressure) / Temp_Gas_Loading[i]
                else:
                    Weighted_Allowable_Gradient = min(self.Allowable_Gradient_He[i], self.Allowable_Gradient_N2[i])

                Allowable_Gas_Loading[i] = Ending_Ambient_Pressure + Weighted_Allowable_Gradient - self.Constant_Pressure_Other_Gases

            end_sub = True
            for j in range(16):
                if Temp_Gas_Loading[j] > Allowable_Gas_Loading[j]:
                    New_Ambient_Pressure = Ending_Ambient_Pressure + Step_Size
                    self.Deco_Stop_Depth = self.Deco_Stop_Depth + Step_Size
                    end_sub = False

                    break

            if not end_sub:
                continue
            else:
                break

    def _calculate_deco_gradient(self, Allowable_Gradient_Molecule, Amb_Press_First_Stop_Pascals, Amb_Press_Next_Stop_Pascals):
        """Calculates the decompression gradient for Boyles_Law_Compensation.

        Side Effects: None

        Returns: A floating point value
        """

        Allow_Grad_First_Stop_Pa = (Allowable_Gradient_Molecule / self.Units_Factor) * self.ATM
        Radius_First_Stop = (2.0 * self.Surface_Tension_Gamma) / Allow_Grad_First_Stop_Pa

        A = Amb_Press_Next_Stop_Pascals
        B = -2.0 * self.Surface_Tension_Gamma
        C = (Amb_Press_First_Stop_Pascals + (2.0 * self.Surface_Tension_Gamma) / Radius_First_Stop) * Radius_First_Stop * (Radius_First_Stop * (Radius_First_Stop))

        Low_Bound = Radius_First_Stop
        High_Bound = Radius_First_Stop * (Amb_Press_First_Stop_Pascals / Amb_Press_Next_Stop_Pascals) ** (1.0 / 3.0)

        Ending_Radius = RADIUS_ROOT_FINDER(A, B, C, Low_Bound, High_Bound)

        Deco_Gradient_Pascals = (2.0 * self.Surface_Tension_Gamma) / Ending_Radius
        return (Deco_Gradient_Pascals / self.ATM) * self.Units_Factor

    def BOYLES_LAW_COMPENSATION(self, First_Stop_Depth, Deco_Stop_Depth, Step_Size):
        """
        Purpose: This subprogram calculates the reduction in allowable gradients
        with decreasing ambient pressure during the decompression profile based
        on Boyle's Law considerations.

        Side Effects: Sets
        `self.Deco_Gradient_He`,
        `self.Deco_Gradient_N2`

        Returns: None
        """

        Next_Stop = Deco_Stop_Depth - Step_Size
        Ambient_Pressure_First_Stop = First_Stop_Depth + self.Barometric_Pressure
        Ambient_Pressure_Next_Stop = Next_Stop + self.Barometric_Pressure

        Amb_Press_First_Stop_Pascals = (Ambient_Pressure_First_Stop / self.Units_Factor) * self.ATM
        Amb_Press_Next_Stop_Pascals = (Ambient_Pressure_Next_Stop / self.Units_Factor) * self.ATM

        for i in range(16):
            # Helium Calculation
            self.Deco_Gradient_He[i] = self._calculate_deco_gradient(self.Allowable_Gradient_He[i], Amb_Press_First_Stop_Pascals, Amb_Press_Next_Stop_Pascals)

            # Nitrogen Calculation
            self.Deco_Gradient_N2[i] = self._calculate_deco_gradient(self.Allowable_Gradient_N2[i], Amb_Press_First_Stop_Pascals, Amb_Press_Next_Stop_Pascals)

    def DECOMPRESSION_STOP(self, Deco_Stop_Depth, Step_Size):
        """
        Purpose: This subprogram calculates the required time at each
        decompression stop.

        Side Effects: Sets
        `self.Ending_Ambient_Pressure`,
        `self.Helium_Pressure`,
        `self.Nitrogen_Pressure`
        `self.Run_Time`,
        `self.Segment_Number`,
        `self.Segment_Time`,

        or

        Raises an OffGassingException

        Returns: None
        """

        Deco_Ceiling_Depth = 0.0
        Last_Run_Time = self.Run_Time
        Round_Up_Operation = round((Last_Run_Time / self.Minimum_Deco_Stop_Time) + 0.5) * self.Minimum_Deco_Stop_Time
        self.Segment_Time = Round_Up_Operation - self.Run_Time
        self.Run_Time = Round_Up_Operation
        Temp_Segment_Time = self.Segment_Time
        self.Segment_Number = self.Segment_Number + 1
        Ambient_Pressure = Deco_Stop_Depth + self.Barometric_Pressure
        self.Ending_Ambient_Pressure = Ambient_Pressure
        Next_Stop = Deco_Stop_Depth - Step_Size

        Inspired_Helium_Pressure = (Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Helium[self.Mix_Number - 1]

        Inspired_Nitrogen_Pressure = (Ambient_Pressure - self.Water_Vapor_Pressure) * self.Fraction_Nitrogen[self.Mix_Number - 1]

        # Check to make sure that program won't lock up if unable to decompress
        # to the next stop.  If so, write error message and terminate program.

        for i in range(16):
            if(Inspired_Helium_Pressure + Inspired_Nitrogen_Pressure) > 0.0:
                Weighted_Allowable_Gradient = (self.Deco_Gradient_He[i] * Inspired_Helium_Pressure + self.Deco_Gradient_N2[i] * Inspired_Nitrogen_Pressure) / (Inspired_Helium_Pressure + Inspired_Nitrogen_Pressure)

                if (Inspired_Helium_Pressure + Inspired_Nitrogen_Pressure + self.Constant_Pressure_Other_Gases - Weighted_Allowable_Gradient) > (Next_Stop + self.Barometric_Pressure):
                    raise OffGassingException("ERROR! OFF-GASSING GRADIENT IS TOO SMALL TO DECOMPRESS AT THE %f STOP. Next stop: %f" % (Deco_Stop_Depth, Next_Stop))

        while(True):
            for i in range(16):

                Initial_Helium_Pressure = self.Helium_Pressure[i]
                Initial_Nitrogen_Pressure = self.Nitrogen_Pressure[i]

                self.Helium_Pressure[i] = HALDANE_EQUATION(Initial_Helium_Pressure, Inspired_Helium_Pressure, self.Helium_Time_Constant[i], self.Segment_Time)

                self.Nitrogen_Pressure[i] = HALDANE_EQUATION(Initial_Nitrogen_Pressure, Inspired_Nitrogen_Pressure, self.Nitrogen_Time_Constant[i], self.Segment_Time)

            Deco_Ceiling_Depth = self.CALC_DECO_CEILING()
            if Deco_Ceiling_Depth > Next_Stop:
                self.Segment_Time = self.Minimum_Deco_Stop_Time
                Time_Counter = Temp_Segment_Time
                Temp_Segment_Time = Time_Counter + self.Minimum_Deco_Stop_Time
                Last_Run_Time = self.Run_Time
                self.Run_Time = Last_Run_Time + self.Minimum_Deco_Stop_Time
                continue
            break

        self.Segment_Time = Temp_Segment_Time

    def CALC_DECO_CEILING(self):
        """
        Purpose: This subprogram calculates the deco ceiling (the safe ascent
        depth) in each compartment, based on the allowable "deco gradients"
        computed in the Boyle's Law Compensation subroutine, and then finds the
        deepest deco ceiling across all compartments.  This deepest value
        (Deco Ceiling Depth) is then used by the Decompression Stop subroutine
        to determine the actual deco schedule.

        Side Effects: None

        Returns: `self.Deco_Ceiling_Depth`
        """

        # Since there are two sets of deco gradients being tracked, one for
        # helium and one for nitrogen, a "weighted allowable gradient" must be
        # computed each time based on the proportions of helium and nitrogen in
        # each compartment.  This proportioning follows the methodology of
        # Buhlmann/Keller.  If there is no helium and nitrogen in the compartment,
        # such as after extended periods of oxygen breathing, then the minimum value
        # across both gases will be used.  It is important to note that if a
        # compartment is empty of helium and nitrogen, then the weighted allowable
        # gradient formula cannot be used since it will result in division by zero.

        Compartment_Deco_Ceiling = [0.0 for i in range(16)]

        for i in range(16):
            Gas_Loading = self.Helium_Pressure[i] + self.Nitrogen_Pressure[i]

            if Gas_Loading > 0.0:
                Weighted_Allowable_Gradient = (self.Deco_Gradient_He[i] * self.Helium_Pressure[i] + self.Deco_Gradient_N2[i] * self.Nitrogen_Pressure[i]) / (self.Helium_Pressure[i] + self.Nitrogen_Pressure[i])

                Tolerated_Ambient_Pressure = (Gas_Loading + self.Constant_Pressure_Other_Gases) - Weighted_Allowable_Gradient
            else:
                Weighted_Allowable_Gradient = min(self.Deco_Gradient_He[i], self.Deco_Gradient_N2[i])
                Tolerated_Ambient_Pressure = self.Constant_Pressure_Other_Gases - Weighted_Allowable_Gradient

            # The tolerated ambient pressure cannot be less than zero absolute, i.e.,
            # the vacuum of outer space!
            if Tolerated_Ambient_Pressure < 0.0:
                Tolerated_Ambient_Pressure = 0.0

            # The Deco Ceiling Depth is computed in a loop after all of the individual
            # compartment deco ceilings have been calculated.  It is important that the
            # Deco Ceiling Depth (max deco ceiling across all compartments) only be
            # extracted from the compartment values and not be compared against some
            # initialization value.  For example, if MAX(Deco_Ceiling_Depth . .) was
            # compared against zero, this could cause a program lockup because sometimes
            # the Deco Ceiling Depth needs to be negative (but not less than absolute
            # zero) in order to decompress to the last stop at zero depth.

            Compartment_Deco_Ceiling[i] = Tolerated_Ambient_Pressure - self.Barometric_Pressure

        Deco_Ceiling_Depth = Compartment_Deco_Ceiling[0]
        # can replace these for loops with just max(Compartment_Deco_Ceiling)
        for i in range(1, 16):
            Deco_Ceiling_Depth = max(Deco_Ceiling_Depth, Compartment_Deco_Ceiling[i])

        return Deco_Ceiling_Depth

    def CRITICAL_VOLUME(self, Deco_Phase_Volume_Time):
        """
        Purpose: This subprogram applies the VPM Critical Volume Algorithm.  This
        algorithm will compute "relaxed" gradients for helium and nitrogen based
        on the setting of the Critical Volume Parameter Lambda.

        Side Effects: Sets
        `self.Allowable_Gradient_He`,
        `self.Allowable_Gradient_N2`

        Returns: None
        """

        # Note:  Since the Critical Volume Parameter Lambda was defined in units of
        # fsw-min in the original papers by Yount and colleauges, the same
        # convention is retained here.  Although Lambda is adjustable only in units
        # of fsw-min in the program settings (range from 6500 to 8300 with default
        # 7500), it will convert to the proper value in Pascals-min in this
        # subroutine regardless of which diving pressure units are being used in
        # the main program - feet of seawater (fsw) or meters of seawater (msw).
        # The allowable gradient is computed using the quadratic formula (refer to
        # separate write-up posted on the Deco List web site).

        Parameter_Lambda_Pascals = (self.Crit_Volume_Parameter_Lambda / 33.0) * self.ATM
        Phase_Volume_Time = [0.0 for i in range(16)]

        for i in range(16):

            Phase_Volume_Time = Deco_Phase_Volume_Time + self.Surface_Phase_Volume_Time[i]

            # Helium Calculations
            Adj_Crush_Pressure_He_Pascals = (self.Adjusted_Crushing_Pressure_He[i] / self.Units_Factor) * self.ATM
            Initial_Allowable_Grad_He_Pa = (self.Initial_Allowable_Gradient_He[i] / self.Units_Factor) * self.ATM

            B = Initial_Allowable_Grad_He_Pa + (Parameter_Lambda_Pascals * self.Surface_Tension_Gamma) / (self.Skin_Compression_GammaC * Phase_Volume_Time)

            C = (self.Surface_Tension_Gamma * (self.Surface_Tension_Gamma * (Parameter_Lambda_Pascals * Adj_Crush_Pressure_He_Pascals))) / (self.Skin_Compression_GammaC * (self.Skin_Compression_GammaC * Phase_Volume_Time))

            New_Allowable_Grad_He_Pascals = (B + sqrt(B ** 2 - 4.0 * C)) / 2.0

            self.Allowable_Gradient_He[i] = (New_Allowable_Grad_He_Pascals / self.ATM) * self.Units_Factor

            # Nitrogen Calculations
            Adj_Crush_Pressure_N2_Pascals = (self.Adjusted_Crushing_Pressure_N2[i] / self.Units_Factor) * self.ATM

            Initial_Allowable_Grad_N2_Pa = (self.Initial_Allowable_Gradient_N2[i] / self.Units_Factor) * self.ATM

            B = Initial_Allowable_Grad_N2_Pa + (Parameter_Lambda_Pascals * self.Surface_Tension_Gamma) / (self.Skin_Compression_GammaC * Phase_Volume_Time)

            C = (self.Surface_Tension_Gamma * (self.Surface_Tension_Gamma * (Parameter_Lambda_Pascals * Adj_Crush_Pressure_N2_Pascals))) / (self.Skin_Compression_GammaC * (self.Skin_Compression_GammaC * Phase_Volume_Time))

            New_Allowable_Grad_N2_Pascals = (B + sqrt(B ** 2 - 4.0 * C)) / 2.0

            self.Allowable_Gradient_N2[i] = (New_Allowable_Grad_N2_Pascals / self.ATM) * self.Units_Factor

    def CALC_SURFACE_PHASE_VOLUME_TIME(self):
        """
        Purpose: This subprogram computes the surface portion of the total phase
        volume time.  This is the time factored out of the integration of
        supersaturation gradient x time over the surface interval.  The VPM
        considers the gradients that allow bubbles to form or to drive bubble
        growth both in the water and on the surface after the dive.

        This subroutine is a new development to the VPM algorithm in that it
        computes the time course of supersaturation gradients on the surface
        when both helium and nitrogen are present.  Refer to separate write-up
        for a more detailed explanation of this algorithm.

        Side Effects: Sets
        `self.Surface_Phase_Volume_Time`

        Returns: None
        """

        Surface_Inspired_N2_Pressure = (self.Barometric_Pressure - self.Water_Vapor_Pressure) * self.fraction_inert_gas
        for i in range(16):
            if self.Nitrogen_Pressure[i] > Surface_Inspired_N2_Pressure:

                self.Surface_Phase_Volume_Time[i] = (self.Helium_Pressure[i] / self.Helium_Time_Constant[i] + (self.Nitrogen_Pressure[i] - Surface_Inspired_N2_Pressure) / self.Nitrogen_Time_Constant[i]) / (self.Helium_Pressure[i] + self.Nitrogen_Pressure[i] - Surface_Inspired_N2_Pressure)

            elif (self.Nitrogen_Pressure[i] <= Surface_Inspired_N2_Pressure) and (self.Helium_Pressure[i] + self.Nitrogen_Pressure[i] >= Surface_Inspired_N2_Pressure):
                Decay_Time_to_Zero_Gradient = 1.0 / (self.Nitrogen_Time_Constant[i] - self.Helium_Time_Constant[i]) * log((Surface_Inspired_N2_Pressure - self.Nitrogen_Pressure[i]) / self.Helium_Pressure[i])

                Integral_Gradient_x_Time = self.Helium_Pressure[i] / self.Helium_Time_Constant[i] * (1.0 - exp(-self.Helium_Time_Constant[i] * Decay_Time_to_Zero_Gradient)) + (self.Nitrogen_Pressure[i] - Surface_Inspired_N2_Pressure) / self.Nitrogen_Time_Constant[i] * (1.0 - exp(-self.Nitrogen_Time_Constant[i] * Decay_Time_to_Zero_Gradient))

                self.Surface_Phase_Volume_Time[i] = Integral_Gradient_x_Time / (self.Helium_Pressure[i] + self.Nitrogen_Pressure[i] - Surface_Inspired_N2_Pressure)

            else:
                self.Surface_Phase_Volume_Time[i] = 0.0

    def validate_data(self):
        """
        Purpose: Check the the data loaded from the input file is valid

        Side Effects: Sets
        `self.Critical_Radius_He_Microns`
        `self.Critical_Radius_N2_Microns`,
        `self.Diver_Acclimatized`,
        `self.units_fsw`,

        or

        Raises InputFileException, AltitudeException, ValueError

        Returns: None
        """
        if self.settings_values["Units"].lower() == "fsw":
            self.units_fsw = True
        elif self.settings_values["Units"].lower() == "msw":
            self.units_fsw = False
        else:
            raise ValueError("Bad Unit of measurement: Units = %s, must be 'fsw' or 'msw'" % self.settings_values["Units"])

        if self.settings_values["Regeneration_Time_Constant"] <= 0:
            raise InputFileException("Regeneration_Time_Constant must be greater than 0")

        if self.units_fsw and self.altitude_values['Altitude_of_Dive'] > 30000.0:
            raise AltitudeException("ERROR! ALTITUDE OF DIVE HIGHER THAN MOUNT EVEREST")

        if (not self.units_fsw) and self.altitude_values['Altitude_of_Dive'] > 9144.0:
            raise AltitudeException("ERROR! ALTITUDE OF DIVE HIGHER THAN MOUNT EVEREST")

        if self.altitude_values['Diver_Acclimatized_at_Altitude'].lower() == 'yes':
            self.Diver_Acclimatized = True
        elif self.altitude_values['Diver_Acclimatized_at_Altitude'].lower() == 'no':
            self.Diver_Acclimatized = False
        else:
            raise AltitudeException("ERROR! DIVER ACCLIMATIZED AT ALTITUDE MUST BE YES OR NO")

        self.Critical_Radius_N2_Microns = self.settings_values['Critical_Radius_N2_Microns']
        self.Critical_Radius_He_Microns = self.settings_values['Critical_Radius_He_Microns']

        # nitrogen
        if self.settings_values['Critical_Radius_N2_Microns'] < 0.2 or self.settings_values['Critical_Radius_N2_Microns'] > 1.35:
            raise ValueError("Bad Critical Radius N2 Microns: Critical_Radius_N2_Microns = %f, must be between '0.2' and '1.35'" % self.settings_values["Critical_Radius_N2_Microns"])

        # helium
        if self.settings_values['Critical_Radius_He_Microns'] < 0.2 or self.settings_values['Critical_Radius_He_Microns'] > 1.35:
            raise ValueError("Bad Critical_Radius_He_Microns: Critical_Radius_He_Microns = %f, must be between '0.2' and '1.35'" % self.settings_values["Critical_Radius_He_Microns"])

    def initialize_data(self):
        """
        Purpose: Initialize the object with the data loaded from the input file

        Side Effects: Sets

        `self.Adjusted_Critical_Radius_He`,
        `self.Adjusted_Critical_Radius_N2`,
        `self.Altitude_Dive_Algorithm_Off`,
        `self.Altitude_of_Dive`,
        `self.Amb_Pressure_Onset_of_Imperm`,
        `self.Barometric_Pressure`,
        `self.Constant_Pressure_Other_Gases`,
        `self.Crit_Volume_Parameter_Lambda`,
        `self.Critical_Radius_He_Microns`,
        `self.Critical_Radius_N2_Microns`,
        `self.Critical_Volume_Algorithm_Off`,
        `self.Gas_Tension_Onset_of_Imperm`,
        `self.Gradient_Onset_of_Imperm_Atm`,
        `self.Helium_Pressure`,
        `self.Helium_Time_Constant`,
        `self.Initial_Critical_Radius_He`,
        `self.Initial_Critical_Radius_N2`,
        `self.Max_Actual_Gradient`,
        `self.Max_Crushing_Pressure_He`,
        `self.Max_Crushing_Pressure_N2`,
        `self.Minimum_Deco_Stop_Time`,
        `self.Minimum_Deco_Stop_Time`,
        `self.Nitrogen_Pressure`,
        `self.Nitrogen_Time_Constant`,
        `self.Pressure_Other_Gases_mmHg`,
        `self.Regeneration_Time_Constant`,
        `self.Run_Time`,
        `self.Segment_Number`,
        `self.Skin_Compression_GammaC`,
        `self.Surface_Phase_Volume_Time`,
        `self.Surface_Tension_Gamma`,
        `self.Surface_Tension_Gamma`,
        `self.Units_Factor`,
        `self.Units_Word1`,
        `self.Units_Word2`,
        `self.Water_Vapor_Pressure`,

        or

        Raises AltitudeException, ValueError

        Returns: None
        """
        self.Surface_Tension_Gamma = self.settings_values['Surface_Tension_Gamma']
        self.Skin_Compression_GammaC = self.settings_values['Skin_Compression_GammaC']
        self.Crit_Volume_Parameter_Lambda = self.settings_values['Crit_Volume_Parameter_Lambda']
        self.Gradient_Onset_of_Imperm_Atm = self.settings_values['Gradient_Onset_of_Imperm_Atm']

        self.Minimum_Deco_Stop_Time = self.settings_values['Minimum_Deco_Stop_Time']
        self.Critical_Radius_N2_Microns = self.settings_values['Critical_Radius_N2_Microns']
        self.Critical_Radius_He_Microns = self.settings_values['Critical_Radius_He_Microns']
        self.Regeneration_Time_Constant = self.settings_values['Regeneration_Time_Constant']

        self.Surface_Tension_Gamma = self.settings_values['Surface_Tension_Gamma']
        self.Minimum_Deco_Stop_Time = self.settings_values['Minimum_Deco_Stop_Time']

        # INITIALIZE CONSTANTS/VARIABLES BASED ON SELECTION OF UNITS - FSW OR MSW
        # fsw = feet of seawater, a unit of pressure
        # msw = meters of seawater, a unit of pressure

        if self.units_fsw:
            self.Units_Word1 = "fswg"
            self.Units_Word2 = "fsw/min"
            self.Units_Factor = 33.0
            self.Water_Vapor_Pressure = 1.607  # based on respiratory quotient of 0.8 (Schreiner value)

        else:
            self.Units_Word1 = "mswg"
            self.Units_Word2 = "msw/min"
            self.Units_Factor = 10.1325
            self.Water_Vapor_Pressure = 0.493     # based on respiratory quotient of 0.8 (Schreiner value)

        # INITIALIZE CONSTANTS/VARIABLES
        self.Constant_Pressure_Other_Gases = (self.settings_values["Pressure_Other_Gases_mmHg"] / 760.0) * self.Units_Factor
        self.Run_Time = 0.0
        self.Segment_Number = 0

        for i in range(16):
            self.Helium_Time_Constant[i] = log(2.0) / self.Helium_Half_Time[i]
            self.Nitrogen_Time_Constant[i] = log(2.0) / self.Nitrogen_Half_Time[i]
            self.Max_Crushing_Pressure_He[i] = 0.0
            self.Max_Crushing_Pressure_N2[i] = 0.0
            self.Max_Actual_Gradient[i] = 0.0
            self.Surface_Phase_Volume_Time[i] = 0.0
            self.Amb_Pressure_Onset_of_Imperm[i] = 0.0
            self.Gas_Tension_Onset_of_Imperm[i] = 0.0
            self.Initial_Critical_Radius_N2[i] = self.settings_values["Critical_Radius_N2_Microns"] * 1.0E-6
            self.Initial_Critical_Radius_He[i] = self.settings_values["Critical_Radius_He_Microns"] * 1.0E-6

        if self.settings_values['Critical_Volume_Algorithm'].lower() == "on":
            self.Critical_Volume_Algorithm_Off = False
        elif self.settings_values['Critical_Volume_Algorithm'].lower() == "off":
            self.Critical_Volume_Algorithm_Off = True
        else:
            raise ValueError("Bad Critical Volume Algorithm: Critical_Volume_Algorithm = %s, must be 'OFF or 'ON''" % self.settings_values["Critical_Volume_Algorithm"])

        if self.settings_values["Altitude_Dive_Algorithm"].lower() == "on":
            self.Altitude_Dive_Algorithm_Off = False
            if self.altitude_values["Ascent_to_Altitude_Hours"] <= 0 and self.Diver_Acclimatized is False:
                raise AltitudeException("If diver is not acclimatized, Ascent_to_Altitude_Time must be greater than 0")

        elif self.settings_values["Altitude_Dive_Algorithm"].lower() == "off":
            self.Altitude_Dive_Algorithm_Off = True
        else:
            raise ValueError("Bad Altitude Dive Algorithm: Altitude_Dive_Algorithm = %s, must be 'OFF or 'ON''" % self.settings_values["Altitude_Dive_Algorithm"])

        #     INITIALIZE VARIABLES FOR SEA LEVEL OR ALTITUDE DIVE
        #     See subroutines for explanation of altitude calculations.  Purposes are
        #     1) to determine barometric pressure and 2) set or adjust the VPM critical
        #     radius variables and gas loadings, as applicable, based on altitude,
        #     ascent to altitude before the dive, and time at altitude before the dive

        if self.Altitude_Dive_Algorithm_Off:
            self.Altitude_of_Dive = 0.0
            self.Barometric_Pressure = CALC_BAROMETRIC_PRESSURE(self.Altitude_of_Dive, self.units_fsw)

            for i in range(16):
                self.Adjusted_Critical_Radius_N2[i] = self.Initial_Critical_Radius_N2[i]
                self.Adjusted_Critical_Radius_He[i] = self.Initial_Critical_Radius_He[i]
                self.Helium_Pressure[i] = 0.0
                self.Nitrogen_Pressure[i] = (self.Barometric_Pressure - self.Water_Vapor_Pressure) * self.fraction_inert_gas
        else:
            self.VPM_ALTITUDE_DIVE_ALGORITHM(self.altitude_values)

    def set_gas_mixes(self, dive):
        """
        Purpose: Checks the given gas mix fractions add up to 1.0, and addes them
        to the object

        Side Effects: Sets
        `self.Fraction_Helium`,
        `self.Fraction_Nitrogen`

        or

        Raises an InputFileException

        Returns: None
        """
        num_gas_mixes = dive["num_gas_mixes"]
        Fraction_Oxygen = [0.0 for i in range(num_gas_mixes)]
        self.Fraction_Helium = [0.0 for i in range(num_gas_mixes)]
        self.Fraction_Nitrogen = [0.0 for i in range(num_gas_mixes)]

        for i in range(num_gas_mixes):
            gasmix_summary = dive["gasmix_summary"]
            Fraction_Oxygen[i] = gasmix_summary[i]["fraction_O2"]
            self.Fraction_Nitrogen[i] = gasmix_summary[i]["fraction_N2"]
            self.Fraction_Helium[i] = gasmix_summary[i]["fraction_He"]
            Sum_of_Fractions = Fraction_Oxygen[i] + self.Fraction_Nitrogen[i] + self.Fraction_Helium[i]

            if Sum_of_Fractions != 1.0:
                raise InputFileException("ERROR IN INPUT FILE (gas mixes don't add up to 1.0")

        for i in range(num_gas_mixes):
            self.output_object.add_gasmix(Fraction_Oxygen[i], self.Fraction_Nitrogen[i], self.Fraction_Helium[i])

    def profile_code_loop(self, dive):
        """
        Purpose:
        PROCESS DIVE AS A SERIES OF ASCENT/DESCENT AND CONSTANT DEPTH
        SEGMENTS. THIS ALLOWS FOR MULTI-LEVEL DIVES AND UNUSUAL PROFILES. UPDATE
        GAS LOADINGS FOR EACH SEGMENT.  IF IT IS A DESCENT SEGMENT, CALC CRUSHING
        PRESSURE ON CRITICAL RADII IN EACH COMPARTMENT.
        "Instantaneous" descents are not used in the VPM.  All ascent/descent
        segments must have a realistic rate of ascent/descent.  Unlike Haldanian
        models, the VPM is actually more conservative when the descent rate is
        slower because the effective crushing pressure is reduced.  Also, a
        realistic actual supersaturation gradient must be calculated during
        ascents as this affects critical radii adjustments for repetitive dives.

        Profile codes: 1 = Ascent/Descent, 2 = Constant Depth, 99 = Decompress

        Side Effects: Sets
        `self.Depth`,
        `self.Ending_Depth`,
        `self.Mix_Number`,
        `self.Rate`,
        `self.Run_Time_End_of_Segment`
        `self.Starting_Depth`,

        or

        Raises an InputFileException

        Returns: None
        """
        for profile in dive["profile_codes"]:
            Profile_Code = profile["profile_code"]

            if Profile_Code == 1:
                self.Starting_Depth = profile["starting_depth"]
                self.Ending_Depth = profile["ending_depth"]
                self.Rate = profile["rate"]
                self.Mix_Number = profile["gasmix"]

                self.GAS_LOADINGS_ASCENT_DESCENT(self.Starting_Depth, self.Ending_Depth, self.Rate)
                if self.Ending_Depth > self.Starting_Depth:
                    self.CALC_CRUSHING_PRESSURE(self.Starting_Depth, self.Ending_Depth, self.Rate)

                # the error seems unnecessary
                Word = ""
                if self.Ending_Depth > self.Starting_Depth:
                    Word = 'Descent'
                elif self.Starting_Depth > self.Ending_Depth:
                    Word = 'Ascent '
                else:
                    Word = 'ERROR'

                self.output_object.add_dive_profile_entry_descent(self.Segment_Number, self.Segment_Time, self.Run_Time, self.Mix_Number, Word, self.Starting_Depth, self.Ending_Depth, self.Rate)

            elif Profile_Code == 2:
                self.Depth = profile["depth"]
                self.Run_Time_End_of_Segment = profile["run_time_at_end_of_segment"]
                self.Mix_Number = profile["gasmix"]

                self.GAS_LOADINGS_CONSTANT_DEPTH(self.Depth, self.Run_Time_End_of_Segment)

                self.output_object.add_dive_profile_entry_ascent(self.Segment_Number, self.Segment_Time, self.Run_Time, self.Mix_Number, self.Depth)
            elif Profile_Code == 99:
                break
            else:
                raise InputFileException("Invalid profile code %d. Valid profile codes are 1 (descent), 2 (constant), and 99 (ascent)" % (Profile_Code))

    def deco_stop_loop_block_within_critical_volume_loop(self):
        """
        Purpose:
        DECO STOP LOOP BLOCK WITHIN CRITICAL VOLUME LOOP
        This loop computes a decompression schedule to the surface during each
        iteration of the critical volume loop.  No output is written from this
        loop, rather it computes a schedule from which the in-water portion of the
        total phase volume time (Deco_Phase_Volume_Time) can be extracted.  Also,
        the gas loadings computed at the end of this loop are used in the subroutine
        which computes the out-of-water portion of the total phase volume time
        (Surface_Phase_Volume_Time) for that schedule.

        Note that exit is made from the loop after last ascent is made to a deco
        stop depth that is less than or equal to zero.  A final deco stop less
        than zero can happen when the user makes an odd step size change during
        ascent - such as specifying a 5 msw step size change at the 3 msw stop!

        Side Effects: Sets
        `self.Deco_Stop_Depth`,
        `self.Last_Run_Time`
        `self.Next_Stop`,
        `self.Starting_Depth`,

        Returns: None
        """
        while(True):
            self.GAS_LOADINGS_ASCENT_DESCENT(self.Starting_Depth, self.Deco_Stop_Depth, self.Rate)

            if self.Deco_Stop_Depth <= 0.0:
                break

            if self.Number_of_Changes > 1:
                for i in range(1, self.Number_of_Changes):
                    if self.Depth_Change[i] >= self.Deco_Stop_Depth:
                        self.Mix_Number = self.Mix_Change[i]
                        self.Rate = self.Rate_Change[i]
                        self.Step_Size = self.Step_Size_Change[i]

            self.BOYLES_LAW_COMPENSATION(self.First_Stop_Depth, self.Deco_Stop_Depth, self.Step_Size)

            self.DECOMPRESSION_STOP(self.Deco_Stop_Depth, self.Step_Size)

            self.Starting_Depth = self.Deco_Stop_Depth
            self.Next_Stop = self.Deco_Stop_Depth - self.Step_Size
            self.Deco_Stop_Depth = self.Next_Stop
            self.Last_Run_Time = self.Run_Time

    def critical_volume_decision_tree(self):
        """Purpose:

        Side Effects: Sets

        `self.Deco_Stop_Depth`,
        `self.Helium_Pressure`,
        `self.Last_Run_Time`,
        `self.Mix_Number`,
        `self.Next_Stop`
        `self.Nitrogen_Pressure`,
        `self.Rate`,
        `self.Run_Time`,
        `self.Segment_Number`,
        `self.Starting_Depth`,
        `self.Step_Size`,
        `self.Stop_Time`,

        Returns: None
        """
        for i in range(16):
            self.Helium_Pressure[i] = self.He_Pressure_Start_of_Ascent[i]
            self.Nitrogen_Pressure[i] = self.N2_Pressure_Start_of_Ascent[i]

        self.Run_Time = self.Run_Time_Start_of_Ascent
        self.Segment_Number = self.Segment_Number_Start_of_Ascent
        self.Starting_Depth = self.Depth_Change[0]
        self.Mix_Number = self.Mix_Change[0]
        self.Rate = self.Rate_Change[0]
        self.Step_Size = self.Step_Size_Change[0]
        self.Deco_Stop_Depth = self.First_Stop_Depth
        self.Last_Run_Time = 0.0

        # DECO STOP LOOP BLOCK FOR FINAL DECOMPRESSION SCHEDULE

        while(True):
            self.GAS_LOADINGS_ASCENT_DESCENT(self.Starting_Depth, self.Deco_Stop_Depth, self.Rate)
            # DURING FINAL DECOMPRESSION SCHEDULE PROCESS, COMPUTE MAXIMUM ACTUAL
            # SUPERSATURATION GRADIENT RESULTING IN EACH COMPARTMENT
            # If there is a repetitive dive, this will be used later in the VPM
            # Repetitive Algorithm to adjust the values for critical radii.
            self.CALC_MAX_ACTUAL_GRADIENT(self.Deco_Stop_Depth)

            self.output_object.add_decompression_profile_ascent(self.Segment_Number, self.Segment_Time, self.Run_Time, self.Mix_Number, self.Deco_Stop_Depth, self.Rate)
            if self.Deco_Stop_Depth <= 0.0:
                break

            if self.Number_of_Changes > 1:
                for i in range(1, self.Number_of_Changes):
                    if self.Depth_Change[i] >= self.Deco_Stop_Depth:
                        self.Mix_Number = self.Mix_Change[i]
                        self.Rate = self.Rate_Change[i]
                        self.Step_Size = self.Step_Size_Change[i]

            self.BOYLES_LAW_COMPENSATION(self.First_Stop_Depth, self.Deco_Stop_Depth, self.Step_Size)
            self.DECOMPRESSION_STOP(self.Deco_Stop_Depth, self.Step_Size)
            # This next bit justs rounds up the stop time at the first stop to be in
            # whole increments of the minimum stop time (to make for a nice deco table).

            if self.Last_Run_Time == 0.0:
                self.Stop_Time = round((self.Segment_Time / self.Minimum_Deco_Stop_Time) + 0.5) * self.Minimum_Deco_Stop_Time
            else:
                self.Stop_Time = self.Run_Time - self.Last_Run_Time

            # DURING FINAL DECOMPRESSION SCHEDULE, IF MINIMUM STOP TIME PARAMETER IS A
            # WHOLE NUMBER (i.e. 1 minute) THEN WRITE DECO SCHEDULE USING INTEGER
            # NUMBERS (looks nicer).  OTHERWISE, USE DECIMAL NUMBERS.
            # Note: per the request of a noted exploration diver(!), program now allows
            # a minimum stop time of less than one minute so that total ascent time can
            # be minimized on very long dives.  In fact, with step size set at 1 fsw or
            # 0.2 msw and minimum stop time set at 0.1 minute (6 seconds), a near
            # continuous decompression schedule can be computed.

            if trunc(self.Minimum_Deco_Stop_Time) == self.Minimum_Deco_Stop_Time:
                self.output_object.add_decompression_profile_constant(self.Segment_Number, self.Segment_Time, self.Run_Time, self.Mix_Number, int(self.Deco_Stop_Depth), int(self.Stop_Time))
            else:
                self.output_object.add_decompression_profile_constant(self.Segment_Number, self.Segment_Time, self.Run_Time, self.Mix_Number, self.Deco_Stop_Depth, self.Stop_Time)

            self.Starting_Depth = self.Deco_Stop_Depth
            self.Next_Stop = self.Deco_Stop_Depth - self.Step_Size
            self.Deco_Stop_Depth = self.Next_Stop
            self.Last_Run_Time = self.Run_Time

    def critical_volume_loop(self):
        """
        Purpose:
        If the Critical Volume
        Algorithm is toggled "off" in the program settings, there will only be
        one pass through this loop.  Otherwise, there will be two or more passes
        through this loop until the deco schedule is "converged" - that is when a
        comparison between the phase volume time of the present iteration and the
        last iteration is less than or equal to one minute.  This implies that
        the volume of released gas in the most recent iteration differs from the
        "critical" volume limit by an acceptably small amount.  The critical
        volume limit is set by the Critical Volume Parameter Lambda in the program
        settings (default setting is 7500 fsw-min with adjustability range from
        from 6500 to 8300 fsw-min according to Bruce Wienke).

        Side Effects:

        `self.Critical_Volume_Comparison`,
        `self.Deco_Phase_Volume_Time`,
        `self.Deco_Phase_Volume_Time`,
        `self.Deco_Stop_Depth`,
        `self.Ending_Depth`,
        `self.First_Stop_Depth`,
        `self.Helium_Pressure`,
        `self.Helium_Pressure`,
        `self.Last_Phase_Volume_Time`,
        `self.Mix_Number`,
        `self.Nitrogen_Pressure`
        `self.Nitrogen_Pressure`,
        `self.Phase_Volume_Time`,
        `self.Rate`,
        `self.Run_Time`,
        `self.Run_Time`,
        `self.Schedule_Converged`,
        `self.Segment_Number`,
        `self.Starting_Depth`,
        `self.Starting_Depth`,
        `self.Starting_Depth`,
        `self.Step_Size`,

        or

        Raises a DecompressionStepException

        Returns: None
        """
        while(True):
            # CALCULATE INITIAL ASCENT CEILING BASED ON ALLOWABLE SUPERSATURATION
            # GRADIENTS AND SET FIRST DECO STOP.  CHECK TO MAKE SURE THAT SELECTED STEP
            # SIZE WILL NOT ROUND UP FIRST STOP TO A DEPTH THAT IS BELOW THE DECO ZONE.

            self.CALC_ASCENT_CEILING()
            if self.Ascent_Ceiling_Depth <= 0.0:
                self.Deco_Stop_Depth = 0.0
            else:
                Rounding_Operation2 = (self.Ascent_Ceiling_Depth / self.Step_Size) + 0.5
                self.Deco_Stop_Depth = round(Rounding_Operation2) * self.Step_Size

            if self.Deco_Stop_Depth > self.Depth_Start_of_Deco_Zone:
                raise DecompressionStepException("ERROR! STEP SIZE IS TOO LARGE TO DECOMPRESS")

            # PERFORM A SEPARATE "PROJECTED ASCENT" OUTSIDE OF THE MAIN PROGRAM TO MAKE
            # SURE THAT AN INCREASE IN GAS LOADINGS DURING ASCENT TO THE FIRST STOP WILL
            # NOT CAUSE A VIOLATION OF THE DECO CEILING.  IF SO, ADJUST THE FIRST STOP
            # DEEPER BASED ON STEP SIZE UNTIL A SAFE ASCENT CAN BE MADE.
            # Note: this situation is a possibility when ascending from extremely deep
            # dives or due to an unusual gas mix selection.
            # CHECK AGAIN TO MAKE SURE THAT ADJUSTED FIRST STOP WILL NOT BE BELOW THE
            # DECO ZONE.

            self.PROJECTED_ASCENT(self.Depth_Start_of_Deco_Zone, self.Rate, self.Step_Size)

            if self.Deco_Stop_Depth > self.Depth_Start_of_Deco_Zone:
                raise DecompressionStepException("ERROR! STEP SIZE IS TOO LARGE TO DECOMPRESS")

            #     HANDLE THE SPECIAL CASE WHEN NO DECO STOPS ARE REQUIRED - ASCENT CAN BE
            #     MADE DIRECTLY TO THE SURFACE
            #     Write ascent data to output file and exit the Critical Volume Loop.

            if self.Deco_Stop_Depth == 0.0:
                for i in range(16):
                    self.Helium_Pressure[i] = self.He_Pressure_Start_of_Ascent[i]
                    self.Nitrogen_Pressure[i] = self.N2_Pressure_Start_of_Ascent[i]

                self.Run_Time = self.Run_Time_Start_of_Ascent
                self.Segment_Number = self.Segment_Number_Start_of_Ascent
                self.Starting_Depth = self.Depth_Change[0]
                self.Ending_Depth = 0.0
                self.GAS_LOADINGS_ASCENT_DESCENT(self.Starting_Depth, self.Ending_Depth, self.Rate)

                self.output_object.add_decompression_profile_ascent(self.Segment_Number, self.Segment_Time, self.Run_Time, self.Mix_Number, self.Deco_Stop_Depth, self.Rate)
                break

            # ASSIGN VARIABLES FOR ASCENT FROM START OF DECO ZONE TO FIRST STOP.  SAVE
            # FIRST STOP DEPTH FOR LATER USE WHEN COMPUTING THE FINAL ASCENT PROFILE

            self.Starting_Depth = self.Depth_Start_of_Deco_Zone
            self.First_Stop_Depth = self.Deco_Stop_Depth
            self.deco_stop_loop_block_within_critical_volume_loop()

            # COMPUTE TOTAL PHASE VOLUME TIME AND MAKE CRITICAL VOLUME COMPARISON
            # The deco phase volume time is computed from the run time.  The surface
            # phase volume time is computed in a subroutine based on the surfacing gas
            # loadings from previous deco loop block.  Next the total phase volume time
            # (in-water + surface) for each compartment is compared against the previous
            # total phase volume time.  The schedule is converged when the difference is
            # less than or equal to 1 minute in any one of the 16 compartments.

            # Note:  the "phase volume time" is somewhat of a mathematical concept.
            # It is the time divided out of a total integration of supersaturation
            # gradient x time (in-water and surface).  This integration is multiplied
            # by the excess bubble number to represent the amount of free-gas released
            # as a result of allowing a certain number of excess bubbles to form.
            self.Deco_Phase_Volume_Time = self.Run_Time - self.Run_Time_Start_of_Deco_Zone

            self.CALC_SURFACE_PHASE_VOLUME_TIME()

            for i in range(16):
                self.Phase_Volume_Time[i] = self.Deco_Phase_Volume_Time + self.Surface_Phase_Volume_Time[i]
                self.Critical_Volume_Comparison = abs(self.Phase_Volume_Time[i] - self.Last_Phase_Volume_Time[i])
                if self.Critical_Volume_Comparison <= 1.0:
                    self.Schedule_Converged = True

            # There are two options here.  If the Critical Volume Agorithm setting is
            # "on" and the schedule is converged, or the Critical Volume Algorithm
            # setting was "off" in the first place, the program will re-assign variables
            # to their values at the start of ascent (end of bottom time) and process
            # a complete decompression schedule once again using all the same ascent
            # parameters and first stop depth.  This decompression schedule will match
            # the last iteration of the Critical Volume Loop and the program will write
            # the final deco schedule to the output file.

            # Note: if the Critical Volume Agorithm setting was "off", the final deco
            # schedule will be based on "Initial Allowable Supersaturation Gradients."
            # If it was "on", the final schedule will be based on "Adjusted Allowable
            # Supersaturation Gradients" (gradients that are "relaxed" as a result of
            # the Critical Volume Algorithm).

            # If the Critical Volume Agorithm setting is "on" and the schedule is not
            # converged, the program will re-assign variables to their values at the
            # start of the deco zone and process another trial decompression schedule.

            if self.Schedule_Converged or self.Critical_Volume_Algorithm_Off:
                self.critical_volume_decision_tree()

            else:
                self.CRITICAL_VOLUME(self.Deco_Phase_Volume_Time)
                self.Deco_Phase_Volume_Time = 0.0
                self.Run_Time = self.Run_Time_Start_of_Deco_Zone
                self.Starting_Depth = self.Depth_Start_of_Deco_Zone
                self.Mix_Number = self.Mix_Change[0]
                self.Rate = self.Rate_Change[0]
                self.Step_Size = self.Step_Size_Change[0]
                for i in range(16):
                    self.Last_Phase_Volume_Time[i] = self.Phase_Volume_Time[i]
                    self.Helium_Pressure[i] = self.He_Pressure_Start_of_Deco_Zone[i]
                    self.Nitrogen_Pressure[i] = self.N2_Pressure_Start_of_Deco_Zone[i]
                continue
            break

    def decompression_loop(self, dive):
        """
        Purpose:
        BEGIN PROCESS OF ASCENT AND DECOMPRESSION

        Side Effects: Sets

        `self.Deco_Phase_Volume_Time`,
        `self.Deepest_Possible_Stop_Depth`,
        `self.Depth_Change`,
        `self.Depth_Change`,
        `self.He_Pressure_Start_of_Ascent`,
        `self.He_Pressure_Start_of_Deco_Zone`,
        `self.Last_Phase_Volume_Time`,
        `self.Last_Run_Time`,
        `self.Max_Actual_Gradient`
        `self.Mix_Change`,
        `self.Mix_Change`,
        `self.Mix_Number`,
        `self.N2_Pressure_Start_of_Ascent`,
        `self.N2_Pressure_Start_of_Deco_Zone`,
        `self.Number_of_Changes`,
        `self.Rate_Change`,
        `self.Rate_Change`,
        `self.Rate`,
        `self.Run_Time_Start_of_Ascent`,
        `self.Run_Time_Start_of_Deco_Zone`,
        `self.Schedule_Converged`,
        `self.Segment_Number_Start_of_Ascent`,
        `self.Starting_Depth`,
        `self.Step_Size_Change`,
        `self.Step_Size_Change`,
        `self.Step_Size`,

        Returns: None
        """
        # First, calculate the regeneration of critical radii that takes place over
        # the dive time.  The regeneration time constant has a time scale of weeks
        # so this will have very little impact on dives of normal length, but will
        # have major impact for saturation dives.
        self.NUCLEAR_REGENERATION(self.Run_Time)

        #   CALCULATE INITIAL ALLOWABLE GRADIENTS FOR ASCENT
        #   This is based on the maximum effective crushing pressure on critical radii
        #   in each compartment achieved during the dive profile.
        self.CALC_INITIAL_ALLOWABLE_GRADIENT()

        #     SAVE VARIABLES AT START OF ASCENT (END OF BOTTOM TIME) SINCE THESE WILL
        #     BE USED LATER TO COMPUTE THE FINAL ASCENT PROFILE THAT IS WRITTEN TO THE
        #     OUTPUT FILE.
        #     The VPM uses an iterative process to compute decompression schedules so
        #     there will be more than one pass through the decompression loop.

        for i in range(16):
            self.He_Pressure_Start_of_Ascent[i] = self.Helium_Pressure[i]
            self.N2_Pressure_Start_of_Ascent[i] = self.Nitrogen_Pressure[i]

        self.Run_Time_Start_of_Ascent = self.Run_Time
        self.Segment_Number_Start_of_Ascent = self.Segment_Number

        #     INPUT PARAMETERS TO BE USED FOR STAGED DECOMPRESSION AND SAVE IN ARRAYS.
        #     ASSIGN INITAL PARAMETERS TO BE USED AT START OF ASCENT
        #     The user has the ability to change mix, ascent rate, and step size in any
        #     combination at any depth during the ascent.

        for profile in dive["profile_codes"]:
            Profile_Code = profile["profile_code"]

            if Profile_Code == 99:
                self.Number_of_Changes = profile["number_of_ascent_parameter_changes"]
                self.Depth_Change = [0.0 for i in range(self.Number_of_Changes)]
                self.Mix_Change = [0.0 for i in range(self.Number_of_Changes)]
                self.Rate_Change = [0.0 for i in range(self.Number_of_Changes)]
                self.Step_Size_Change = [0.0 for i in range(self.Number_of_Changes)]

                for i, ascents in enumerate(profile["ascent_summary"]):
                    self.Depth_Change[i] = ascents["starting_depth"]
                    self.Mix_Change[i] = ascents["gasmix"]
                    self.Rate_Change[i] = ascents["rate"]
                    self.Step_Size_Change[i] = ascents["step_size"]

                self.Starting_Depth = self.Depth_Change[0]
                self.Mix_Number = self.Mix_Change[0]
                self.Rate = self.Rate_Change[0]
                self.Step_Size = self.Step_Size_Change[0]

        # CALCULATE THE DEPTH WHERE THE DECOMPRESSION ZONE BEGINS FOR THIS PROFILE
        # BASED ON THE INITIAL ASCENT PARAMETERS AND WRITE THE DEEPEST POSSIBLE
        # DECOMPRESSION STOP DEPTH TO THE OUTPUT FILE
        # Knowing where the decompression zone starts is very important.  Below
        # that depth there is no possibility for bubble formation because there
        # will be no supersaturation gradients.  Deco stops should never start
        # below the deco zone.  The deepest possible stop deco stop depth is
        # defined as the next "standard" stop depth above the point where the
        # leading compartment enters the deco zone.  Thus, the program will not
        # base this calculation on step sizes larger than 10 fsw or 3 msw.  The
        # deepest possible stop depth is not used in the program, per se, rather
        # it is information to tell the diver where to start putting on the brakes
        # during ascent.  This should be prominently displayed by any deco program.

        self.CALC_START_OF_DECO_ZONE(self.Starting_Depth, self.Rate)

        if self.units_fsw:
            if self.Step_Size < 10.0:
                rounding_op = (self.Depth_Start_of_Deco_Zone / self.Step_Size) - 0.5
                self.Deepest_Possible_Stop_Depth = round(rounding_op) * self.Step_Size
            else:
                rounding_op = (self.Depth_Start_of_Deco_Zone / 10.0) - 0.5
                self.Deepest_Possible_Stop_Depth = round(rounding_op) * 10.0

        else:
            if self.Step_Size < 3.0:
                rounding_op = (self.Depth_Start_of_Deco_Zone / self.Step_Size) - 0.5
                self.Deepest_Possible_Stop_Depth = round(rounding_op) * self.Step_Size
            else:
                rounding_op = (self.Depth_Start_of_Deco_Zone / 3.0) - 0.5
                self.Deepest_Possible_Stop_Depth = round(rounding_op) * 3.0

        #     TEMPORARILY ASCEND PROFILE TO THE START OF THE DECOMPRESSION ZONE, SAVE
        #     VARIABLES AT THIS POINT, AND INITIALIZE VARIABLES FOR CRITICAL VOLUME LOOP
        #     The iterative process of the VPM Critical Volume Algorithm will operate
        #     only in the decompression zone since it deals with excess gas volume
        #     released as a result of supersaturation gradients (not possible below the
        #     decompression zone).

        self.GAS_LOADINGS_ASCENT_DESCENT(self.Starting_Depth, self.Depth_Start_of_Deco_Zone, self.Rate)
        self.Run_Time_Start_of_Deco_Zone = self.Run_Time
        self.Deco_Phase_Volume_Time = 0.0
        self.Last_Run_Time = 0.0
        self.Schedule_Converged = False

        for i in range(16):
            self.Last_Phase_Volume_Time[i] = 0.0
            self.He_Pressure_Start_of_Deco_Zone[i] = self.Helium_Pressure[i]
            self.N2_Pressure_Start_of_Deco_Zone[i] = self.Nitrogen_Pressure[i]
            self.Max_Actual_Gradient[i] = 0.0

        self.critical_volume_loop()

    def main(self):
        """
        Purpose:
        Main decompression loop. Checks that validates the input file, initializes the data
        and loops over the dives.

        Side Effects: Sets

        `self.Max_Actual_Gradient`,
        `self.Max_Crushing_Pressure_He`,
        `self.Max_Crushing_Pressure_N2`,
        `self.Run_Time`,
        `self.Segment_Number`

        or

        Raises an InputFileException

        Returns: None
        """
        # pycallgraph.start_trace()

        self.validate_data()
        self.initialize_data()

        # START OF REPETITIVE DIVE LOOP
        # This is the largest loop in the main program and operates between Lines
        # 30 and 330.  If there is one or more repetitive dives, the program will
        # return to this point to process each repetitive dive.
        for dive in self.input_values:
            self.output_object.new_dive(dive["desc"])

            self.set_gas_mixes(dive)
            self.profile_code_loop(dive)
            self.decompression_loop(dive)

            # PROCESSING OF DIVE COMPLETE.  READ INPUT FILE TO DETERMINE IF THERE IS A
            # REPETITIVE DIVE.  IF NONE, THEN EXIT REPETITIVE LOOP.
            Repetitive_Dive_Flag = dive["repetitive_code"]

            if Repetitive_Dive_Flag == 0:
                continue

            # IF THERE IS A REPETITIVE DIVE, COMPUTE GAS LOADINGS (OFF-GASSING) DURING
            # SURFACE INTERVAL TIME.  ADJUST CRITICAL RADII USING VPM REPETITIVE
            # ALGORITHM.  RE-INITIALIZE SELECTED VARIABLES AND RETURN TO START OF
            # REPETITIVE LOOP AT LINE 30.

            elif Repetitive_Dive_Flag == 1:
                self.Surface_Interval_Time = dive["surface_interval_time_minutes"]

                self.GAS_LOADINGS_SURFACE_INTERVAL(self.Surface_Interval_Time)
                self.VPM_REPETITIVE_ALGORITHM(self.Surface_Interval_Time)

                for i in range(16):
                    self.Max_Crushing_Pressure_He[i] = 0.0
                    self.Max_Crushing_Pressure_N2[i] = 0.0
                    self.Max_Actual_Gradient[i] = 0.0

                self.Run_Time = 0.0
                self.Segment_Number = 0

                # may not be needed anymore
                continue

            else:
                raise InputFileException("Invalid repetitive dive flag %d. Must be 0 (don't repeat) or 1 (repeat)" % (Repetitive_Dive_Flag))

        # pycallgraph.make_dot_graph('pycallgraph.png')


class Output(object):
    """Provides a convenient dumping ground for output. Can be exported to json and html
    """

    def __init__(self, state_object):
        self.output = []
        self.current_dive = None
        self.state_object = state_object

    def new_dive(self, description):
        """Creates a new dive hash table and updates the current dive to
        point at the newly created hash"""
        self.output.append({"desc": description,
                            "time": str(datetime.datetime.now().strftime("%B %d, %Y at %H:%M")),  # make this cofigurable
                            "gasmix": [],
                            "dive_profile": [],
                            "decompression_profile": []})

        if self.current_dive is None:
            self.current_dive = 0
        else:
            self.current_dive += 1

    def add_gasmix(self, oxygen, nitrogen, helium):
        """Adds a new gasmix to the current dive"""
        self.output[self.current_dive]["gasmix"].append({"oxygen": oxygen,
                                                         "nitrogen": nitrogen,
                                                         "helium": helium})

    def add_dive_profile_entry_ascent(self, Segment_Number, Segment_Time, Run_Time, Mix_Number, Depth):
        """Adds a new ascent entry to the dive profile table"""
        split = "%d|%6.1f |%6.1f | %d |  | | | | %6.1f| " \
                % (Segment_Number, Segment_Time, Run_Time, Mix_Number, Depth)
        self.output[self.current_dive]["dive_profile"].append(split.split("|"))

    def add_dive_profile_entry_descent(self, Segment_Number, Segment_Time, Run_Time, Mix_Number, Word, Starting_Depth, Ending_Depth, Rate):
        """Adds a new descent entry to the dive profile table"""
        split = "%d|%6.1f |%6.1f | %d | %s | %6.1f| %6.1f| %6.1f| |" \
                % (Segment_Number, Segment_Time, Run_Time, Mix_Number, Word, Starting_Depth, Ending_Depth, Rate)
        self.output[self.current_dive]["dive_profile"].append(split.split("|"))

    def add_decompression_profile_ascent(self, Segment_Number, Segment_Time, Run_Time, Mix_Number, Deco_Stop_Depth, Rate):
        """Adds a new ascent entry to the decompression table"""
        split = "%d|%5.1f|%6.1f|%d|%d|%d | | |" \
                % (Segment_Number, Segment_Time, Run_Time, Mix_Number, Deco_Stop_Depth, Rate)
        self.output[self.current_dive]["decompression_profile"].append(split.split("|"))

    def add_decompression_profile_constant(self, Segment_Number, Segment_Time, Run_Time, Mix_Number, Deco_Stop_Depth, Stop_Time):
        """Adds a new constant depth entry to the decompression table"""

        split = "%d|%5.1f|%6.1f|%d| | |%d|%d|%d|" \
                % (Segment_Number, Segment_Time, Run_Time, Mix_Number, Deco_Stop_Depth, Stop_Time, int(Run_Time))
        self.output[self.current_dive]["decompression_profile"].append(split.split("|"))

    def to_html(self, filename=None):
        """Export the dive output object to an HTML file or the console"""

        output = ""
        output += self.html_header()

        for dive in self.output:
            output += self.html_description_time_and_gasmix(dive)
            output += self.html_dive_table(dive)
            output += self.html_decompression_table(dive)

        output += self.html_footer()
        if filename:
            f = open(filename, "w")
            f.write(output)
            f.close()
        else:
            print output

    def html_description_time_and_gasmix(self, dive):
        return_string = ""
        return_string += "<p>"

        return_string += dive["desc"]
        return_string += "<br/>"
        return_string += dive["time"]
        return_string += "<br/>"
        return_string += "Gasmix Summary"
        return_string += "</p>"
        return_string += """
            <table border="2" >
            <tr>
            <th>Gasmix number</th>
            <th>O2</th>
            <th>He</th>
            <th>N2</th>
            </tr>
            """

        for i, gas in enumerate(dive["gasmix"]):
            return_string += """<tr>
            <td>%d</td>
            <td>%5.3f</td>
            <td>%5.3f</td>
            <td>%5.3f</td>
            </tr>
            """ % (i + 1, gas["oxygen"], gas["helium"], gas["nitrogen"])
        return_string += "</table>"

        return return_string

    def html_dive_table(self, dive):
        return_string = ""
        return_string += "<p>Dive Profile</p>"
        return_string += """<table border="2" >
        <tr>
        <th>Segment</th>
        <th>Segment Time (min)</th>
        <th>Run Time (min)</th>
        <th>Gasmix Used #</th>
        <th>Ascent or Descent</th>
        <th>From Depth (%s)</th>
        <th>To Depth (%s)</th>
        <th>Rate +Dn/-Up (%s)</th>
        <th>Constant Depth (%s)</th>
        </tr>
        """ % (self.state_object.Units_Word1, self.state_object.Units_Word1, self.state_object.Units_Word2, self.state_object.Units_Word1)

        for d in dive["dive_profile"]:
            return_string += "<tr>"
            for elem in d:
                return_string += "<td>%s</td>" % (str(elem))
            return_string += "</tr>"

        return_string += "</table>"
        return return_string

    def html_decompression_table(self, dive):
        return_string = ""
        return_string += "<p>"

        return_string += "DECOMPRESSION PROFILE"
        return_string += "<br/>"
        return_string += "Leading compartment enters the decompression zone at, %.1f %s" % (self.state_object.Depth_Start_of_Deco_Zone, self.state_object.Units_Word1)
        return_string += "<br/>"

        return_string += "Deepest possible decompression stop is %.1f %s" % (self.state_object.Deepest_Possible_Stop_Depth, self.state_object.Units_Word1)
        return_string += "<br/>"
        return_string += "</p>"

        return_string += """<table border="2">
        <tr>
        <th>Segment #</th>
        <th>Segment Time (min)</th>
        <th>Run Time(min)</th>
        <th>Gasmix Used #</th>
        <th>Ascent to (%s)</th>
        <th>Ascent Rate(%s)</th>
        <th>Deco Stop(%s)</th>
        <th>Stop Time (min)</th>
        <th>Run Time(min)</th>
        </tr>
        """ % (self.state_object.Units_Word1, self.state_object.Units_Word2, self.state_object.Units_Word1)

        for d in dive["decompression_profile"]:
            return_string += "<tr>"
            for elem in d:
                return_string += "<td>%s</td>" % (str(elem))
            return_string += "</tr>"

        return_string += "</table>"
        return return_string

    def html_header(self):
        header = '''<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN"
        "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
        <html xmlns="http://www.w3.org/1999/xhtml"
        lang="en" xml:lang="en">
        <head>
        <meta http-equiv="content-type" content="text/html;charset=utf-8" />
        <title>VPMB Dive chart</title>
        <style type="text/css">

        html { font-family: Times, serif; font-size: 12pt; }
        .title  { text-align: center; }
        .todo   { color: red; }
        .done   { color: green; }
        .tag    { background-color:lightblue; font-weight:normal }
        .target { }
        .timestamp { color: grey }
        .timestamp-kwd { color: CadetBlue }
        p.verse { margin-left: 3% }
        pre {
          border: 1pt solid #AEBDCC;
          background-color: #F3F5F7;
          padding: 5pt;
          font-family: courier, monospace;
          font-size: 90%;
          overflow:auto;
        }
        table { border-collapse:collapse;  }
        td, th { vertical-align: top; }
        dt { font-weight: bold; }
        div.figure { padding: 0.5em; }
        div.figure p { text-align: center; }
        .linenr { font-size:smaller }
        .code-highlighted {background-color:#ffff00;}

        </style>
        </head><body>
        '''
        return header

    def html_footer(self):
        return "</body></html>"

    def to_json(self, file_pointer):
        """Export the json formatted output to a file
        """
        json.dump(self.output, file_pointer)

    def get_json(self):
        """Return the output JSON"""
        return self.output


# functions
def SCHREINER_EQUATION(Initial_Inspired_Gas_Pressure, Rate_Change_Insp_Gas_Pressure, Interval_Time, Gas_Time_Constant, Initial_Gas_Pressure):
    """Function for ascent and descent gas loading calculations"""
    return Initial_Inspired_Gas_Pressure + Rate_Change_Insp_Gas_Pressure * (Interval_Time - 1.0 / Gas_Time_Constant) - (Initial_Inspired_Gas_Pressure - Initial_Gas_Pressure - Rate_Change_Insp_Gas_Pressure / Gas_Time_Constant) * exp(-Gas_Time_Constant * Interval_Time)


def HALDANE_EQUATION(Initial_Gas_Pressure, Inspired_Gas_Pressure, Gas_Time_Constant, Interval_Time):
    """Function for gas loading calculations at a constant depth"""
    return Initial_Gas_Pressure + (Inspired_Gas_Pressure - Initial_Gas_Pressure) * (1.0 - exp(-Gas_Time_Constant * Interval_Time))


def RADIUS_ROOT_FINDER(A, B, C, Low_Bound, High_Bound):
    """
    Purpose: This subroutine is a "fail-safe" routine that combines the
    Bisection Method and the Newton-Raphson Method to find the desired root.
    This hybrid algorithm takes a bisection step whenever Newton-Raphson would
    take the solution out of bounds, or whenever Newton-Raphson is not
    converging fast enough.  Source:  "Numerical Recipes in Fortran 77",
    Cambridge University Press, 1992.

    Side Effects: None

    or

    Raises a RootException, MaxIterationException

    Returns: A floating point value
    """
    # BEGIN CALCULATIONS BY MAKING SURE THAT THE ROOT LIES WITHIN BOUNDS
    # In this case we are solving for radius in a cubic equation of the form,
    # Ar^3 - Br^2 - C = 0.  The coefficients A, B, and C were passed to this
    # subroutine as arguments.
    Function_at_Low_Bound = Low_Bound * (Low_Bound * (A * Low_Bound - B)) - C

    Function_at_High_Bound = High_Bound * (High_Bound * (A * High_Bound - B)) - C

    if Function_at_Low_Bound > 0.0 and Function_at_High_Bound > 0.0:
        raise RootException("ERROR! ROOT IS NOT WITHIN BRACKETS")

    # Next the algorithm checks for special conditions and then prepares for
    # the first bisection.
    if Function_at_Low_Bound < 0.0 and Function_at_High_Bound < 0.0:
        raise RootException("ERROR! ROOT IS NOT WITHIN BRACKETS")

    if Function_at_Low_Bound == 0.0:
        return Low_Bound
    elif Function_at_High_Bound == 0.0:
        return High_Bound
    elif Function_at_Low_Bound < 0.0:
        Radius_at_Low_Bound = Low_Bound
        Radius_at_High_Bound = High_Bound
    else:
        Radius_at_High_Bound = Low_Bound
        Radius_at_Low_Bound = High_Bound

    Ending_Radius = 0.5 * (Low_Bound + High_Bound)
    Last_Diff_Change = abs(High_Bound - Low_Bound)
    Differential_Change = Last_Diff_Change

    # At this point, the Newton-Raphson Method is applied which uses a function
    # and its first derivative to rapidly converge upon a solution.
    # Note: the program allows for up to 100 iterations.  Normally an exit will
    # be made from the loop well before that number.  If, for some reason, the
    # program exceeds 100 iterations, there will be a pause to alert the user.
    # When a solution with the desired accuracy is found, exit is made from the
    # loop by returning to the calling program.  The last value of ending
    # radius has been assigned as the solution.

    Function = Ending_Radius * (Ending_Radius * (A * Ending_Radius - B)) - C

    Derivative_of_Function = Ending_Radius * (Ending_Radius * 3.0 * A - 2.0 * B)

    for i in range(100):
        if((((Ending_Radius - Radius_at_High_Bound) * Derivative_of_Function - Function) *
            ((Ending_Radius - Radius_at_Low_Bound) * Derivative_of_Function - Function) >= 0.0)
           or (abs(2.0 * Function) > (abs(Last_Diff_Change * Derivative_of_Function)))):

            Last_Diff_Change = Differential_Change
            Differential_Change = 0.5 * (Radius_at_High_Bound - Radius_at_Low_Bound)

            Ending_Radius = Radius_at_Low_Bound + Differential_Change
            if Radius_at_Low_Bound == Ending_Radius:
                return Ending_Radius
        else:
            Last_Diff_Change = Differential_Change
            Differential_Change = Function / Derivative_of_Function
            Last_Ending_Radius = Ending_Radius
            Ending_Radius = Ending_Radius - Differential_Change
            if Last_Ending_Radius == Ending_Radius:
                return Ending_Radius
        if abs(Differential_Change) < 1.0E-12:
            return Ending_Radius
        Function = Ending_Radius * (Ending_Radius * (A * Ending_Radius - B)) - C

        Derivative_of_Function = Ending_Radius * (Ending_Radius * 3.0 * A - 2.0 * B)

        if Function < 0.0:
            Radius_at_Low_Bound = Ending_Radius
        else:
            Radius_at_High_Bound = Ending_Radius

    raise MaxIterationException('ERROR! ROOT SEARCH EXCEEDED MAXIMUM ITERATIONS')


def CALC_BAROMETRIC_PRESSURE(Altitude, units_fsw):
    """
    Purpose: This function calculates barometric pressure at altitude based on the
    publication "U.S. Standard Atmosphere, 1976", U.S. Government Printing
    Office, Washington, D.C. The basis for this code is a Fortran 90 program
    written by Ralph L. Carmichael (retired NASA researcher) and endorsed by
    the National Geophysical Data Center of the National Oceanic and
    Atmospheric Administration.  It is available for download free from
    Public Domain Aeronautical Software at:  http://www.pdas.com/atmos.htm

    Side Effects: None

    Returns: A floating point value
    """

    Radius_of_Earth = 6369.0  # kilometers
    Acceleration_of_Gravity = 9.80665  # meters/second^2
    Molecular_weight_of_Air = 28.9644  # mols
    Gas_Constant_R = 8.31432  # Joules/mol*deg Kelvin
    Temp_at_Sea_Level = 288.15  # degrees Kelvin

    Pressure_at_Sea_Level_Fsw = 33.0  # feet of seawater based on 101325 Pa at sea level (Standard Atmosphere)

    Pressure_at_Sea_Level_Msw = 10.0  # meters of seawater based on 100000 Pa at sea level (European System)

    # Change in Temp deg Kelvin with
    # change in geopotential altitude,
    # valid for first layer of atmosphere
    # up to 11 kilometers or 36,000 feet
    Temp_Gradient = -6.5

    GMR_Factor = Acceleration_of_Gravity * Molecular_weight_of_Air / Gas_Constant_R

    if units_fsw:
        Altitude_Feet = Altitude
        Altitude_Kilometers = Altitude_Feet / 3280.839895
        Pressure_at_Sea_Level = Pressure_at_Sea_Level_Fsw
    else:
        Altitude_Meters = Altitude
        Altitude_Kilometers = Altitude_Meters / 1000.0
        Pressure_at_Sea_Level = Pressure_at_Sea_Level_Msw

    Geopotential_Altitude = (Altitude_Kilometers * Radius_of_Earth) / (Altitude_Kilometers + Radius_of_Earth)
    Temp_at_Geopotential_Altitude = Temp_at_Sea_Level + Temp_Gradient * Geopotential_Altitude

    Barometric_Pressure = Pressure_at_Sea_Level * exp(log(Temp_at_Sea_Level / Temp_at_Geopotential_Altitude) * GMR_Factor / Temp_Gradient)
    return Barometric_Pressure


def parse_settings():
    """
    Purpose:
    Use OptParse to parse the command line switches and return the
    'options' and 'args' objects

    Side Effects: None

    Returns: `options` and `args` objects containg the command line arguments
    """

    # Load settings and inputs
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-i", action="store", dest="input_file_name",
                      default="vpm_decompression_input.json",
                      help="Input file containing dive information")

    parser.add_option("-o", action="store", dest="output_file_name",
                      default="output.html",
                      help="Output file for dive log")

    parser.add_option("-j", action="store", dest="json_output", default=None,
                      help="Output json instead of html")

    (options, args) = parser.parse_args()

    return options, args

# if they ran this at the command line, output the parse the command line
# options and output the results
if __name__ == '__main__':
    options, args = parse_settings()
    program_state = DiveState(input_file_name=options.input_file_name)
    program_state.main()

    if options.json_output:
        json_out = open(options.json_output, "w")
        program_state.output_object.to_json(json_out)
        json_out.close()
    else:
        program_state.output_object.to_html(options.output_file_name)
