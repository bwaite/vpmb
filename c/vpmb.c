#include "vpmb.h"
/* Copyright 2010, Bryan Waite, Erik C. Baker. All rights reserved.
 * Redistribution and use in source and binary forms, with or without modification, are
 * permitted provided that the following conditions are met:
 *
 *    1. Redistributions of source code must retain the above copyright notice, this list of
 *       conditions and the following disclaimer.
 *
 *    2. Redistributions in binary form must reproduce the above copyright notice, this list
 *       of conditions and the following disclaimer in the documentation and/or other materials
 *       provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY Bryan Waite, Erik C. Baker ``AS IS'' AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
 * FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL Bryan Waite, or Erik C. Baker OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
 * ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 * The views and conclusions contained in the software and documentation are those of the
 * authors and should not be interpreted as representing official policies, either expressed
 * or implied, of Bryan Waite, or Erik C. Baker.
 *
 */

const short TRUE = 1;
const short FALSE = 0;

/* errors */
const int ALLGOOD = 256; /* My Lucky Number */

const int BADJSON = -1;
const int GOODJSON = 1;

const int INVALIDDATA = -2;
const int VALIDDATA = 2;

const int ROOTERROR = -3;
const int ROOTFOUND = 3;
const int BADALTITUDE = -4;

const int OFFGASSINGERROR = -5;
const int BADDECOSTOP = -6;
const int NOFILE = -7;

/* Constants */
const double ATM = 101325.0;
const double fraction_inert_gas = 0.79;
const double Helium_Half_Time[] = {1.88, 3.02, 4.72, 6.99, 10.21, 14.48, 20.53, 29.11, 41.20, 55.19, 70.69, 90.34, 115.29, 147.42, 188.24, 240.03};
const double Nitrogen_Half_Time[] = {5.0, 8.0, 12.5, 18.5, 27.0, 38.3, 54.3, 77.0, 109.0, 146.0, 187.0, 239.0, 305.0, 390.0, 498.0, 635.0};
const int Buhlmann_Compartments = 16;

const int DONE_PROFILE = 1;
const int NOT_DONE_PROFILE = -1;
const int DEPTH_CHANGED = 0;

const double meters_per_minute_change = 2.5; /* when deciding whether to ascend or descend (schreiner_equation) or staying at a constant depth (haldane_equation), see if the depth change over the last minute is +- this value */



void vpmb_failure()
{
    /* turn this into a better way to handle failures for a real program */
    exit(EXIT_FAILURE);
}

double vpmb_get_setpoint_or_null(const cJSON *setpoint)
{
    if (setpoint == NULL) {
        return 0.0;
    }
    return setpoint->valuedouble;
}

int vpmb_input_add_dive_profiles(dive_profile *in, cJSON *profile)
{
    int i;
    cJSON *current_ascent, *ascent_fields;

    in->profile_code = cJSON_GetObjectItem(profile, "profile_code")->valueint;
    if (in->profile_code == 1) {
        in->starting_depth = cJSON_GetObjectItem(profile, "starting_depth")->valuedouble;
        in->ending_depth  = cJSON_GetObjectItem(profile, "ending_depth")->valuedouble;
        in->rate  = cJSON_GetObjectItem(profile, "rate")->valuedouble;
        in->gasmix  = cJSON_GetObjectItem(profile, "gasmix")->valueint;
        in->setpoint = vpmb_get_setpoint_or_null(cJSON_GetObjectItem(profile, "setpoint"));
    } else if (in->profile_code == 2) {
        in->depth = cJSON_GetObjectItem(profile, "depth")->valuedouble;
        in->run_time_at_end_of_segment  = cJSON_GetObjectItem(profile, "run_time_at_end_of_segment")->valuedouble;
        in->gasmix = cJSON_GetObjectItem(profile, "gasmix")->valueint;
        in->setpoint  = vpmb_get_setpoint_or_null(cJSON_GetObjectItem(profile, "setpoint"));

    } else if (in->profile_code == 99) {
        in->number_of_ascent_parameter_changes = cJSON_GetObjectItem(profile, "number_of_ascent_parameter_changes")->valueint;
        in->ascents = calloc(in->number_of_ascent_parameter_changes, sizeof(ascent_summary));

        if(in->ascents == NULL) {
            vpmb_failure();
        }

        ascent_fields = cJSON_GetObjectItem(profile, "ascent_summary");
        for(i = 0; i < in->number_of_ascent_parameter_changes; i++) {
            current_ascent = cJSON_GetArrayItem(ascent_fields, i);

            in->ascents[i].starting_depth = cJSON_GetObjectItem(current_ascent, "starting_depth")->valuedouble;
            in->ascents[i].gasmix = cJSON_GetObjectItem(current_ascent, "gasmix")->valueint;
            in->ascents[i].rate = cJSON_GetObjectItem(current_ascent, "rate")->valuedouble;
            in->ascents[i].step_size = cJSON_GetObjectItem(current_ascent, "step_size")->valuedouble;
            in->ascents[i].setpoint = vpmb_get_setpoint_or_null(cJSON_GetObjectItem(current_ascent, "setpoint"));
        }
    } else {
        return BADJSON;
    }
    return GOODJSON;
}

int vpmb_input_add_dive(single_dive *in, cJSON *dive)
{
    int i;
    cJSON *current_gasmix, *gasmix_summary_fields, *profile_code_fields, *current_profile;

    strlcpy(in->desc, cJSON_GetObjectItem(dive, "desc")->valuestring, sizeof(in->desc));
    in->num_gas_mixes = cJSON_GetObjectItem(dive, "num_gas_mixes")->valueint;
    in->repetitive_code = cJSON_GetObjectItem(dive, "repetitive_code")->valueint;

    if(in->repetitive_code == 1) {
        in->surface_interval_time_minutes = cJSON_GetObjectItem(dive, "surface_interval_time_minutes")->valuedouble;
    }

    in->gasmixes = calloc(in->num_gas_mixes, sizeof(gasmix_summary));
    gasmix_summary_fields = cJSON_GetObjectItem(dive, "gasmix_summary");
    for(i = 0; i < in->num_gas_mixes; i++) {

        current_gasmix = cJSON_GetArrayItem(gasmix_summary_fields, i);
        in->gasmixes[i].fraction_O2 = cJSON_GetObjectItem(current_gasmix, "fraction_O2")->valuedouble;
        in->gasmixes[i].fraction_He = cJSON_GetObjectItem(current_gasmix, "fraction_He")->valuedouble;
        in->gasmixes[i].fraction_N2 = cJSON_GetObjectItem(current_gasmix, "fraction_N2")->valuedouble;
    }

    profile_code_fields = cJSON_GetObjectItem(dive, "profile_codes");
    in->num_profile_codes = cJSON_GetArraySize(profile_code_fields);
    in->dive_profiles = calloc(in->num_profile_codes, sizeof(dive_profile));
    for(i = 0; i < in->num_profile_codes; i++) {

        current_profile = cJSON_GetArrayItem(profile_code_fields, i);

        vpmb_input_add_dive_profiles(&(in->dive_profiles[i]), current_profile);
    }
    return GOODJSON;
}

void vpmb_free_dives(json_input *in)
{
    int i, j;
    for(i = 0; i < in->number_of_dives; i++) {
        for (j = 0; j < in->dives[i].num_profile_codes; j++) {
            if(in->dives[i].dive_profiles[j].profile_code == 99) {
                free(in->dives[i].dive_profiles[j].ascents);
            }
        }
        free(in->dives[i].gasmixes);
        free(in->dives[i].dive_profiles);
    }
    free(in->dives);
}

int vpmb_load_from_json(json_input *in, const char* filename)
{
    int i;
    size_t len;
    char *data;
    cJSON *json, *current_item, *current_dive, *setpoint_is_bar;
    /* find out the size of the file and allocate enough memory to hold it */
    FILE *f=fopen(filename,"rb");
    if (f == NULL) {
        return NOFILE;
    }

    fseek(f, 0, SEEK_END);
    len=ftell(f);
    fseek(f,0,SEEK_SET);

    data=malloc(len+1);
    if (fread(data,1,len,f) < len) {
        free(data);
        fclose(f);
        return BADJSON;
    }
    fclose(f);

    json=cJSON_Parse(data);
    current_item = cJSON_GetObjectItem(json, "input");

    in->number_of_dives = cJSON_GetArraySize(current_item);
    in->dives = calloc(in->number_of_dives, sizeof(single_dive));

    for (i = 0; i < cJSON_GetArraySize(current_item); i++) {
        current_dive = cJSON_GetArrayItem(current_item, i);

        if (current_dive != NULL) {
            vpmb_input_add_dive(&(in->dives[i]), current_dive);
        } else {
            free(data);
            return BADJSON;
        }
    }

    /* copy altitude */
    current_item = cJSON_GetObjectItem(json, "altitude");

    in->Altitude_of_Dive = cJSON_GetObjectItem(current_item, "Altitude_of_Dive")->valuedouble;

    strlcpy(in->Diver_Acclimatized_at_Altitude, cJSON_GetObjectItem(current_item, "Diver_Acclimatized_at_Altitude")->valuestring, sizeof(in->Diver_Acclimatized_at_Altitude));

    in->Starting_Acclimatized_Altitude = cJSON_GetObjectItem(current_item, "Starting_Acclimatized_Altitude")->valuedouble;
    in->Ascent_to_Altitude_Hours = cJSON_GetObjectItem(current_item, "Ascent_to_Altitude_Hours")->valuedouble;
    in->Hours_at_Altitude_Before_Dive = cJSON_GetObjectItem(current_item, "Hours_at_Altitude_Before_Dive")->valuedouble;

    /* copy settings */
    current_item = cJSON_GetObjectItem(json, "settings");

    strlcpy(in->Units, cJSON_GetObjectItem(current_item, "Units")->valuestring, sizeof(in->Units));

    setpoint_is_bar = cJSON_GetObjectItem(current_item, "SetPoint_Is_Bar");
    if (setpoint_is_bar != NULL) {
        if (cJSON_True == setpoint_is_bar->type) {
            in->SetPoint_Is_Bar = TRUE;
        } else if (cJSON_False == setpoint_is_bar->type) {
            in->SetPoint_Is_Bar = FALSE;
        } else {
            free(data);
            return BADJSON;
        }
    }
    strlcpy(in->Altitude_Dive_Algorithm, cJSON_GetObjectItem(current_item, "Altitude_Dive_Algorithm")->valuestring, sizeof(in->Altitude_Dive_Algorithm));

    in->Minimum_Deco_Stop_Time = cJSON_GetObjectItem(current_item, "Minimum_Deco_Stop_Time")->valuedouble;
    in->Critical_Radius_N2_Microns = cJSON_GetObjectItem(current_item, "Critical_Radius_N2_Microns")->valuedouble;
    in->Critical_Radius_He_Microns = cJSON_GetObjectItem(current_item, "Critical_Radius_He_Microns")->valuedouble;
    strlcpy(in->Critical_Volume_Algorithm, cJSON_GetObjectItem(current_item, "Critical_Volume_Algorithm")->valuestring, sizeof(in->Critical_Volume_Algorithm));
    in->Crit_Volume_Parameter_Lambda = cJSON_GetObjectItem(current_item, "Crit_Volume_Parameter_Lambda")->valuedouble;
    in->Gradient_Onset_of_Imperm_Atm = cJSON_GetObjectItem(current_item, "Gradient_Onset_of_Imperm_Atm")->valuedouble;
    in->Surface_Tension_Gamma = cJSON_GetObjectItem(current_item, "Surface_Tension_Gamma")->valuedouble;
    in->Skin_Compression_GammaC = cJSON_GetObjectItem(current_item, "Skin_Compression_GammaC")->valuedouble;
    in->Regeneration_Time_Constant = cJSON_GetObjectItem(current_item, "Regeneration_Time_Constant")->valuedouble;
    in->Pressure_Other_Gases_mmHg = cJSON_GetObjectItem(current_item, "Pressure_Other_Gases_mmHg")->valuedouble;

    cJSON_Delete(json);
    free(data);

    return GOODJSON;
}

/* VPM RELATED STUFF */

/**
 * \brief Function for ascent and descent gas loading calculations.
 *
 * Based on the derivations done by Erik Baker (Note: ascents must be given in negative numbers (ie: -10 for 10 feet up), as descents use positive values).
 *
 * Calculates
 * \f$P(t) = P_{Io} + R(t - 1 / k) - (P_{Io} − P_O − R / k) e^{kt}\f$
 * \return Pressure at interval time
 */
double vpmb_schreiner_equation(double Initial_Inspired_Gas_Pressure, double Rate_Change_Insp_Gas_Pressure,
                               double Interval_Time, double Gas_Time_Constant, double Initial_Gas_Pressure)
{

    return Initial_Inspired_Gas_Pressure + Rate_Change_Insp_Gas_Pressure * (Interval_Time - 1.0 / Gas_Time_Constant) - (Initial_Inspired_Gas_Pressure - Initial_Gas_Pressure - Rate_Change_Insp_Gas_Pressure/Gas_Time_Constant) * exp (-Gas_Time_Constant * Interval_Time);
}

/**
 * \brief Function for gas loading calculations at a constant depth.
 *
 * Based on the derivations by Erik Baker.
 *
 * Calculates
 * \f$P = P_0 + (P_I - P_0)(1 - e^{kt})\f$
 * \return Pressure at interval time.
 */
double vpmb_haldane_equation(double Initial_Gas_Pressure, double Inspired_Gas_Pressure, double Gas_Time_Constant, double Interval_Time)
{

    return Initial_Gas_Pressure + (Inspired_Gas_Pressure - Initial_Gas_Pressure) * (1.0 - exp(-Gas_Time_Constant * Interval_Time));
}



int vpmb_radius_root_finder(double A, double B, double C, double Low_Bound, double High_Bound, double *result)
{
    int i;
    const double iterations = 100;
    double Function_at_Low_Bound = Low_Bound * (Low_Bound * (A * Low_Bound - B)) - C;
    double Function_at_High_Bound = High_Bound * (High_Bound * (A * High_Bound - B)) - C;
    double Radius_at_Low_Bound, Radius_at_High_Bound;
    double Ending_Radius, Last_Diff_Change, Differential_Change;
    double Function, Derivative_of_Function, Last_Ending_Radius;

    if ((Function_at_Low_Bound > 0.0) && (Function_at_High_Bound > 0.0)) {
        return ROOTERROR;
    }

    if ((Function_at_Low_Bound < 0.0) && (Function_at_High_Bound < 0.0)) {
        return ROOTERROR;
    }
    if (Function_at_Low_Bound == 0.0) {
        *result = Low_Bound;
        return ROOTFOUND;
    } else if (Function_at_High_Bound == 0.0) {
        *result = High_Bound;
        return ROOTFOUND;
    } else if (Function_at_Low_Bound < 0.0) {
        Radius_at_Low_Bound = Low_Bound;
        Radius_at_High_Bound = High_Bound;
    } else {
        Radius_at_High_Bound = Low_Bound;
        Radius_at_Low_Bound = High_Bound;
    }

    Ending_Radius = 0.5 * (Low_Bound + High_Bound);
    Last_Diff_Change = fabs(High_Bound - Low_Bound);
    Differential_Change = Last_Diff_Change;

    /*At this point, the Newton-Raphson Method is applied which uses a function
      and its first derivative to rapidly converge upon a solution.
      Note: the program allows for up to 100 iterations.  Normally an exit will
      be made from the loop well before that number.  If, for some reason, the
      program exceeds 100 iterations, there will be a pause to alert the user.
      When a solution with the desired accuracy is found, exit is made from the
      loop by returning to the calling program.  The last value of ending
      radius has been assigned as the solution.*/

    Function = Ending_Radius * (Ending_Radius * (A * Ending_Radius - B)) - C;
    Derivative_of_Function = Ending_Radius * (Ending_Radius * 3.0 * A - 2.0 * B);

    for(i = 0; i < iterations; i++) {
        if((((Ending_Radius-Radius_at_High_Bound) * Derivative_of_Function - Function) *
            ((Ending_Radius-Radius_at_Low_Bound) * Derivative_of_Function - Function) >= 0.0)
           || (fabs(2.0 * Function) > (fabs(Last_Diff_Change*Derivative_of_Function)))) {

            Last_Diff_Change = Differential_Change;
            Differential_Change = 0.5 * (Radius_at_High_Bound - Radius_at_Low_Bound);

            Ending_Radius = Radius_at_Low_Bound + Differential_Change;
            if (Radius_at_Low_Bound == Ending_Radius) {
                *result = Ending_Radius;
                return ROOTFOUND;
            }
        }

        else {
            Last_Diff_Change = Differential_Change;
            Differential_Change = Function / Derivative_of_Function;
            Last_Ending_Radius = Ending_Radius;
            Ending_Radius = Ending_Radius - Differential_Change;
            if (Last_Ending_Radius == Ending_Radius) {
                *result = Ending_Radius;
                return ROOTFOUND;
            }
            if (fabs(Differential_Change) < 1.0E-12) {
                *result = Ending_Radius;
                return ROOTFOUND;
            }
        }

        Function = Ending_Radius * (Ending_Radius * (A * Ending_Radius - B)) - C;
        Derivative_of_Function = Ending_Radius * (Ending_Radius * 3.0 * A - 2.0 * B);

        if (Function < 0.0) {
            Radius_at_Low_Bound = Ending_Radius;
        } else {
            Radius_at_High_Bound = Ending_Radius;
        }
    }
    return ROOTERROR;
}

/** \brief Calculate the barometric pressure at the current altitude
 *
 * Based on the FORTRAN code by Ralph L. Carmichael available at http://www.pdas.com/atmos.html ,
 * and the version provided by Erik Baker.
 * \returns The barometric pressure for the given altitude, in the units specified by the user.
 */
double vpmb_calc_barometric_pressure(double Altitude, BOOL units_fsw)
{
    double Radius_of_Earth = 6369.0;
    double Acceleration_of_Gravity = 9.80665;
    double Molecular_weight_of_Air = 28.9644;
    double Gas_Constant_R = 8.31432;
    double Temp_at_Sea_Level = 288.15;
    double Pressure_at_Sea_Level_Fsw = 33.0;
    double Pressure_at_Sea_Level_Msw = 10.0;
    double Temp_Gradient = -6.5;
    double GMR_Factor = Acceleration_of_Gravity * Molecular_weight_of_Air / Gas_Constant_R;
    double Pressure_at_Sea_Level, Geopotential_Altitude, Temp_at_Geopotential_Altitude, Barometric_Pressure;
    double Altitude_Kilometers;

    if (units_fsw == TRUE) {
        double Altitude_Feet = Altitude;
        Altitude_Kilometers = Altitude_Feet / 3280.839895;
        Pressure_at_Sea_Level = Pressure_at_Sea_Level_Fsw;
    } else {
        double Altitude_Meters = Altitude;
        Altitude_Kilometers = Altitude_Meters / 1000.0;
        Pressure_at_Sea_Level = Pressure_at_Sea_Level_Msw;
    }

    Geopotential_Altitude =  (Altitude_Kilometers * Radius_of_Earth) / (Altitude_Kilometers + Radius_of_Earth);
    Temp_at_Geopotential_Altitude = Temp_at_Sea_Level + Temp_Gradient * Geopotential_Altitude;

    Barometric_Pressure = Pressure_at_Sea_Level * exp(log(Temp_at_Sea_Level / Temp_at_Geopotential_Altitude) * GMR_Factor / Temp_Gradient);
    return Barometric_Pressure;
}

double vpmb_calc_deco_ceiling(dive_state *dive)
{
    int i;
    double Compartment_Deco_Ceiling[16];
    double Gas_Loading, Tolerated_Ambient_Pressure, Deco_Ceiling_Depth;

    for(i=0; i < 16; i++) {
        Gas_Loading = dive->Helium_Pressure[i] + dive->Nitrogen_Pressure[i];

        if (Gas_Loading > 0.0) {
            double Weighted_Allowable_Gradient = (dive->Deco_Gradient_He[i] * dive->Helium_Pressure[i] + dive->Deco_Gradient_N2[i] * dive->Nitrogen_Pressure[i]) / (dive->Helium_Pressure[i] + dive->Nitrogen_Pressure[i]);

            Tolerated_Ambient_Pressure = (Gas_Loading + dive->Constant_Pressure_Other_Gases) - Weighted_Allowable_Gradient;
        } else {
            double Weighted_Allowable_Gradient = min(dive->Deco_Gradient_He[i], dive->Deco_Gradient_N2[i]);
            Tolerated_Ambient_Pressure = dive->Constant_Pressure_Other_Gases - Weighted_Allowable_Gradient;
        }
        if(Tolerated_Ambient_Pressure < 0.0) {
            Tolerated_Ambient_Pressure = 0.0;
        }


        Compartment_Deco_Ceiling[i] = Tolerated_Ambient_Pressure - dive->Barometric_Pressure;
    }

    Deco_Ceiling_Depth = Compartment_Deco_Ceiling[0];

    for(i=1; i < 16; i++) {
        Deco_Ceiling_Depth = max(Deco_Ceiling_Depth, Compartment_Deco_Ceiling[i]);
    }

    return Deco_Ceiling_Depth;
}



int vpmb_vpm_altitude_dive_algorithm(json_input *input, dive_state *dive)
{
    int i;
    double Ascent_to_Altitude_Time = input->Ascent_to_Altitude_Hours * 60.0;
    double Time_at_Altitude_Before_Dive = input->Hours_at_Altitude_Before_Dive * 60.0;

    if (dive->Diver_Acclimatized == TRUE) {
        dive->Barometric_Pressure = vpmb_calc_barometric_pressure(dive->Altitude_of_Dive, dive->units_fsw);

        for(i = 0; i < 16; i++) {
            dive->Adjusted_Critical_Radius_N2[i] = dive->Initial_Critical_Radius_N2[i];
            dive->Adjusted_Critical_Radius_He[i] = dive->Initial_Critical_Radius_He[i];
            dive->Helium_Pressure[i] = 0.0;
            dive->Nitrogen_Pressure[i] = (dive->Barometric_Pressure - dive->Water_Vapor_Pressure) * fraction_inert_gas;
        }
    } else {
        double Starting_Ambient_Pressure;
        double Ending_Ambient_Pressure;
        double Initial_Inspired_N2_Pressure;
        double Rate;
        double Nitrogen_Rate;
        double Inspired_Nitrogen_Pressure;

        if ((input->Starting_Acclimatized_Altitude >= dive->Altitude_of_Dive) || (input->Starting_Acclimatized_Altitude < 0.0)) {
            return BADALTITUDE;
        }

        dive->Barometric_Pressure = vpmb_calc_barometric_pressure(input->Starting_Acclimatized_Altitude, dive->units_fsw);

        Starting_Ambient_Pressure = dive->Barometric_Pressure;

        for(i = 0; i < 16; i++) {
            dive->Helium_Pressure[i] = 0.0;
            dive->Nitrogen_Pressure[i] = (dive->Barometric_Pressure - dive->Water_Vapor_Pressure) * fraction_inert_gas;
        }

        dive->Barometric_Pressure = vpmb_calc_barometric_pressure(dive->Altitude_of_Dive, dive->units_fsw);

        Ending_Ambient_Pressure = dive->Barometric_Pressure;
        Initial_Inspired_N2_Pressure = (Starting_Ambient_Pressure - dive->Water_Vapor_Pressure) * fraction_inert_gas;
        Rate = (Ending_Ambient_Pressure - Starting_Ambient_Pressure) / Ascent_to_Altitude_Time;
        Nitrogen_Rate = Rate * fraction_inert_gas;

        for (i=0; i < 16; i++) {
            double Initial_Nitrogen_Pressure = dive->Nitrogen_Pressure[i];

            dive->Nitrogen_Pressure[i] = vpmb_schreiner_equation(Initial_Inspired_N2_Pressure, Nitrogen_Rate, Ascent_to_Altitude_Time, dive->Nitrogen_Time_Constant[i], Initial_Nitrogen_Pressure);

            double Compartment_Gradient = (dive->Nitrogen_Pressure[i] + dive->Constant_Pressure_Other_Gases) - Ending_Ambient_Pressure;

            double Compartment_Gradient_Pascals = (Compartment_Gradient / dive->Units_Factor) * ATM;

            double Gradient_He_Bubble_Formation = ((2.0 * dive->Surface_Tension_Gamma * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma)) / (dive->Initial_Critical_Radius_He[i] * dive->Skin_Compression_GammaC));

            if (Compartment_Gradient_Pascals > Gradient_He_Bubble_Formation) {

                double New_Critical_Radius_He = ((2.0 * dive->Surface_Tension_Gamma * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma))) / (Compartment_Gradient_Pascals * dive->Skin_Compression_GammaC);

                dive->Adjusted_Critical_Radius_He[i] = dive->Initial_Critical_Radius_He[i] + (dive->Initial_Critical_Radius_He[i] -  New_Critical_Radius_He) * exp(-Time_at_Altitude_Before_Dive / dive->Regeneration_Time_Constant);

                dive->Initial_Critical_Radius_He[i] = dive->Adjusted_Critical_Radius_He[i];
            } else {
                double Ending_Radius_He = 1.0/(Compartment_Gradient_Pascals/ (2.0 * (dive->Surface_Tension_Gamma - dive->Skin_Compression_GammaC))  + 1.0/ dive->Initial_Critical_Radius_He[i]);

                double Regenerated_Radius_He = dive->Initial_Critical_Radius_He[i] + (Ending_Radius_He - dive->Initial_Critical_Radius_He[i]) * exp(-Time_at_Altitude_Before_Dive / dive->Regeneration_Time_Constant);

                dive->Initial_Critical_Radius_He[i] = Regenerated_Radius_He;
                dive->Adjusted_Critical_Radius_He[i] = dive->Initial_Critical_Radius_He[i];
            }
            double Gradient_N2_Bubble_Formation = ((2.0 * dive->Surface_Tension_Gamma * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma)) / (dive->Initial_Critical_Radius_N2[i] * dive->Skin_Compression_GammaC));

            if (Compartment_Gradient_Pascals > Gradient_N2_Bubble_Formation) {

                double New_Critical_Radius_N2 = ((2.0 * dive->Surface_Tension_Gamma * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma))) / (Compartment_Gradient_Pascals * dive->Skin_Compression_GammaC);

                dive->Adjusted_Critical_Radius_N2[i] = dive->Initial_Critical_Radius_N2[i] + (dive->Initial_Critical_Radius_N2[i]- New_Critical_Radius_N2) * exp(-Time_at_Altitude_Before_Dive / dive->Regeneration_Time_Constant);

                dive->Initial_Critical_Radius_N2[i] = dive->Adjusted_Critical_Radius_N2[i];
            } else {
                double Ending_Radius_N2 = 1.0/(Compartment_Gradient_Pascals/ (2.0*(dive->Surface_Tension_Gamma - dive->Skin_Compression_GammaC)) + 1.0/dive->Initial_Critical_Radius_N2[i]);

                double Regenerated_Radius_N2 = dive->Initial_Critical_Radius_N2[i] + (Ending_Radius_N2 - dive->Initial_Critical_Radius_N2[i]) * exp(-Time_at_Altitude_Before_Dive/ dive->Regeneration_Time_Constant);

                dive->Initial_Critical_Radius_N2[i] = Regenerated_Radius_N2;
                dive->Adjusted_Critical_Radius_N2[i] = dive->Initial_Critical_Radius_N2[i];
            }

        }
        Inspired_Nitrogen_Pressure = (dive->Barometric_Pressure - dive->Water_Vapor_Pressure) * fraction_inert_gas;

        for(i=0; i < 16; i++) {
            double Initial_Nitrogen_Pressure = dive->Nitrogen_Pressure[i];

            dive->Nitrogen_Pressure[i] = vpmb_haldane_equation(Initial_Nitrogen_Pressure, Inspired_Nitrogen_Pressure, dive->Nitrogen_Time_Constant[i], Time_at_Altitude_Before_Dive);
        }
    }
    return VALIDDATA;
}

/** \brief Updates the gas pressures when staying at a constant depth.
 *
 * Uses the ::HALDANE_EQUATION to update the gas pressures (dive_state.Helium_Pressure and dive_state.Nitrogen_Pressure)
 * for a time segment at constant depth.
 *
 * Side Effects: Sets
 * - dive_state.Ending_Ambient_Pressure
 * - dive_state.Helium_Pressure
 * - dive_state.Nitrogen_Pressure
 * - dive_state.Run_Time
 * - dive_state.Segment_Number
 * - dive_state.Segment_Time
 *
 */
void vpmb_gas_loadings_constant_depth(dive_state *dive, double Depth, double Run_Time_End_of_Segment)
{
    int i;
    double Ambient_Pressure = Depth + dive->Barometric_Pressure;

    dive->Segment_Time = Run_Time_End_of_Segment - dive->Run_Time;
    dive->Run_Time = Run_Time_End_of_Segment;
    dive->Segment_Number = dive->Segment_Number + 1;

    double Inspired_Helium_Pressure = (Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Helium[dive->Mix_Number - 1];

    double Inspired_Nitrogen_Pressure = (Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Nitrogen[dive->Mix_Number - 1];

    dive->Ending_Ambient_Pressure = Ambient_Pressure;

    double Temp_Helium_Pressure = 0.0;
    double Temp_Nitrogen_Pressure = 0.0;

    for(i = 0; i < 16; i++) {
        Temp_Helium_Pressure = dive->Helium_Pressure[i];
        Temp_Nitrogen_Pressure = dive->Nitrogen_Pressure[i];

        dive->Helium_Pressure[i] = vpmb_haldane_equation(Temp_Helium_Pressure, Inspired_Helium_Pressure, dive->Helium_Time_Constant[i], dive->Segment_Time);

        dive->Nitrogen_Pressure[i] = vpmb_haldane_equation(Temp_Nitrogen_Pressure, Inspired_Nitrogen_Pressure, dive->Nitrogen_Time_Constant[i], dive->Segment_Time);
    }
}

void vpmb_nuclear_regeneration(dive_state *dive, double Dive_Time)
{
    int i;

    for (i = 0; i < 16; i++) {
        double Crushing_Pressure_Pascals_He = (dive->Max_Crushing_Pressure_He[i] / dive->Units_Factor) * ATM;
        double Crushing_Pressure_Pascals_N2 = (dive->Max_Crushing_Pressure_N2[i]/ dive->Units_Factor) * ATM;

        double Ending_Radius_He = 1.0/(Crushing_Pressure_Pascals_He / (2.0 * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma)) + 1.0/dive->Adjusted_Critical_Radius_He[i]);

        double Ending_Radius_N2 = 1.0/(Crushing_Pressure_Pascals_N2/ (2.0 * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma)) + 1.0/ dive->Adjusted_Critical_Radius_N2[i]);

        dive->Regenerated_Radius_He[i] = dive->Adjusted_Critical_Radius_He[i] + (Ending_Radius_He - dive->Adjusted_Critical_Radius_He[i]) * exp(-Dive_Time / dive->Regeneration_Time_Constant);

        dive->Regenerated_Radius_N2[i] = dive->Adjusted_Critical_Radius_N2[i] + (Ending_Radius_N2 - dive->Adjusted_Critical_Radius_N2[i]) *exp(-Dive_Time / dive->Regeneration_Time_Constant);


        double Crush_Pressure_Adjust_Ratio_He = (Ending_Radius_He * (dive->Adjusted_Critical_Radius_He[i] - dive->Regenerated_Radius_He[i])) / (dive->Regenerated_Radius_He[i] * (dive->Adjusted_Critical_Radius_He[i] - Ending_Radius_He));

        double Crush_Pressure_Adjust_Ratio_N2 = (Ending_Radius_N2 * (dive->Adjusted_Critical_Radius_N2[i] - dive->Regenerated_Radius_N2[i])) / (dive->Regenerated_Radius_N2[i] * (dive->Adjusted_Critical_Radius_N2[i] - Ending_Radius_N2));

        double Adj_Crush_Pressure_He_Pascals = Crushing_Pressure_Pascals_He * Crush_Pressure_Adjust_Ratio_He;
        double Adj_Crush_Pressure_N2_Pascals = Crushing_Pressure_Pascals_N2 * Crush_Pressure_Adjust_Ratio_N2;

        dive->Adjusted_Crushing_Pressure_He[i] = (Adj_Crush_Pressure_He_Pascals / ATM) * dive->Units_Factor;
        dive->Adjusted_Crushing_Pressure_N2[i] = (Adj_Crush_Pressure_N2_Pascals / ATM) * dive->Units_Factor;
    }
}


double vpmb_crushing_pressure_helper(dive_state *dive, double Radius_Onset_of_Imperm_Molecule, double Ending_Ambient_Pressure_Pa, double Amb_Press_Onset_of_Imperm_Pa, double Gas_Tension_Onset_of_Imperm_Pa, double Gradient_Onset_of_Imperm_Pa)
{

    double A = Ending_Ambient_Pressure_Pa - Amb_Press_Onset_of_Imperm_Pa + Gas_Tension_Onset_of_Imperm_Pa + (2.0 * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma)) / Radius_Onset_of_Imperm_Molecule;
    double B = 2.0 * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma);
    double C = Gas_Tension_Onset_of_Imperm_Pa * pow(Radius_Onset_of_Imperm_Molecule, 3);

    double High_Bound = Radius_Onset_of_Imperm_Molecule;
    double Low_Bound = B / A;

    double Ending_Radius;
    double Crushing_Pressure_Pascals;

    if (vpmb_radius_root_finder(A, B, C, Low_Bound, High_Bound, &Ending_Radius) < 0) {
        vpmb_failure();
    }
    Crushing_Pressure_Pascals = Gradient_Onset_of_Imperm_Pa + Ending_Ambient_Pressure_Pa - Amb_Press_Onset_of_Imperm_Pa +  Gas_Tension_Onset_of_Imperm_Pa * (1.0 - pow(Radius_Onset_of_Imperm_Molecule, 3) / pow(Ending_Radius, 3));

    return (Crushing_Pressure_Pascals / ATM) * dive->Units_Factor;
}
int vpmb_onset_of_impermeability(dive_state *dive, double Starting_Ambient_Pressure, double Ending_Ambient_Pressure, double Rate, int i)
{
    int j;

    double Gradient_Onset_of_Imperm = dive->Gradient_Onset_of_Imperm_Atm * dive->Units_Factor;

    double Initial_Inspired_He_Pressure = (Starting_Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Helium[dive->Mix_Number - 1];

    double Initial_Inspired_N2_Pressure = (Starting_Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Nitrogen[dive->Mix_Number - 1];

    double Helium_Rate = Rate * dive->Fraction_Helium[dive->Mix_Number - 1];
    double Nitrogen_Rate = Rate * dive->Fraction_Nitrogen[dive->Mix_Number - 1];
    double Low_Bound = 0.0;

    double High_Bound = (Ending_Ambient_Pressure - Starting_Ambient_Pressure) / Rate;

    double Starting_Gas_Tension = dive->Initial_Helium_Pressure[i] + dive->Initial_Nitrogen_Pressure[i] + dive->Constant_Pressure_Other_Gases;

    double Function_at_Low_Bound = Starting_Ambient_Pressure - Starting_Gas_Tension - Gradient_Onset_of_Imperm;

    double High_Bound_Helium_Pressure = vpmb_schreiner_equation(Initial_Inspired_He_Pressure, Helium_Rate, High_Bound, dive->Helium_Time_Constant[i], dive->Initial_Helium_Pressure[i]);

    double High_Bound_Nitrogen_Pressure = vpmb_schreiner_equation(Initial_Inspired_N2_Pressure, Nitrogen_Rate, High_Bound, dive->Nitrogen_Time_Constant[i], dive->Initial_Nitrogen_Pressure[i]);

    double Ending_Gas_Tension = High_Bound_Helium_Pressure + High_Bound_Nitrogen_Pressure + dive->Constant_Pressure_Other_Gases;

    double Function_at_High_Bound = Ending_Ambient_Pressure -  Ending_Gas_Tension - Gradient_Onset_of_Imperm;

    double Time, Differential_Change;
    double Mid_Range_Ambient_Pressure, Gas_Tension_at_Mid_Range;

    if((Function_at_High_Bound * Function_at_Low_Bound) >= 0.0) {
        return  ROOTERROR;
    }


    if( Function_at_Low_Bound < 0.0) {
        Time = Low_Bound;
        Differential_Change = High_Bound - Low_Bound;
    } else {
        Time = High_Bound;
        Differential_Change = Low_Bound - High_Bound;
    }

    for(j = 0; j < 100; j++) {
        double Last_Diff_Change = Differential_Change;
        Differential_Change = Last_Diff_Change * 0.5;
        double Mid_Range_Time = Time + Differential_Change;

        Mid_Range_Ambient_Pressure = (Starting_Ambient_Pressure + Rate * Mid_Range_Time);

        double Mid_Range_Helium_Pressure = vpmb_schreiner_equation(Initial_Inspired_He_Pressure, Helium_Rate, Mid_Range_Time, dive->Helium_Time_Constant[i], dive->Initial_Helium_Pressure[i]);

        double Mid_Range_Nitrogen_Pressure = vpmb_schreiner_equation(Initial_Inspired_N2_Pressure, Nitrogen_Rate, Mid_Range_Time, dive->Nitrogen_Time_Constant[i], dive->Initial_Nitrogen_Pressure[i]);

        Gas_Tension_at_Mid_Range = Mid_Range_Helium_Pressure +  Mid_Range_Nitrogen_Pressure + dive->Constant_Pressure_Other_Gases;

        double Function_at_Mid_Range = Mid_Range_Ambient_Pressure - Gas_Tension_at_Mid_Range - Gradient_Onset_of_Imperm;

        if(Function_at_Mid_Range <= 0.0) {
            Time = Mid_Range_Time;
        }

        if ((fabs(Differential_Change) < 1.0E-3) || (Function_at_Mid_Range == 0.0)) {
            break;
        }

        if(j == 100) {
            return ROOTERROR;
        }
    }

    dive->Amb_Pressure_Onset_of_Imperm[i] = Mid_Range_Ambient_Pressure;
    dive->Gas_Tension_Onset_of_Imperm[i] = Gas_Tension_at_Mid_Range;

    return ROOTFOUND;
}


void vpmb_calc_crushing_pressure(dive_state *dive, double Starting_Depth, double Ending_Depth, double Rate)
{
    int i;

    double Gradient_Onset_of_Imperm = dive->Gradient_Onset_of_Imperm_Atm * dive->Units_Factor;
    double Gradient_Onset_of_Imperm_Pa = dive->Gradient_Onset_of_Imperm_Atm * ATM;

    double Starting_Ambient_Pressure = Starting_Depth + dive->Barometric_Pressure;
    double Ending_Ambient_Pressure = Ending_Depth + dive->Barometric_Pressure;

    for(i = 0; i < 16; i++) {
        double Starting_Gas_Tension = dive->Initial_Helium_Pressure[i] + dive->Initial_Nitrogen_Pressure[i] + dive->Constant_Pressure_Other_Gases;

        double Starting_Gradient = Starting_Ambient_Pressure - Starting_Gas_Tension;

        double Ending_Gas_Tension = dive->Helium_Pressure[i] + dive->Nitrogen_Pressure[i] + dive->Constant_Pressure_Other_Gases;

        double Ending_Gradient = Ending_Ambient_Pressure - Ending_Gas_Tension;

        double Radius_Onset_of_Imperm_He = 1.0 / (Gradient_Onset_of_Imperm_Pa / (2.0 * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma)) + 1.0 / dive->Adjusted_Critical_Radius_He[i]);

        double Radius_Onset_of_Imperm_N2 = 1.0/(Gradient_Onset_of_Imperm_Pa/ (2.0 * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma)) + 1.0 / dive->Adjusted_Critical_Radius_N2[i]);

        double Crushing_Pressure_He;
        double Crushing_Pressure_N2;
        if(Ending_Gradient <= Gradient_Onset_of_Imperm) {
            Crushing_Pressure_He = Ending_Ambient_Pressure - Ending_Gas_Tension;
            Crushing_Pressure_N2 = Ending_Ambient_Pressure - Ending_Gas_Tension;
        } else {
            double Ending_Ambient_Pressure_Pa;
            double Amb_Press_Onset_of_Imperm_Pa;
            double Gas_Tension_Onset_of_Imperm_Pa;

            if (Starting_Gradient == Gradient_Onset_of_Imperm) {
                dive->Amb_Pressure_Onset_of_Imperm[i] = Starting_Ambient_Pressure;
                dive->Gas_Tension_Onset_of_Imperm[i] = Starting_Gas_Tension;
            }

            if (Starting_Gradient < Gradient_Onset_of_Imperm) {
                vpmb_onset_of_impermeability(dive, Starting_Ambient_Pressure, Ending_Ambient_Pressure, Rate, i);
            }

            Ending_Ambient_Pressure_Pa = (Ending_Ambient_Pressure / dive->Units_Factor) * ATM;

            Amb_Press_Onset_of_Imperm_Pa = (dive->Amb_Pressure_Onset_of_Imperm[i] / dive->Units_Factor) * ATM;

            Gas_Tension_Onset_of_Imperm_Pa = (dive->Gas_Tension_Onset_of_Imperm[i] / dive->Units_Factor) * ATM;

            Crushing_Pressure_He = vpmb_crushing_pressure_helper(dive, Radius_Onset_of_Imperm_He, Ending_Ambient_Pressure_Pa, Amb_Press_Onset_of_Imperm_Pa, Gas_Tension_Onset_of_Imperm_Pa, Gradient_Onset_of_Imperm_Pa);

            Crushing_Pressure_N2 = vpmb_crushing_pressure_helper(dive, Radius_Onset_of_Imperm_N2, Ending_Ambient_Pressure_Pa, Amb_Press_Onset_of_Imperm_Pa, Gas_Tension_Onset_of_Imperm_Pa, Gradient_Onset_of_Imperm_Pa);

        }

        dive->Max_Crushing_Pressure_He[i] = max(dive->Max_Crushing_Pressure_He[i], Crushing_Pressure_He);
        dive->Max_Crushing_Pressure_N2[i] = max(dive->Max_Crushing_Pressure_N2[i], Crushing_Pressure_N2);
    }
}

void vpmb_gas_loadings_ascent_descent(dive_state *dive, double Starting_Depth, double Ending_Depth, double Rate)
{
    int i;
    double Starting_Ambient_Pressure;
    double Initial_Inspired_He_Pressure;
    double Initial_Inspired_N2_Pressure;
    double Helium_Rate;
    double Nitrogen_Rate;

    dive->Segment_Time = (Ending_Depth - Starting_Depth) / Rate;

    dive->Run_Time = dive->Run_Time + dive->Segment_Time;
    dive->Segment_Number = dive->Segment_Number + 1;
    dive->Ending_Ambient_Pressure = Ending_Depth + dive->Barometric_Pressure;

    Starting_Ambient_Pressure = Starting_Depth + dive->Barometric_Pressure;
    Initial_Inspired_He_Pressure = (Starting_Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Helium[dive->Mix_Number - 1];

    Initial_Inspired_N2_Pressure = (Starting_Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Nitrogen[dive->Mix_Number - 1];
    Helium_Rate = Rate * dive->Fraction_Helium[dive->Mix_Number - 1];
    Nitrogen_Rate = Rate * dive->Fraction_Nitrogen[dive->Mix_Number - 1];

    for(i=0; i < 16; i++) {
        dive->Initial_Helium_Pressure[i] = dive->Helium_Pressure[i];
        dive->Initial_Nitrogen_Pressure[i] = dive->Nitrogen_Pressure[i];

        dive->Helium_Pressure[i] = vpmb_schreiner_equation(Initial_Inspired_He_Pressure, Helium_Rate, dive->Segment_Time, dive->Helium_Time_Constant[i], dive->Initial_Helium_Pressure[i]);

        dive->Nitrogen_Pressure[i] = vpmb_schreiner_equation(Initial_Inspired_N2_Pressure, Nitrogen_Rate, dive->Segment_Time, dive->Nitrogen_Time_Constant[i], dive->Initial_Nitrogen_Pressure[i]);
    }
}

int vpmb_decompression_stop(dive_state *dive, double Deco_Stop_Depth, double Step_Size)
{
    int i;
    double Deco_Ceiling_Depth = 0.0;
    double Last_Run_Time = dive->Run_Time;
    double Round_Up_Operation = round((Last_Run_Time / dive->Minimum_Deco_Stop_Time) + 0.5) * dive->Minimum_Deco_Stop_Time;
    dive->Segment_Time = Round_Up_Operation - dive->Run_Time;
    dive->Run_Time = Round_Up_Operation;
    double Temp_Segment_Time = dive->Segment_Time;
    dive->Segment_Number = dive->Segment_Number + 1;
    double Ambient_Pressure = Deco_Stop_Depth + dive->Barometric_Pressure;
    dive->Ending_Ambient_Pressure = Ambient_Pressure;
    double Next_Stop = Deco_Stop_Depth - Step_Size;

    double Inspired_Helium_Pressure = (Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Helium[dive->Mix_Number - 1];

    double Inspired_Nitrogen_Pressure = (Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Nitrogen[dive->Mix_Number - 1];


    for(i=0; i <16; i++) {
        if((Inspired_Helium_Pressure + Inspired_Nitrogen_Pressure) > 0.0) {
            double Weighted_Allowable_Gradient = (dive->Deco_Gradient_He[i] * Inspired_Helium_Pressure + dive->Deco_Gradient_N2[i]* Inspired_Nitrogen_Pressure) / (Inspired_Helium_Pressure + Inspired_Nitrogen_Pressure);

            if ((Inspired_Helium_Pressure + Inspired_Nitrogen_Pressure + dive->Constant_Pressure_Other_Gases - Weighted_Allowable_Gradient) > (Next_Stop + dive->Barometric_Pressure)) {
                return OFFGASSINGERROR;
            }
        }
    }
    while(TRUE) {
        for(i=0; i < 16; i++) {
            double Initial_Helium_Pressure = dive->Helium_Pressure[i];
            double Initial_Nitrogen_Pressure = dive->Nitrogen_Pressure[i];

            dive->Helium_Pressure[i] = vpmb_haldane_equation(Initial_Helium_Pressure, Inspired_Helium_Pressure, dive->Helium_Time_Constant[i], dive->Segment_Time);

            dive->Nitrogen_Pressure[i] = vpmb_haldane_equation(Initial_Nitrogen_Pressure, Inspired_Nitrogen_Pressure, dive->Nitrogen_Time_Constant[i], dive->Segment_Time);
        }
        Deco_Ceiling_Depth = vpmb_calc_deco_ceiling(dive);
        if (Deco_Ceiling_Depth > Next_Stop) {
            dive->Segment_Time = dive->Minimum_Deco_Stop_Time;
            double Time_Counter = Temp_Segment_Time;
            Temp_Segment_Time =  Time_Counter + dive->Minimum_Deco_Stop_Time;
            double Last_Run_Time = dive->Run_Time;
            dive->Run_Time = Last_Run_Time + dive->Minimum_Deco_Stop_Time;
            continue;
        }
        break;
    }
    dive->Segment_Time = Temp_Segment_Time;
    return ALLGOOD;
}

double vpmb_calculate_deco_gradient(dive_state *dive, double Allowable_Gradient_Molecule, double Amb_Press_First_Stop_Pascals, double Amb_Press_Next_Stop_Pascals)
{

    double Allow_Grad_First_Stop_Pa = (Allowable_Gradient_Molecule / dive->Units_Factor) * ATM;
    double Radius_First_Stop = (2.0 * dive->Surface_Tension_Gamma) / Allow_Grad_First_Stop_Pa;

    double A = Amb_Press_Next_Stop_Pascals;
    double B = -2.0 * dive->Surface_Tension_Gamma;
    double C = (Amb_Press_First_Stop_Pascals + (2.0 * dive->Surface_Tension_Gamma)/ Radius_First_Stop) * Radius_First_Stop * (Radius_First_Stop * (Radius_First_Stop));

    double Low_Bound = Radius_First_Stop;
    double High_Bound = Radius_First_Stop * pow((Amb_Press_First_Stop_Pascals / Amb_Press_Next_Stop_Pascals), (1.0/3.0));

    double Ending_Radius;
    double Deco_Gradient_Pascals;
    if(vpmb_radius_root_finder(A, B, C, Low_Bound, High_Bound, &Ending_Radius) < 0) {
        vpmb_failure();
    }

    Deco_Gradient_Pascals = (2.0 * dive->Surface_Tension_Gamma) / Ending_Radius;
    return (Deco_Gradient_Pascals / ATM) * dive->Units_Factor;
}

void vpmb_boyles_law_compensation(dive_state *dive, double First_Stop_Depth, double Deco_Stop_Depth, double Step_Size)
{

    int i;
    double Next_Stop = Deco_Stop_Depth - Step_Size;
    double Ambient_Pressure_First_Stop = First_Stop_Depth + dive->Barometric_Pressure;
    double Ambient_Pressure_Next_Stop = Next_Stop + dive->Barometric_Pressure;

    double Amb_Press_First_Stop_Pascals = (Ambient_Pressure_First_Stop / dive->Units_Factor) * ATM;
    double Amb_Press_Next_Stop_Pascals = (Ambient_Pressure_Next_Stop / dive->Units_Factor) * ATM;

    for(i=0; i < 16; i++) {
        dive->Deco_Gradient_He[i] = vpmb_calculate_deco_gradient(dive, dive->Allowable_Gradient_He[i], Amb_Press_First_Stop_Pascals, Amb_Press_Next_Stop_Pascals);

        dive->Deco_Gradient_N2[i] = vpmb_calculate_deco_gradient(dive, dive->Allowable_Gradient_N2[i], Amb_Press_First_Stop_Pascals, Amb_Press_Next_Stop_Pascals);
    }
}

void vpmb_calc_max_actual_gradient(dive_state *dive, double Deco_Stop_Depth)
{
    int i;
    for(i = 0; i <16; i++) {
        double Compartment_Gradient = (dive->Helium_Pressure[i] + dive->Nitrogen_Pressure[i] + dive->Constant_Pressure_Other_Gases) - (Deco_Stop_Depth + dive->Barometric_Pressure);
        if( Compartment_Gradient <= 0.0) {
            Compartment_Gradient = 0.0;
        }

        dive->Max_Actual_Gradient[i] = max(dive->Max_Actual_Gradient[i], Compartment_Gradient);
    }
}

void vpmb_projected_ascent(dive_state *dive, double Starting_Depth, double Rate, double Step_Size)
{
    int i, j;

    double New_Ambient_Pressure = dive->Deco_Stop_Depth + dive->Barometric_Pressure;
    double Starting_Ambient_Pressure = Starting_Depth + dive->Barometric_Pressure;

    double Initial_Inspired_He_Pressure = (Starting_Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Helium[dive->Mix_Number - 1];

    double Initial_Inspired_N2_Pressure = (Starting_Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Nitrogen[dive->Mix_Number - 1];

    double Helium_Rate = Rate * dive->Fraction_Helium[dive->Mix_Number - 1];
    double Nitrogen_Rate = Rate * dive->Fraction_Nitrogen[dive->Mix_Number - 1];

    double Temp_Gas_Loading[16];
    double Allowable_Gas_Loading[16];
    double Initial_Helium_Pressure[16];
    double Initial_Nitrogen_Pressure[16];

    for(i=0; i < 16; i++) {
        Initial_Helium_Pressure[i] = dive->Helium_Pressure[i];
        Initial_Nitrogen_Pressure[i] = dive->Nitrogen_Pressure[i];
    }

    while(TRUE) {
        double Ending_Ambient_Pressure = New_Ambient_Pressure;

        double Segment_Time = (Ending_Ambient_Pressure - Starting_Ambient_Pressure) / Rate;

        for(i=0; i < 16; i++) {
            double Temp_Helium_Pressure = vpmb_schreiner_equation(Initial_Inspired_He_Pressure, Helium_Rate, Segment_Time, dive->Helium_Time_Constant[i], Initial_Helium_Pressure[i]);

            double Temp_Nitrogen_Pressure = vpmb_schreiner_equation(Initial_Inspired_N2_Pressure, Nitrogen_Rate, Segment_Time, dive->Nitrogen_Time_Constant[i], Initial_Nitrogen_Pressure[i]);
            double Weighted_Allowable_Gradient;

            Temp_Gas_Loading[i] = Temp_Helium_Pressure + Temp_Nitrogen_Pressure;
            if (Temp_Gas_Loading[i] > 0.0) {
                Weighted_Allowable_Gradient = (dive->Allowable_Gradient_He[i] * Temp_Helium_Pressure + dive->Allowable_Gradient_N2[i] * Temp_Nitrogen_Pressure) / Temp_Gas_Loading[i];
            } else {
                Weighted_Allowable_Gradient = min(dive->Allowable_Gradient_He[i], dive->Allowable_Gradient_N2[i]);
            }

            Allowable_Gas_Loading[i] = Ending_Ambient_Pressure + Weighted_Allowable_Gradient - dive->Constant_Pressure_Other_Gases;
        }

        BOOL end_sub = TRUE;
        for (j=0; j < 16; j++) {
            if(Temp_Gas_Loading[j] > Allowable_Gas_Loading[j]) {
                New_Ambient_Pressure = Ending_Ambient_Pressure + Step_Size;
                dive->Deco_Stop_Depth = dive->Deco_Stop_Depth + Step_Size;
                end_sub = FALSE;
                break;
            }
        }
        if(end_sub != TRUE) {
            continue;
        } else {
            break;
        }
    }
}

void vpmb_calc_ascent_ceiling(dive_state *dive)
{
    int i;
    double Gas_Loading = 0.0;
    double Compartment_Ascent_Ceiling[16] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

    for(i=0; i < 16; i++) {
        double Weighted_Allowable_Gradient, Tolerated_Ambient_Pressure;

        Gas_Loading = dive->Helium_Pressure[i] + dive->Nitrogen_Pressure[i];

        if(Gas_Loading > 0.0) {
            Weighted_Allowable_Gradient = (dive->Allowable_Gradient_He[i] * dive->Helium_Pressure[i] +  dive->Allowable_Gradient_N2[i] * dive->Nitrogen_Pressure[i]) / (dive->Helium_Pressure[i] + dive->Nitrogen_Pressure[i]);
            Tolerated_Ambient_Pressure = (Gas_Loading + dive->Constant_Pressure_Other_Gases) - Weighted_Allowable_Gradient;
        } else {
            Weighted_Allowable_Gradient = min(dive->Allowable_Gradient_He[i], dive->Allowable_Gradient_N2[i]);
            Tolerated_Ambient_Pressure = dive->Constant_Pressure_Other_Gases - Weighted_Allowable_Gradient;
        }
        if (Tolerated_Ambient_Pressure < 0.0) {
            Tolerated_Ambient_Pressure = 0.0;
        }

        Compartment_Ascent_Ceiling[i] = Tolerated_Ambient_Pressure - dive->Barometric_Pressure;
    }


    dive->Ascent_Ceiling_Depth = Compartment_Ascent_Ceiling[0];

    for(i=1; i <16; i++) {
        dive->Ascent_Ceiling_Depth = max(dive->Ascent_Ceiling_Depth, Compartment_Ascent_Ceiling[i]);
    }
}


void vpmb_calc_surface_phase_volume_time(dive_state *dive)
{
    int i;
    double Surface_Inspired_N2_Pressure = (dive->Barometric_Pressure - dive->Water_Vapor_Pressure) * fraction_inert_gas;
    for(i=0; i <16; i++) {
        if(dive->Nitrogen_Pressure[i] > Surface_Inspired_N2_Pressure) {
            dive->Surface_Phase_Volume_Time[i]= (dive->Helium_Pressure[i] / dive->Helium_Time_Constant[i] + (dive->Nitrogen_Pressure[i] - Surface_Inspired_N2_Pressure) / dive->Nitrogen_Time_Constant[i]) /(dive->Helium_Pressure[i] + dive->Nitrogen_Pressure[i] - Surface_Inspired_N2_Pressure);
        }

        else if((dive->Nitrogen_Pressure[i] <= Surface_Inspired_N2_Pressure) && (dive->Helium_Pressure[i] + dive->Nitrogen_Pressure[i] >= Surface_Inspired_N2_Pressure)) {
            double Decay_Time_to_Zero_Gradient = 1.0 / (dive->Nitrogen_Time_Constant[i] - dive->Helium_Time_Constant[i]) * log((Surface_Inspired_N2_Pressure - dive->Nitrogen_Pressure[i])/dive->Helium_Pressure[i]);

            double Integral_Gradient_x_Time = dive->Helium_Pressure[i] / dive->Helium_Time_Constant[i] * (1.0 - exp(-dive->Helium_Time_Constant[i] * Decay_Time_to_Zero_Gradient)) + (dive->Nitrogen_Pressure[i] - Surface_Inspired_N2_Pressure) / dive->Nitrogen_Time_Constant[i] * (1.0 - exp(-dive->Nitrogen_Time_Constant[i] * Decay_Time_to_Zero_Gradient));

            dive->Surface_Phase_Volume_Time[i] = Integral_Gradient_x_Time/(dive->Helium_Pressure[i] + dive->Nitrogen_Pressure[i] - Surface_Inspired_N2_Pressure);
        }

        else {
            dive->Surface_Phase_Volume_Time[i] = 0.0;
        }
    }
}

void vpmb_critical_volume(dive_state *dive, double Deco_Phase_Volume_Time)
{

    int i;

    double Parameter_Lambda_Pascals = (dive->Crit_Volume_Parameter_Lambda / 33.0) * ATM;
    double B, C;

    for(i=0; i <16; i++) {
        double Phase_Volume_Time = Deco_Phase_Volume_Time +  dive->Surface_Phase_Volume_Time[i];
        double Adj_Crush_Pressure_He_Pascals = (dive->Adjusted_Crushing_Pressure_He[i]/dive->Units_Factor) * ATM;

        double Initial_Allowable_Grad_He_Pa = (dive->Initial_Allowable_Gradient_He[i]/dive->Units_Factor) * ATM;

        B = Initial_Allowable_Grad_He_Pa + (Parameter_Lambda_Pascals * dive->Surface_Tension_Gamma) / (dive->Skin_Compression_GammaC * Phase_Volume_Time);

        C = (dive->Surface_Tension_Gamma * (dive->Surface_Tension_Gamma * (Parameter_Lambda_Pascals * Adj_Crush_Pressure_He_Pascals))) /(dive->Skin_Compression_GammaC * (dive->Skin_Compression_GammaC * Phase_Volume_Time));

        double New_Allowable_Grad_He_Pascals = (B + sqrt(pow(B, 2) - 4.0 * C))/2.0;

        dive->Allowable_Gradient_He[i] = (New_Allowable_Grad_He_Pascals / ATM) * dive->Units_Factor;


        double Adj_Crush_Pressure_N2_Pascals = (dive->Adjusted_Crushing_Pressure_N2[i] / dive->Units_Factor) * ATM;

        double Initial_Allowable_Grad_N2_Pa = (dive->Initial_Allowable_Gradient_N2[i] / dive->Units_Factor) * ATM;

        B = Initial_Allowable_Grad_N2_Pa + (Parameter_Lambda_Pascals * dive->Surface_Tension_Gamma)/ (dive->Skin_Compression_GammaC * Phase_Volume_Time);

        C = (dive->Surface_Tension_Gamma * (dive->Surface_Tension_Gamma * (Parameter_Lambda_Pascals * Adj_Crush_Pressure_N2_Pascals))) /(dive->Skin_Compression_GammaC * (dive->Skin_Compression_GammaC * Phase_Volume_Time));

        double New_Allowable_Grad_N2_Pascals = (B + sqrt(pow(B, 2) - 4.0 * C))/2.0;

        dive->Allowable_Gradient_N2[i] = (New_Allowable_Grad_N2_Pascals / ATM) * dive->Units_Factor;
    }
}

int vpmb_calc_start_of_deco_zone(dive_state *dive, double Starting_Depth, double Rate)
{
    int i, j;

    dive->Depth_Start_of_Deco_Zone = 0.0;
    double Starting_Ambient_Pressure = Starting_Depth + dive->Barometric_Pressure;

    double Initial_Inspired_He_Pressure = (Starting_Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Helium[dive->Mix_Number - 1];

    double Initial_Inspired_N2_Pressure = (Starting_Ambient_Pressure - dive->Water_Vapor_Pressure) * dive->Fraction_Nitrogen[dive->Mix_Number - 1];

    double Helium_Rate = Rate * dive->Fraction_Helium[dive->Mix_Number - 1];
    double Nitrogen_Rate = Rate * dive->Fraction_Nitrogen[dive->Mix_Number - 1];

    double Low_Bound = 0.0;
    double High_Bound = -1.0 * (Starting_Ambient_Pressure / Rate);

    for(i = 0; i < 16; i++) {
        double Initial_Helium_Pressure = dive->Helium_Pressure[i];
        double Initial_Nitrogen_Pressure = dive->Nitrogen_Pressure[i];

        double Function_at_Low_Bound = Initial_Helium_Pressure + Initial_Nitrogen_Pressure + dive->Constant_Pressure_Other_Gases - Starting_Ambient_Pressure;

        double High_Bound_Helium_Pressure = vpmb_schreiner_equation(Initial_Inspired_He_Pressure, Helium_Rate, High_Bound, dive->Helium_Time_Constant[i], Initial_Helium_Pressure);

        double High_Bound_Nitrogen_Pressure = vpmb_schreiner_equation(Initial_Inspired_N2_Pressure, Nitrogen_Rate, High_Bound, dive->Nitrogen_Time_Constant[i], Initial_Nitrogen_Pressure);

        double Function_at_High_Bound = High_Bound_Helium_Pressure + High_Bound_Nitrogen_Pressure + dive->Constant_Pressure_Other_Gases;

        if((Function_at_High_Bound * Function_at_Low_Bound) >= 0.0) {
            return ROOTERROR;
        }

        double Time_to_Start_of_Deco_Zone, Differential_Change;
        if (Function_at_Low_Bound < 0.0) {
            Time_to_Start_of_Deco_Zone = Low_Bound;
            Differential_Change = High_Bound - Low_Bound;
        } else {
            Time_to_Start_of_Deco_Zone = High_Bound;
            Differential_Change = Low_Bound - High_Bound;
        }
        double Last_Diff_Change;
        for(j = 0; j < 100; j++) {
            Last_Diff_Change = Differential_Change;
            Differential_Change = Last_Diff_Change * 0.5;

            double Mid_Range_Time = Time_to_Start_of_Deco_Zone + Differential_Change;

            double Mid_Range_Helium_Pressure = vpmb_schreiner_equation(Initial_Inspired_He_Pressure, Helium_Rate, Mid_Range_Time, dive->Helium_Time_Constant[i], Initial_Helium_Pressure);

            double Mid_Range_Nitrogen_Pressure = vpmb_schreiner_equation(Initial_Inspired_N2_Pressure, Nitrogen_Rate, Mid_Range_Time, dive->Nitrogen_Time_Constant[i], Initial_Nitrogen_Pressure);

            double Function_at_Mid_Range = Mid_Range_Helium_Pressure + Mid_Range_Nitrogen_Pressure + dive->Constant_Pressure_Other_Gases - (Starting_Ambient_Pressure + Rate * Mid_Range_Time);

            if (Function_at_Mid_Range <= 0.0) {
                Time_to_Start_of_Deco_Zone = Mid_Range_Time;
            }

            if ((fabs(Differential_Change) < 1.0E-3) || (Function_at_Mid_Range == 0.0)) {
                break;
            }

            if(j == 100) {
                return ROOTERROR;
            }
        }

        double Cpt_Depth_Start_of_Deco_Zone = (Starting_Ambient_Pressure + Rate * Time_to_Start_of_Deco_Zone) - dive->Barometric_Pressure;

        dive->Depth_Start_of_Deco_Zone = max(dive->Depth_Start_of_Deco_Zone, Cpt_Depth_Start_of_Deco_Zone);
    }
    return ALLGOOD;
}

void vpmb_calc_initial_allowable_gradient(dive_state *dive)
{
    int i;
    for(i=0; i < 16; i++) {
        double Initial_Allowable_Grad_N2_Pa = ((2.0 * dive->Surface_Tension_Gamma * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma)) / (dive->Regenerated_Radius_N2[i] * dive->Skin_Compression_GammaC));

        double Initial_Allowable_Grad_He_Pa = ((2.0 * dive->Surface_Tension_Gamma * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma)) / (dive->Regenerated_Radius_He[i] * dive->Skin_Compression_GammaC));

        dive->Initial_Allowable_Gradient_N2[i] = (Initial_Allowable_Grad_N2_Pa / ATM) * dive->Units_Factor;
        dive->Initial_Allowable_Gradient_He[i] = (Initial_Allowable_Grad_He_Pa / ATM) * dive->Units_Factor;

        dive->Allowable_Gradient_He[i] = dive->Initial_Allowable_Gradient_He[i];
        dive->Allowable_Gradient_N2[i] = dive->Initial_Allowable_Gradient_N2[i];
    }
}


void vpmb_gas_loadings_surface_interval(dive_state *dive, double Surface_Interval_Time)
{
    int i;
    double Inspired_Helium_Pressure = 0.0;
    double Inspired_Nitrogen_Pressure = (dive->Barometric_Pressure - dive->Water_Vapor_Pressure) * fraction_inert_gas;

    for (i=0; i < 16; i++) {
        double Temp_Helium_Pressure = dive->Helium_Pressure[i];
        double Temp_Nitrogen_Pressure = dive->Nitrogen_Pressure[i];

        dive->Helium_Pressure[i] = vpmb_haldane_equation(Temp_Helium_Pressure, Inspired_Helium_Pressure, dive->Helium_Time_Constant[i], Surface_Interval_Time);

        dive->Nitrogen_Pressure[i] = vpmb_haldane_equation(Temp_Nitrogen_Pressure, Inspired_Nitrogen_Pressure, dive->Nitrogen_Time_Constant[i], Surface_Interval_Time);
    }
}

double vpmb_new_critical_radius(dive_state *dive, double Max_Actual_Gradient_Pascals, double Adj_Crush_Pressure_Pascals)
{
    return ((2.0 * dive->Surface_Tension_Gamma * (dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma))) / (Max_Actual_Gradient_Pascals * dive->Skin_Compression_GammaC - dive->Surface_Tension_Gamma * Adj_Crush_Pressure_Pascals);
}

void vpmb_vpm_repetitive_algorithm(dive_state *dive, double Surface_Interval_Time)
{
    int i;
    for (i=0; i <16; i++) {
        double Max_Actual_Gradient_Pascals = (dive->Max_Actual_Gradient[i] / dive->Units_Factor) * ATM;

        double Adj_Crush_Pressure_He_Pascals = (dive->Adjusted_Crushing_Pressure_He[i] / dive->Units_Factor) * ATM;
        double Adj_Crush_Pressure_N2_Pascals = (dive->Adjusted_Crushing_Pressure_N2[i] / dive->Units_Factor) * ATM;

        if (dive->Max_Actual_Gradient[i] > dive->Initial_Allowable_Gradient_N2[i]) {
            double New_Critical_Radius_N2 = vpmb_new_critical_radius(dive, Max_Actual_Gradient_Pascals, Adj_Crush_Pressure_N2_Pascals);

            dive->Adjusted_Critical_Radius_N2[i] = dive->Initial_Critical_Radius_N2[i] + (dive->Initial_Critical_Radius_N2[i] - New_Critical_Radius_N2) * exp(-Surface_Interval_Time / dive->Regeneration_Time_Constant);
        }

        else {
            dive->Adjusted_Critical_Radius_N2[i] = dive->Initial_Critical_Radius_N2[i];
        }

        if( dive->Max_Actual_Gradient[i] > dive->Initial_Allowable_Gradient_He[i]) {
            double New_Critical_Radius_He = vpmb_new_critical_radius(dive, Max_Actual_Gradient_Pascals, Adj_Crush_Pressure_He_Pascals);

            dive->Adjusted_Critical_Radius_He[i] = dive->Initial_Critical_Radius_He[i] + (dive->Initial_Critical_Radius_He[i] - New_Critical_Radius_He) * exp(-Surface_Interval_Time / dive->Regeneration_Time_Constant);
        } else {
            dive->Adjusted_Critical_Radius_He[i] = dive->Initial_Critical_Radius_He[i];
        }

    }
}

int vpmb_validate_data(json_input *input, dive_state *dive)
{

    lowercase_string(input->Units);

    if (strcmp(input->Units, "fsw") == 0) {
        dive->units_fsw = TRUE;
    } else if (strcmp(input->Units, "msw") == 0) {
        dive->units_fsw = FALSE;
    } else {
        return INVALIDDATA;
    }

    if (input->Regeneration_Time_Constant <= 0) {
        return INVALIDDATA;
    }

    if ((dive->units_fsw == TRUE) && (input->Altitude_of_Dive > 30000.0)) {
        return INVALIDDATA;
    }

    if ((dive->units_fsw == FALSE) && (input->Altitude_of_Dive > 9144.0)) {
        return INVALIDDATA;
    }

    lowercase_string(input->Diver_Acclimatized_at_Altitude);

    if (strcmp(input->Diver_Acclimatized_at_Altitude, "yes") == 0) {
        dive->Diver_Acclimatized = TRUE;
    } else if (strcmp(input->Diver_Acclimatized_at_Altitude, "no") == 0) {
        dive->Diver_Acclimatized = FALSE;
    } else {
        return INVALIDDATA;
    }

    dive->Critical_Radius_N2_Microns = input->Critical_Radius_N2_Microns;
    dive->Critical_Radius_He_Microns = input->Critical_Radius_He_Microns;

    /* nitrogen */
    if ((input->Critical_Radius_N2_Microns < 0.2) || (input->Critical_Radius_N2_Microns) > 1.35) {
        return INVALIDDATA;
    }

    /* helium */
    if ((input->Critical_Radius_He_Microns < 0.2) || (input->Critical_Radius_He_Microns > 1.35)) {
        return INVALIDDATA;
    }

    return VALIDDATA;
}

int vpmb_initialize_data(json_input *input, dive_state *dive)
{
    int i;

    dive->Depth = 0;
    dive->Decompressing = FALSE;
    dive->decomp_stops=NULL;
    dive->Start_of_Decompression_Zone = 0;
    dive->Wait_Time = 0;
    /* This is to make sure the values are initialized before they get used for anything */
    dive->Fraction_Helium = NULL;
    dive->Fraction_Nitrogen = NULL;
    strlcpy(dive->Units, "" , sizeof(dive->Units));
    strlcpy(dive->Units_Word1, "" , sizeof(dive->Units_Word1));
    strlcpy(dive->Units_Word2, "" , sizeof(dive->Units_Word2));

    dive->Surface_Tension_Gamma = input->Surface_Tension_Gamma;
    dive->Skin_Compression_GammaC = input->Skin_Compression_GammaC;
    dive->Crit_Volume_Parameter_Lambda = input->Crit_Volume_Parameter_Lambda;
    dive->Gradient_Onset_of_Imperm_Atm = input->Gradient_Onset_of_Imperm_Atm;
    dive->Minimum_Deco_Stop_Time = input->Minimum_Deco_Stop_Time;
    dive->Critical_Radius_N2_Microns = input->Critical_Radius_N2_Microns;
    dive->Critical_Radius_He_Microns = input->Critical_Radius_He_Microns;
    dive->Regeneration_Time_Constant = input->Regeneration_Time_Constant;
    dive->Minimum_Deco_Stop_Time = input->Minimum_Deco_Stop_Time;

    /* INITIALIZE CONSTANTS/VARIABLES BASED ON SELECTION OF UNITS - FSW OR MSW */
    /* fsw = feet of seawater, a unit of pressure */
    /* msw = meters of seawater, a unit of pressure */

    if(dive->units_fsw == TRUE) {
        strlcpy(dive->Units_Word1,"fswg", strlen(dive->Units_Word1));
        strlcpy(dive->Units_Word2,"fsw/min", strlen(dive->Units_Word2));
        dive->Units_Factor = 33.0;
        dive->Water_Vapor_Pressure = 1.607; /* based on respiratory quotient of 0.8  (Schreiner value) */
    } else {
        size_t t= strlen(dive->Units_Word1);
        strlcpy(dive->Units_Word1, "mswg", t);
        strlcpy(dive->Units_Word2, "msw/min", strlen(dive->Units_Word2));
        dive->Units_Factor = 10.1325;
        dive->Water_Vapor_Pressure = 0.493; /* based on respiratory quotient of 0.8  (Schreiner value) */
    }

    /* INITIALIZE CONSTANTS/VARIABLES */
    dive->Constant_Pressure_Other_Gases = (input->Pressure_Other_Gases_mmHg / 760.0) * dive->Units_Factor;
    dive->Run_Time = 0.0;
    dive->Segment_Number = 0;
    dive->Last_Direction_Depth=dive->Depth;
    dive->Last_Direction_Time=dive->Run_Time;

    for (i = 0; i < Buhlmann_Compartments; i++) {
        dive->Helium_Time_Constant[i] = log(2.0) / Helium_Half_Time[i];
        dive->Nitrogen_Time_Constant[i] = log(2.0) / Nitrogen_Half_Time[i];
        dive->Max_Crushing_Pressure_He[i] = 0.0;
        dive->Max_Crushing_Pressure_N2[i] = 0.0;
        dive->Max_Actual_Gradient[i] = 0.0;
        dive->Surface_Phase_Volume_Time[i] = 0.0;
        dive->Amb_Pressure_Onset_of_Imperm[i] = 0.0;
        dive->Gas_Tension_Onset_of_Imperm[i] = 0.0;
        dive->Initial_Critical_Radius_N2[i] = input->Critical_Radius_N2_Microns * 1.0E-6;
        dive->Initial_Critical_Radius_He[i] = input->Critical_Radius_He_Microns * 1.0E-6;
    }

    lowercase_string(input->Critical_Volume_Algorithm);

    if (strcmp(input->Critical_Volume_Algorithm, "on") == 0) {
        dive->Critical_Volume_Algorithm_Off = FALSE;
    } else if (strcmp(input->Critical_Volume_Algorithm, "off") == 0) {
        dive->Critical_Volume_Algorithm_Off = TRUE;
    } else {
        return INVALIDDATA;
    }

    lowercase_string(input->Altitude_Dive_Algorithm);

    if (strcmp(input->Altitude_Dive_Algorithm, "on") == 0) {
        dive->Altitude_Dive_Algorithm_Off = FALSE;
        if ((input->Ascent_to_Altitude_Hours <= 0) && (dive->Diver_Acclimatized == FALSE)) {
            return INVALIDDATA;
        }
    } else if (strcmp(input->Altitude_Dive_Algorithm, "off") == 0) {
        dive->Altitude_Dive_Algorithm_Off = TRUE;
    } else {
        return INVALIDDATA;
    }

    /*INITIALIZE VARIABLES FOR SEA LEVEL OR ALTITUDE DIVE
      See subroutines for explanation of altitude calculations.  Purposes are
      1) to determine barometric pressure and 2) set or adjust the VPM critical
      radius variables and gas loadings, as applicable, based on altitude,
      ascent to altitude before the dive, and time at altitude before the dive*/

    if (dive->Altitude_Dive_Algorithm_Off == TRUE) {
        dive->Altitude_of_Dive = 0.0;
        dive->Barometric_Pressure = vpmb_calc_barometric_pressure(dive->Altitude_of_Dive, dive->units_fsw);

        for (i = 0 ; i < Buhlmann_Compartments; i++) {
            dive->Adjusted_Critical_Radius_N2[i] = dive->Initial_Critical_Radius_N2[i];
            dive->Adjusted_Critical_Radius_He[i] = dive->Initial_Critical_Radius_He[i];
            dive->Helium_Pressure[i] = 0.0;
            dive->Nitrogen_Pressure[i] = (dive->Barometric_Pressure - dive->Water_Vapor_Pressure) * fraction_inert_gas;
        }
    } else {
        vpmb_vpm_altitude_dive_algorithm(input, dive);
    }

    return VALIDDATA;
}


void vpmb_profile_code_loop(dive_state *dive, single_dive *current_dive)
{
    int i;

    for(i=0; i < current_dive->num_profile_codes; i++) {
        dive_profile *current_profile = &(current_dive->dive_profiles[i]);
        int Profile_Code = current_profile->profile_code;

        if(Profile_Code == 1) {
            char Word[10] = "TEMP";
            dive->Starting_Depth = current_profile->starting_depth;
            dive->Ending_Depth = current_profile->ending_depth;
            dive->Rate = current_profile->rate;
            dive->Mix_Number = current_profile->gasmix;

            vpmb_gas_loadings_ascent_descent(dive, dive->Starting_Depth, dive->Ending_Depth, dive->Rate);
            if(dive->Ending_Depth > dive->Starting_Depth) {
                vpmb_calc_crushing_pressure(dive, dive->Starting_Depth, dive->Ending_Depth, dive->Rate);
            }
            /* the error seems unnecessary */
            if( dive->Ending_Depth > dive->Starting_Depth) {
                strlcpy(Word, "Descent", strlen(Word));
            } else if (dive->Starting_Depth <  dive->Ending_Depth) {
                strlcpy(Word, "Ascent", strlen(Word));
            } else {
                strlcpy(Word, "ERROR", strlen(Word));
            }

            /* dive->output_object.add_dive_profile_entry_descent(dive->Segment_Number, dive->Segment_Time, dive->Run_Time, dive->Mix_Number, Word, dive->Starting_Depth, dive->Ending_Depth, dive->Rate) */
        } else if( Profile_Code == 2) {
            dive->Depth = current_profile->depth;
            dive->Run_Time_End_of_Segment = current_profile->run_time_at_end_of_segment;
            dive->Mix_Number = current_profile->gasmix;
            vpmb_gas_loadings_constant_depth(dive, dive->Depth, dive->Run_Time_End_of_Segment);

            /* dive->output_object.add_dive_profile_entry_ascent(dive->Segment_Number, dive->Segment_Time, dive->Run_Time, dive->Mix_Number, dive->Depth) */
        } else if(Profile_Code == 99) {
            break;
        }

    }
}

void vpmb_deco_stop_loop_block_within_critical_volume_loop(dive_state *dive)
{
    int i;
    while(TRUE) {
        vpmb_gas_loadings_ascent_descent(dive, dive->Starting_Depth, dive->Deco_Stop_Depth, dive->Rate);

        if(dive->Deco_Stop_Depth <= 0.0) {
            break;
        }

        if(dive->Number_of_Changes > 1) {
            for (i=1; i < dive->Number_of_Changes; i++) {
                if (dive->Depth_Change[i] >= dive->Deco_Stop_Depth) {
                    dive->Mix_Number = dive->Mix_Change[i];
                    dive->Rate = dive->Rate_Change[i];
                    dive->Step_Size = dive->Step_Size_Change[i];
                }
            }
        }
        vpmb_boyles_law_compensation(dive, dive->First_Stop_Depth, dive->Deco_Stop_Depth, dive->Step_Size);
        vpmb_decompression_stop(dive, dive->Deco_Stop_Depth, dive->Step_Size);

        dive->Starting_Depth = dive->Deco_Stop_Depth;
        dive->Next_Stop = dive->Deco_Stop_Depth - dive->Step_Size;
        dive->Deco_Stop_Depth = dive->Next_Stop;
        dive->Last_Run_Time = dive->Run_Time;
    }
}

void add_decomp_stop(dive_state *dive, double time, double depth, direction dir)
{
    if(dive->Real_Time_Decompression == TRUE) {
        /* printf("WOOOOOT %d", dive->decomp_stop_index); */
        dive->decomp_stops[dive->decomp_stop_index].time = time;
        dive->decomp_stops[dive->decomp_stop_index].depth = depth;
        dive->decomp_stops[dive->decomp_stop_index].ascent_or_const = dir;
        dive->decomp_stop_index++;
    }
}

void vpmb_critical_volume_decision_tree(dive_state *dive)
{
    int i;

    for(i=0; i < 16; i++) {
        dive->Helium_Pressure[i] = dive->He_Pressure_Start_of_Ascent[i];
        dive->Nitrogen_Pressure[i] = dive->N2_Pressure_Start_of_Ascent[i];
    }
    dive->Run_Time = dive->Run_Time_Start_of_Ascent;
    dive->Segment_Number = dive->Segment_Number_Start_of_Ascent;
    dive->Starting_Depth = dive->Depth_Change[0];
    dive->Mix_Number = dive->Mix_Change[0];
    dive->Rate = dive->Rate_Change[0];
    dive->Step_Size = dive->Step_Size_Change[0];
    dive->Deco_Stop_Depth = dive->First_Stop_Depth;
    dive->Last_Run_Time = 0.0;

    if(dive->Real_Time_Decompression == TRUE) {
        int size = (ceil(dive->First_Stop_Depth / dive->Step_Size) + 1) * 2;
        dive->decomp_stops = calloc(size, sizeof(decompression_stops));
        dive->decomp_stop_index = 0;
    }
    while(TRUE) {
        vpmb_gas_loadings_ascent_descent(dive, dive->Starting_Depth, dive->Deco_Stop_Depth, dive->Rate);
        vpmb_calc_max_actual_gradient(dive, dive->Deco_Stop_Depth);

        /* dive->output_object.add_decompression_profile_ascent(dive->Segment_Number, dive->Segment_Time, dive->Run_Time, dive->Mix_Number, dive->Deco_Stop_Depth, dive->Rate) */
        /* dive->decomp_stops[dive->decomp_stop_index].time = dive->Run_Time; */
        /* dive->decomp_stops[dive->decomp_stop_index].depth = dive->dive->Deco_Stop_Depth; */
        /* dive->decomp_stop_index++; */
        add_decomp_stop(dive, dive->Run_Time, dive->Deco_Stop_Depth, ASCENT);
        if (dive->Deco_Stop_Depth <= 0.0) {
            break;
        }

        if( dive->Number_of_Changes > 1) {
            for(i=1; i <  dive->Number_of_Changes; i++) {
                if(dive->Depth_Change[i] >= dive->Deco_Stop_Depth) {
                    dive->Mix_Number = dive->Mix_Change[i];
                    dive->Rate = dive->Rate_Change[i];
                    dive->Step_Size = dive->Step_Size_Change[i];
                }
            }
        }

        vpmb_boyles_law_compensation(dive, dive->First_Stop_Depth, dive->Deco_Stop_Depth, dive->Step_Size);
        vpmb_decompression_stop(dive, dive->Deco_Stop_Depth, dive->Step_Size);

        if ( dive->Last_Run_Time == 0.0) {
            dive->Stop_Time = round((dive->Segment_Time / dive->Minimum_Deco_Stop_Time) + 0.5) * dive->Minimum_Deco_Stop_Time;
        } else {
            dive->Stop_Time = dive->Run_Time - dive->Last_Run_Time;
        }

        if (trunc(dive->Minimum_Deco_Stop_Time) == dive->Minimum_Deco_Stop_Time) {
            /* dive->output_object.add_decompression_profile_constant(dive->Segment_Number, dive->Segment_Time, dive->Run_Time, dive->Mix_Number, int(dive->Deco_Stop_Depth), int(dive->Stop_Time)) */
            /* dive->decomp_stops[dive->decomp_stop_index].time = dive->Run_Time; */
            /* dive->decomp_stops[dive->decomp_stop_index].depth = dive->dive->Deco_Stop_Depth; */
            /* dive->decomp_stop_index++; */
            add_decomp_stop(dive, dive->Run_Time, dive->Deco_Stop_Depth, CONSTANT);

        } else {
            /* dive->output_object.add_decompression_profile_constant(dive->Segment_Number, dive->Segment_Time, dive->Run_Time, dive->Mix_Number, dive->Deco_Stop_Depth, dive->Stop_Time) */
            /* dive->decomp_stops[dive->decomp_stop_index].time = dive->Run_Time; */
            /* dive->decomp_stops[dive->decomp_stop_index].depth = dive->dive->Deco_Stop_Depth; */
            /* dive->decomp_stop_index++; */
            add_decomp_stop(dive, dive->Run_Time, dive->Deco_Stop_Depth, CONSTANT);

        }
        dive->Starting_Depth = dive->Deco_Stop_Depth;
        dive->Next_Stop = dive->Deco_Stop_Depth - dive->Step_Size;
        dive->Deco_Stop_Depth = dive->Next_Stop;
        dive->Last_Run_Time = dive->Run_Time;
    }
}

int critical_volume_loop(dive_state *dive)
{
    int i;
    while(TRUE) {

        vpmb_calc_ascent_ceiling(dive);
        if(dive->Ascent_Ceiling_Depth <= 0.0) {
            dive->Deco_Stop_Depth = 0.0;
        } else {
            double Rounding_Operation2 = (dive->Ascent_Ceiling_Depth / dive->Step_Size) + 0.5;
            dive->Deco_Stop_Depth = round(Rounding_Operation2) * dive->Step_Size;
        }
        if( dive->Deco_Stop_Depth > dive->Depth_Start_of_Deco_Zone) {
            return BADDECOSTOP;
        }


        vpmb_projected_ascent(dive, dive->Depth_Start_of_Deco_Zone, dive->Rate, dive->Step_Size);

        if (dive->Deco_Stop_Depth > dive->Depth_Start_of_Deco_Zone) {
            return BADDECOSTOP;
        }

        if(dive->Deco_Stop_Depth == 0.0) {
            for(i = 0; i < 16; i++) {
                dive->Helium_Pressure[i] = dive->He_Pressure_Start_of_Ascent[i];
                dive->Nitrogen_Pressure[i] = dive->N2_Pressure_Start_of_Ascent[i];
            }
            dive->Run_Time = dive-> Run_Time_Start_of_Ascent;
            dive->Segment_Number = dive->Segment_Number_Start_of_Ascent;
            dive->Starting_Depth = dive->Depth_Change[0];
            dive->Ending_Depth = 0.0;
            vpmb_gas_loadings_ascent_descent(dive, dive->Starting_Depth, dive->Ending_Depth, dive->Rate);

            /* self->output_object.add_decompression_profile_ascent(dive->Segment_Number, dive->Segment_Time, dive->Run_Time, dive->Mix_Number, dive->Deco_Stop_Depth, dive->Rate); */
            add_decomp_stop(dive, dive->Run_Time, dive->Deco_Stop_Depth, ASCENT);

            break;
        }

        dive->Starting_Depth = dive->Depth_Start_of_Deco_Zone;
        dive->First_Stop_Depth = dive->Deco_Stop_Depth;
        vpmb_deco_stop_loop_block_within_critical_volume_loop(dive);

        dive->Deco_Phase_Volume_Time = dive->Run_Time - dive->Run_Time_Start_of_Deco_Zone;

        vpmb_calc_surface_phase_volume_time(dive);

        for(i=0; i < 16; i++) {
            dive->Phase_Volume_Time[i] = dive->Deco_Phase_Volume_Time + dive->Surface_Phase_Volume_Time[i];
            dive->Critical_Volume_Comparison = fabs(dive->Phase_Volume_Time[i] - dive->Last_Phase_Volume_Time[i]);
            if(dive->Critical_Volume_Comparison <= 1.0) {
                dive->Schedule_Converged = TRUE;
            }
        }


        if ((dive->Schedule_Converged ==TRUE)|| (dive->Critical_Volume_Algorithm_Off== TRUE)) {
            vpmb_critical_volume_decision_tree(dive);
        }

        else {
            vpmb_critical_volume(dive, dive->Deco_Phase_Volume_Time);
            dive->Deco_Phase_Volume_Time = 0.0;
            dive->Run_Time = dive->Run_Time_Start_of_Deco_Zone;
            dive->Starting_Depth = dive->Depth_Start_of_Deco_Zone;
            dive->Mix_Number = dive->Mix_Change[0];
            dive->Rate = dive->Rate_Change[0];
            dive->Step_Size = dive->Step_Size_Change[0];

            for(i=0; i < 16; i++) {
                dive->Last_Phase_Volume_Time[i] = dive->Phase_Volume_Time[i];
                dive->Helium_Pressure[i] = dive->He_Pressure_Start_of_Deco_Zone[i];
                dive->Nitrogen_Pressure[i] = dive->N2_Pressure_Start_of_Deco_Zone[i];
            }
            continue;
        }
        break;
    }
    return ALLGOOD;
}
void vpmb_decompression_loop(dive_state *dive, single_dive *current_dive)
{
    int i,j;
    vpmb_nuclear_regeneration(dive, dive->Run_Time);

    vpmb_calc_initial_allowable_gradient(dive);

    for(i=0; i < 16; i++) {
        dive->He_Pressure_Start_of_Ascent[i] = dive->Helium_Pressure[i];
        dive->N2_Pressure_Start_of_Ascent[i] = dive->Nitrogen_Pressure[i];
    }

    dive->Run_Time_Start_of_Ascent = dive->Run_Time;
    dive->Segment_Number_Start_of_Ascent = dive->Segment_Number;

    for(i=0; i < current_dive->num_profile_codes; i++) {
        dive_profile *current_profile = &(current_dive->dive_profiles[i]);

        int Profile_Code = current_profile->profile_code;

        if(Profile_Code == 99) {
            dive->Number_of_Changes = current_profile->number_of_ascent_parameter_changes;

            // TODO remove hard coding of 16 and malloc this
            //dive->Depth_Change = [0.0 for i in range(dive->Number_of_Changes)]
            //        dive->Mix_Change = [0.0 for i in range(dive->Number_of_Changes)]
            //dive->Rate_Change = [0.0 for i in range(dive->Number_of_Changes)]
            //dive->Step_Size_Change = [0.0 for i in range(dive->Number_of_Changes)]

            for(j=0; j < current_profile->number_of_ascent_parameter_changes; j++) {
                ascent_summary *current_ascent = &(current_profile->ascents[j]);
                dive->Depth_Change[j] = current_ascent->starting_depth;
                dive->Mix_Change[j] = current_ascent->gasmix;
                dive->Rate_Change[j] = current_ascent->rate;
                dive->Step_Size_Change[j] = current_ascent->step_size;
            }
        }
    }

    dive->Starting_Depth = dive->Depth_Change[0];
    dive->Mix_Number = dive->Mix_Change[0];
    dive->Rate = dive->Rate_Change[0];
    dive->Step_Size = dive->Step_Size_Change[0];

    vpmb_calc_start_of_deco_zone(dive, dive->Starting_Depth, dive->Rate);

    if(dive->units_fsw == TRUE) {
        if (dive->Step_Size < 10.0) {
            double rounding_op = (dive->Depth_Start_of_Deco_Zone / dive->Step_Size) - 0.5;
            dive->Deepest_Possible_Stop_Depth = round(rounding_op) * dive->Step_Size;
        } else {
            double rounding_op = (dive->Depth_Start_of_Deco_Zone/10.0) - 0.5;
            dive->Deepest_Possible_Stop_Depth = round(rounding_op) * 10.0;
        }
    } else {
        if(dive->Step_Size < 3.0) {
            double rounding_op = (dive->Depth_Start_of_Deco_Zone / dive->Step_Size) - 0.5;
            dive->Deepest_Possible_Stop_Depth = round(rounding_op) * dive->Step_Size;
        } else {
            double rounding_op = (dive->Depth_Start_of_Deco_Zone / 3.0)  - 0.5;
            dive->Deepest_Possible_Stop_Depth = round(rounding_op) * 3.0;
        }
    }

    vpmb_gas_loadings_ascent_descent(dive, dive->Starting_Depth, dive->Depth_Start_of_Deco_Zone, dive->Rate);
    dive->Run_Time_Start_of_Deco_Zone = dive->Run_Time;
    dive->Deco_Phase_Volume_Time = 0.0;
    dive->Last_Run_Time = 0.0;
    dive->Schedule_Converged = FALSE;

    for(i=0; i < 16; i++) {
        dive->Last_Phase_Volume_Time[i] = 0.0;
        dive->He_Pressure_Start_of_Deco_Zone[i] = dive->Helium_Pressure[i];
        dive->N2_Pressure_Start_of_Deco_Zone[i] = dive->Nitrogen_Pressure[i];
        dive->Max_Actual_Gradient[i] = 0.0;
    }

    critical_volume_loop(dive);
}



int vpmb_set_gas_mixes(dive_state *dive, single_dive *current_dive)
{
    int i;
    int num_gas_mixes = current_dive->num_gas_mixes;
    double Fraction_Oxygen;
    gasmix_summary *current_gasmix;

    if(dive->Fraction_Helium != NULL) {
        free(dive->Fraction_Helium);
    }
    if(dive->Fraction_Nitrogen != NULL) {
        free(dive->Fraction_Nitrogen);
    }

    dive->Fraction_Helium = calloc(num_gas_mixes, sizeof(double));
    dive->Fraction_Nitrogen = calloc(num_gas_mixes, sizeof(double));

    for(i=0; i < num_gas_mixes; i++) {

        current_gasmix = &(current_dive->gasmixes[i]);
        Fraction_Oxygen = current_gasmix->fraction_O2;
        dive->Fraction_Nitrogen[i] = current_gasmix->fraction_N2;
        dive->Fraction_Helium[i] = current_gasmix->fraction_He;

        /* make sure the fractions add up */
        if((Fraction_Oxygen + dive->Fraction_Nitrogen[i] + dive->Fraction_Helium[i]) != 1.0) {
            return BADJSON;
        }

        //for i in range(num_gas_mixes):
        //  dive->output_object.add_gasmix(Fraction_Oxygen[i], self.Fraction_Nitrogen[i], self.Fraction_Helium[i])
    }
    return GOODJSON;
}

void vpmb_repetitive_dive_loop(dive_state *dive, json_input *input)
{
    int i, j;
    single_dive *current_dive;
    dive->Real_Time_Decompression=FALSE;
    for(i=0; i < input->number_of_dives; i++) {
        int Repetitive_Dive_Flag;

        current_dive = &(input->dives[i]);
        //self.output_object.new_dive(dive["desc"])

        vpmb_set_gas_mixes(dive, current_dive);
        vpmb_profile_code_loop(dive, current_dive);
        vpmb_decompression_loop(dive, current_dive);

        Repetitive_Dive_Flag = current_dive->repetitive_code;

        if(Repetitive_Dive_Flag == 0) {
            continue;
        }

        else if (Repetitive_Dive_Flag == 1) {
            dive->Surface_Interval_Time = current_dive->surface_interval_time_minutes;

            vpmb_gas_loadings_surface_interval(dive, dive->Surface_Interval_Time);
            vpmb_vpm_repetitive_algorithm(dive, dive->Surface_Interval_Time);

            for(j = 0; j < 16; j++) {
                dive->Max_Crushing_Pressure_He[j] = 0.0;
                dive->Max_Crushing_Pressure_N2[j] = 0.0;
                dive->Max_Actual_Gradient[j] = 0.0;
            }
            dive->Run_Time = 0.0;
            dive->Segment_Number = 0;

            /* may not be needed anymore */
            continue;
        }
    }
}


direction vpmb_current_direction(dive_state *dive, double increment_time)
{
    double high = dive->Last_Direction_Depth + (meters_per_minute_change * increment_time), low = dive->Last_Direction_Depth - (meters_per_minute_change * increment_time);

    /* double rate = (dive->Depth - dive->Last_Direction_Depth ) / increment_time; */
    if((dive->Depth <= high) && (dive->Depth >= low)) {
        dive->Last_Direction_Depth = dive->Depth;
        return CONSTANT;
    } else if(dive->Depth > dive->Last_Direction_Depth) {
        dive->Last_Direction_Depth = dive->Depth;
        return DESCENT;
    } else if(dive->Depth < dive->Last_Direction_Depth) {
        dive->Last_Direction_Depth = dive->Depth;
        return ASCENT;
    }

    return ERROR;
}

/* int msw_driver(dive_state *dive){ */

/*         if(dive->Run_Time >= 30 && dive->Depth > 0){ */
/*                 return 99; */
/*         } */
/*         else if(dive->Depth < 80.0) */
/*                 return 1; */
/*         else if (dive->Run_Time < 30){ */
/*                 return 2; */
/*         } */
/*         return 100; */
/* } */


int vpmb_finished_constant_depth_profile(dive_state *dive, dive_profile *current_profile)
{

    if(dive->Run_Time >= current_profile->run_time_at_end_of_segment) {
        return DONE_PROFILE;
    }
    return NOT_DONE_PROFILE;
}

double vpmb_find_start_of_decompression_zone(dive_state *dive, dive_profile *current_profile)
{
    int i;//, j;

    ascent_summary *current_ascent;
    dive_state *dive_cpy = malloc(sizeof(dive_state));
    memcpy (dive_cpy, dive,sizeof(dive_state));

    vpmb_nuclear_regeneration(dive_cpy, dive_cpy->Run_Time);
    vpmb_calc_initial_allowable_gradient(dive_cpy);

    for(i=0; i < 16; i++) {
        dive_cpy->He_Pressure_Start_of_Ascent[i] = dive_cpy->Helium_Pressure[i];
        dive_cpy->N2_Pressure_Start_of_Ascent[i] = dive_cpy->Nitrogen_Pressure[i];
    }

    dive_cpy->Run_Time_Start_of_Ascent = dive_cpy->Run_Time;
    dive_cpy->Segment_Number_Start_of_Ascent = dive_cpy->Segment_Number;

    /* for(i=0; i < current_dive->num_profile_codes; i++) { */
    /*         dive_cpy->Number_of_Changes = current_profile->number_of_ascent_parameter_changes; */

    /*         for(j=0; j < current_profile->number_of_ascent_parameter_changes; j++){ */
    /*                 ascent_summary *current_ascent = &(current_profile->ascents[j]); */
    /*                 dive_cpy->Depth_Change[j] = current_ascent->starting_depth; */
    /*                 dive_cpy->Mix_Change[j] = current_ascent->gasmix; */
    /*                 dive_cpy->Rate_Change[j] = current_ascent->rate; */
    /*                 dive_cpy->Step_Size_Change[j] = current_ascent->step_size; */
    /*         } */
    /* } */

    //dive_profile *current_profile = &(current_dive->dive_profiles[2]);
    current_ascent = &(current_profile->ascents[0]);

    dive_cpy->Starting_Depth = dive_cpy->Depth; //dive_cpy->Depth_Change[0];
    dive_cpy->Mix_Number = current_ascent->gasmix; //dive_cpy->Mix_Change[0];
    dive_cpy->Rate = current_ascent->rate;//dive_cpy->Rate_Change[0];
    dive_cpy->Step_Size = current_ascent->step_size;//dive_cpy->Step_Size_Change[0];

    vpmb_calc_start_of_deco_zone(dive_cpy, dive_cpy->Depth, dive_cpy->Rate);

    if(dive_cpy->units_fsw == TRUE) {
        if (dive_cpy->Step_Size < 10.0) {
            double rounding_op = (dive_cpy->Depth_Start_of_Deco_Zone / dive_cpy->Step_Size) - 0.5;
            dive_cpy->Deepest_Possible_Stop_Depth = round(rounding_op) * dive_cpy->Step_Size;
        } else {
            double rounding_op = (dive_cpy->Depth_Start_of_Deco_Zone/10.0) - 0.5;
            dive_cpy->Deepest_Possible_Stop_Depth = round(rounding_op) * 10.0;
        }
    } else {
        if(dive_cpy->Step_Size < 3.0) {
            double rounding_op = (dive_cpy->Depth_Start_of_Deco_Zone / dive_cpy->Step_Size) - 0.5;
            dive_cpy->Deepest_Possible_Stop_Depth = round(rounding_op) * dive_cpy->Step_Size;
        } else {
            double rounding_op = (dive_cpy->Depth_Start_of_Deco_Zone / 3.0)  - 0.5;
            dive_cpy->Deepest_Possible_Stop_Depth = round(rounding_op) * 3.0;
        }
    }

    double deepest_depth = dive_cpy->Deepest_Possible_Stop_Depth;
    free(dive_cpy);
    if(deepest_depth < 0) {
        deepest_depth = 0;
    }
    return deepest_depth;
}




void vpmb_free_dive_state(dive_state *dive)
{
    free(dive->Fraction_Helium);
    free(dive->Fraction_Nitrogen);

    if(dive->decomp_stops!=NULL) {
        free(dive->decomp_stops);
    }
    free(dive);
}

void vpmb_calculate_decompression_stops(dive_state *dive, dive_profile *current_profile)
{
    int i, j;
    dive_state *dive_cpy = malloc(sizeof(dive_state));
    memcpy(dive_cpy, dive, sizeof(dive_state));

    vpmb_nuclear_regeneration(dive_cpy, dive_cpy->Run_Time);

    vpmb_calc_initial_allowable_gradient(dive_cpy);

    for(i=0; i < 16; i++) {
        dive_cpy->He_Pressure_Start_of_Ascent[i] = dive_cpy->Helium_Pressure[i];
        dive_cpy->N2_Pressure_Start_of_Ascent[i] = dive_cpy->Nitrogen_Pressure[i];
    }

    dive_cpy->Run_Time_Start_of_Ascent = dive_cpy->Run_Time;
    dive_cpy->Segment_Number_Start_of_Ascent = dive_cpy->Segment_Number;

    //                for(i=0; i < current_dive_cpy->num_profile_codes; i++) {
    //dive_cpy_profile *current_profile = &(current_dive_cpy->dive_cpy_profiles[i]);

    //int Profile_Code = current_profile->profile_code;

    //if(Profile_Code == 99){
    dive_cpy->Number_of_Changes = current_profile->number_of_ascent_parameter_changes;

    // TODO remove hard coding of 16 and malloc this
    //dive->Depth_Change = [0.0 for i in range(dive->Number_of_Changes)]
    //        dive->Mix_Change = [0.0 for i in range(dive->Number_of_Changes)]
    //dive->Rate_Change = [0.0 for i in range(dive->Number_of_Changes)]
    //dive->Step_Size_Change = [0.0 for i in range(dive->Number_of_Changes)]

    for(j=0; j < current_profile->number_of_ascent_parameter_changes; j++) {
        ascent_summary *current_ascent = &(current_profile->ascents[j]);
        dive_cpy->Depth_Change[j] = current_ascent->starting_depth;
        dive_cpy->Mix_Change[j] = current_ascent->gasmix;
        dive_cpy->Rate_Change[j] = current_ascent->rate;
        dive_cpy->Step_Size_Change[j] = current_ascent->step_size;
    }
    //      }
    //}

    dive_cpy->Starting_Depth = dive_cpy->Depth_Change[0];
    dive_cpy->Mix_Number = dive_cpy->Mix_Change[0];
    dive_cpy->Rate = dive_cpy->Rate_Change[0];
    dive_cpy->Step_Size = dive_cpy->Step_Size_Change[0];

    vpmb_calc_start_of_deco_zone(dive_cpy, dive_cpy->Starting_Depth, dive_cpy->Rate);

    if(dive_cpy->units_fsw == TRUE) {
        if (dive_cpy->Step_Size < 10.0) {
            double rounding_op = (dive_cpy->Depth_Start_of_Deco_Zone / dive_cpy->Step_Size) - 0.5;
            dive_cpy->Deepest_Possible_Stop_Depth = round(rounding_op) * dive_cpy->Step_Size;
        } else {
            double rounding_op = (dive_cpy->Depth_Start_of_Deco_Zone/10.0) - 0.5;
            dive_cpy->Deepest_Possible_Stop_Depth = round(rounding_op) * 10.0;
        }
    } else {
        if(dive_cpy->Step_Size < 3.0) {
            double rounding_op = (dive_cpy->Depth_Start_of_Deco_Zone / dive_cpy->Step_Size) - 0.5;
            dive_cpy->Deepest_Possible_Stop_Depth = round(rounding_op) * dive_cpy->Step_Size;
        } else {
            double rounding_op = (dive_cpy->Depth_Start_of_Deco_Zone / 3.0)  - 0.5;
            dive_cpy->Deepest_Possible_Stop_Depth = round(rounding_op) * 3.0;
        }
    }

    vpmb_gas_loadings_ascent_descent(dive_cpy, dive_cpy->Starting_Depth, dive_cpy->Depth_Start_of_Deco_Zone, dive_cpy->Rate);
    dive_cpy->Run_Time_Start_of_Deco_Zone = dive_cpy->Run_Time;
    dive_cpy->Deco_Phase_Volume_Time = 0.0;
    dive_cpy->Last_Run_Time = 0.0;
    dive_cpy->Schedule_Converged = FALSE;

    for(i=0; i < 16; i++) {
        dive_cpy->Last_Phase_Volume_Time[i] = 0.0;
        dive_cpy->He_Pressure_Start_of_Deco_Zone[i] = dive_cpy->Helium_Pressure[i];
        dive_cpy->N2_Pressure_Start_of_Deco_Zone[i] = dive_cpy->Nitrogen_Pressure[i];
        dive_cpy->Max_Actual_Gradient[i] = 0.0;
    }

    critical_volume_loop(dive_cpy);

    dive->decomp_stops = (decompression_stops*) calloc(dive_cpy->decomp_stop_index, sizeof(decompression_stops));
    //int i;
    for(i = 0; i < dive_cpy->decomp_stop_index; i++) {
        dive->decomp_stops[i] = dive_cpy->decomp_stops[i];
    }
    dive->decomp_stop_index = 0;
    free(dive_cpy->decomp_stops);
    free(dive_cpy);
}


void vpmb_critical_volume_decision_tree_to_depth(dive_state *dive, double stop_depth)
{
    int i;

    for(i=0; i < 16; i++) {
        dive->Helium_Pressure[i] = dive->He_Pressure_Start_of_Ascent[i];
        dive->Nitrogen_Pressure[i] = dive->N2_Pressure_Start_of_Ascent[i];
    }
    dive->Run_Time = dive->Run_Time_Start_of_Ascent;
    dive->Segment_Number = dive->Segment_Number_Start_of_Ascent;
    dive->Starting_Depth = dive->Depth_Change[0];
    dive->Mix_Number = dive->Mix_Change[0];
    dive->Rate = dive->Rate_Change[0];
    dive->Step_Size = dive->Step_Size_Change[0];
    dive->Deco_Stop_Depth = dive->First_Stop_Depth;
    dive->Last_Run_Time = 0.0;

    while(TRUE) {
        vpmb_gas_loadings_ascent_descent(dive, dive->Starting_Depth, dive->Deco_Stop_Depth, dive->Rate);
        vpmb_calc_max_actual_gradient(dive, dive->Deco_Stop_Depth);

        //dive->output_object.add_decompression_profile_ascent(dive->Segment_Number, dive->Segment_Time, dive->Run_Time, dive->Mix_Number, dive->Deco_Stop_Depth, dive->Rate)
        /* dive->decomp_stops[dive->decomp_stop_index].time = dive->Run_Time; */
        /* dive->decomp_stops[dive->decomp_stop_index].depth = dive->dive->Deco_Stop_Depth; */
        /* dive->decomp_stop_index++; */
        //add_decomp_stop(dive, dive->Run_Time, dive->Deco_Stop_Depth, ASCENT);
        if (dive->Deco_Stop_Depth <= stop_depth) {
            break;
        }

        if( dive->Number_of_Changes > 1) {
            for(i=1; i <  dive->Number_of_Changes; i++) {
                if(dive->Depth_Change[i] >= dive->Deco_Stop_Depth) {
                    dive->Mix_Number = dive->Mix_Change[i];
                    dive->Rate = dive->Rate_Change[i];
                    dive->Step_Size = dive->Step_Size_Change[i];
                }
            }
        }

        vpmb_boyles_law_compensation(dive, dive->First_Stop_Depth, dive->Deco_Stop_Depth, dive->Step_Size);
        vpmb_decompression_stop(dive, dive->Deco_Stop_Depth, dive->Step_Size);

        if ( dive->Last_Run_Time == 0.0) {
            dive->Stop_Time = round((dive->Segment_Time / dive->Minimum_Deco_Stop_Time) + 0.5) * dive->Minimum_Deco_Stop_Time;
        } else {
            dive->Stop_Time = dive->Run_Time - dive->Last_Run_Time;
        }

        if (trunc(dive->Minimum_Deco_Stop_Time) == dive->Minimum_Deco_Stop_Time) {
            //dive->output_object.add_decompression_profile_constant(dive->Segment_Number, dive->Segment_Time, dive->Run_Time, dive->Mix_Number, int(dive->Deco_Stop_Depth), int(dive->Stop_Time))
            /* dive->decomp_stops[dive->decomp_stop_index].time = dive->Run_Time; */
            /* dive->decomp_stops[dive->decomp_stop_index].depth = dive->dive->Deco_Stop_Depth; */
            /* dive->decomp_stop_index++; */
            //add_decomp_stop(dive, dive->Run_Time, dive->Deco_Stop_Depth, CONSTANT);

        } else {
            //dive->output_object.add_decompression_profile_constant(dive->Segment_Number, dive->Segment_Time, dive->Run_Time, dive->Mix_Number, dive->Deco_Stop_Depth, dive->Stop_Time)
            /* dive->decomp_stops[dive->decomp_stop_index].time = dive->Run_Time; */
            /* dive->decomp_stops[dive->decomp_stop_index].depth = dive->dive->Deco_Stop_Depth; */
            /* dive->decomp_stop_index++; */
            //add_decomp_stop(dive, dive->Run_Time, dive->Deco_Stop_Depth, CONSTANT);

        }
        dive->Starting_Depth = dive->Deco_Stop_Depth;
        dive->Next_Stop = dive->Deco_Stop_Depth - dive->Step_Size;
        dive->Deco_Stop_Depth = dive->Next_Stop;
        dive->Last_Run_Time = dive->Run_Time;
    }
}


int vpmb_critical_volume_loop_init(dive_state *dive)
{
    int i;
    while(TRUE) {

        vpmb_calc_ascent_ceiling(dive);
        if(dive->Ascent_Ceiling_Depth <= 0.0) {
            dive->Deco_Stop_Depth = 0.0;
        } else {
            double Rounding_Operation2 = (dive->Ascent_Ceiling_Depth / dive->Step_Size) + 0.5;
            dive->Deco_Stop_Depth = round(Rounding_Operation2) * dive->Step_Size;
        }
        if( dive->Deco_Stop_Depth > dive->Depth_Start_of_Deco_Zone) {
            return BADDECOSTOP;
        }


        vpmb_projected_ascent(dive, dive->Depth_Start_of_Deco_Zone, dive->Rate, dive->Step_Size);

        if (dive->Deco_Stop_Depth > dive->Depth_Start_of_Deco_Zone) {
            return BADDECOSTOP;
        }

        if(dive->Deco_Stop_Depth == 0.0) {
            for(i = 0; i < 16; i++) {
                dive->Helium_Pressure[i] = dive->He_Pressure_Start_of_Ascent[i];
                dive->Nitrogen_Pressure[i] = dive->N2_Pressure_Start_of_Ascent[i];
            }
            dive->Run_Time = dive-> Run_Time_Start_of_Ascent;
            dive->Segment_Number = dive->Segment_Number_Start_of_Ascent;
            dive->Starting_Depth = dive->Depth_Change[0];
            dive->Ending_Depth = 0.0;
            vpmb_gas_loadings_ascent_descent(dive, dive->Starting_Depth, dive->Ending_Depth, dive->Rate);

            //self->output_object.add_decompression_profile_ascent(dive->Segment_Number, dive->Segment_Time, dive->Run_Time, dive->Mix_Number, dive->Deco_Stop_Depth, dive->Rate);
            add_decomp_stop(dive, dive->Run_Time, dive->Deco_Stop_Depth, ASCENT);

            break;
        }

        dive->Starting_Depth = dive->Depth_Start_of_Deco_Zone;
        dive->First_Stop_Depth = dive->Deco_Stop_Depth;
        vpmb_deco_stop_loop_block_within_critical_volume_loop(dive);

        dive->Deco_Phase_Volume_Time = dive->Run_Time - dive->Run_Time_Start_of_Deco_Zone;

        vpmb_calc_surface_phase_volume_time(dive);

        for(i=0; i < 16; i++) {
            dive->Phase_Volume_Time[i] = dive->Deco_Phase_Volume_Time + dive->Surface_Phase_Volume_Time[i];
            dive->Critical_Volume_Comparison = fabs(dive->Phase_Volume_Time[i] - dive->Last_Phase_Volume_Time[i]);
            if(dive->Critical_Volume_Comparison <= 1.0) {
                dive->Schedule_Converged = TRUE;
            }
        }


        if ((dive->Schedule_Converged ==TRUE)|| (dive->Critical_Volume_Algorithm_Off== TRUE))
            //vpmb_critical_volume_decision_tree(dive);
        {
            return ALLGOOD;
        }

        else {
            vpmb_critical_volume(dive, dive->Deco_Phase_Volume_Time);
            dive->Deco_Phase_Volume_Time = 0.0;
            dive->Run_Time = dive->Run_Time_Start_of_Deco_Zone;
            dive->Starting_Depth = dive->Depth_Start_of_Deco_Zone;
            dive->Mix_Number = dive->Mix_Change[0];
            dive->Rate = dive->Rate_Change[0];
            dive->Step_Size = dive->Step_Size_Change[0];

            for(i=0; i < 16; i++) {
                dive->Last_Phase_Volume_Time[i] = dive->Phase_Volume_Time[i];
                dive->Helium_Pressure[i] = dive->He_Pressure_Start_of_Deco_Zone[i];
                dive->Nitrogen_Pressure[i] = dive->N2_Pressure_Start_of_Deco_Zone[i];
            }
            continue;
        }
        break;
    }
    return ALLGOOD;
}

void vpmb_decompress_init(dive_state *dive, dive_profile *current_profile)
{
    int i,j;
    vpmb_nuclear_regeneration(dive, dive->Run_Time);

    vpmb_calc_initial_allowable_gradient(dive);

    for(i=0; i < 16; i++) {
        dive->He_Pressure_Start_of_Ascent[i] = dive->Helium_Pressure[i];
        dive->N2_Pressure_Start_of_Ascent[i] = dive->Nitrogen_Pressure[i];
    }

    dive->Run_Time_Start_of_Ascent = dive->Run_Time;
    dive->Segment_Number_Start_of_Ascent = dive->Segment_Number;

    dive->Number_of_Changes = current_profile->number_of_ascent_parameter_changes;

    for(j=0; j < current_profile->number_of_ascent_parameter_changes; j++) {
        ascent_summary *current_ascent = &(current_profile->ascents[j]);
        dive->Depth_Change[j] = current_ascent->starting_depth;
        dive->Mix_Change[j] = current_ascent->gasmix;
        dive->Rate_Change[j] = current_ascent->rate;
        dive->Step_Size_Change[j] = current_ascent->step_size;
    }

    dive->Starting_Depth = dive->Depth_Change[0];
    dive->Mix_Number = dive->Mix_Change[0];
    dive->Rate = dive->Rate_Change[0];
    dive->Step_Size = dive->Step_Size_Change[0];

    vpmb_calc_start_of_deco_zone(dive, dive->Starting_Depth, dive->Rate);

    if(dive->units_fsw == TRUE) {
        if (dive->Step_Size < 10.0) {
            double rounding_op = (dive->Depth_Start_of_Deco_Zone / dive->Step_Size) - 0.5;
            dive->Deepest_Possible_Stop_Depth = round(rounding_op) * dive->Step_Size;
        } else {
            double rounding_op = (dive->Depth_Start_of_Deco_Zone/10.0) - 0.5;
            dive->Deepest_Possible_Stop_Depth = round(rounding_op) * 10.0;
        }
    } else {
        if(dive->Step_Size < 3.0) {
            double rounding_op = (dive->Depth_Start_of_Deco_Zone / dive->Step_Size) - 0.5;
            dive->Deepest_Possible_Stop_Depth = round(rounding_op) * dive->Step_Size;
        } else {
            double rounding_op = (dive->Depth_Start_of_Deco_Zone / 3.0)  - 0.5;
            dive->Deepest_Possible_Stop_Depth = round(rounding_op) * 3.0;
        }
    }

    vpmb_gas_loadings_ascent_descent(dive, dive->Starting_Depth, dive->Depth_Start_of_Deco_Zone, dive->Rate);
    dive->Run_Time_Start_of_Deco_Zone = dive->Run_Time;
    dive->Deco_Phase_Volume_Time = 0.0;
    dive->Last_Run_Time = 0.0;
    dive->Schedule_Converged = FALSE;

    for(i=0; i < 16; i++) {
        dive->Last_Phase_Volume_Time[i] = 0.0;
        dive->He_Pressure_Start_of_Deco_Zone[i] = dive->Helium_Pressure[i];
        dive->N2_Pressure_Start_of_Deco_Zone[i] = dive->Nitrogen_Pressure[i];
        dive->Max_Actual_Gradient[i] = 0.0;
    }
    vpmb_critical_volume_loop_init(dive);
}

