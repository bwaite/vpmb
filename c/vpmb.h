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
 */

#ifndef VPMB_H
#define VPMB_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include "strlcpy.h"
#include "cJSON.h"

//structs

typedef short BOOL;
#ifndef _NDS_BUILD_
extern const short TRUE;
extern const short FALSE;
#endif

/* errors */
extern const int ALLGOOD;
extern const int BADJSON;
extern const int GOODJSON;
extern const int INVALIDDATA;
extern const int VALIDDATA;

extern const int ROOTERROR;
extern const int ROOTFOUND;
extern const int BADALTITUDE;

extern const int OFFGASSINGERROR;
extern const int BADDECOSTOP;
extern const int NOFILE;

/* Constants */
extern const double ATM;
extern const double fraction_inert_gas;
extern const double Helium_Half_Time[];
extern const double Nitrogen_Half_Time[];
extern const int Buhlmann_Compartments;

extern const int DONE_PROFILE;
extern const int NOT_DONE_PROFILE;
extern const int DEPTH_CHANGED;

typedef enum  {ASCENT, DESCENT, CONSTANT,  ERROR} direction;

/* Structs */

/**
 * \struct decompression_stops
 * \brief Holds the time, depth and whether the diver is ascending or at constant depth
 */
typedef struct {
        double time; /**< Run time that the decompression stop should take place at */
        double depth; /**< Depth that the decompression stop should take place at */
        direction ascent_or_const; /**< Either ASCENT or CONSTANT, based on whether the diver should be ascending or at constant depth */
} decompression_stops;

/**
 * \struct dive_state
 * \brief Contains all the state needed to calculate the VPM-B algorithm
 */
typedef struct {
        char Units[5];
        char Units_Word1[10];
        char Units_Word2[10];

        double Water_Vapor_Pressure;
        double Surface_Tension_Gamma;
        double Skin_Compression_GammaC;
        double Crit_Volume_Parameter_Lambda;
        double Minimum_Deco_Stop_Time;
        double Regeneration_Time_Constant;
        double Constant_Pressure_Other_Gases;
        double Gradient_Onset_of_Imperm_Atm;

        int Number_of_Changes;
        int Segment_Number_Start_of_Ascent;
        int Repetitive_Dive_Flag;
        BOOL Schedule_Converged;
        BOOL Critical_Volume_Algorithm_Off;
        BOOL Altitude_Dive_Algorithm_Off;
        double Ascent_Ceiling_Depth;
        double Deco_Stop_Depth;
        double Step_Size;
        double Depth;
        double Last_Depth;
        double Ending_Depth;
        double Starting_Depth;
        double Rate;
        double Run_Time_End_of_Segment;
        double Last_Run_Time;
        double Stop_Time;
        double Depth_Start_of_Deco_Zone;
        double Deepest_Possible_Stop_Depth;
        double First_Stop_Depth;
        double Critical_Volume_Comparison;
        double Next_Stop;
        double Run_Time_Start_of_Deco_Zone;
        double Critical_Radius_N2_Microns;
        double Critical_Radius_He_Microns;
        double Run_Time_Start_of_Ascent;
        double Altitude_of_Dive;
        double Deco_Phase_Volume_Time;
        double Surface_Interval_Time;
        double Regenerated_Radius_He[16];
        double Regenerated_Radius_N2[16];

        int Mix_Change[16];
        int Depth_Change[16];
        int Rate_Change[16];
        int Step_Size_Change[16];
        double He_Pressure_Start_of_Ascent[16];
        double N2_Pressure_Start_of_Ascent[16];
        double He_Pressure_Start_of_Deco_Zone[16];
        double N2_Pressure_Start_of_Deco_Zone[16];
        double Phase_Volume_Time[16];
        double Last_Phase_Volume_Time[16];
        double Allowable_Gradient_He[16];
        double Allowable_Gradient_N2[16];
        double Adjusted_Crushing_Pressure_He[16];
        double Adjusted_Crushing_Pressure_N2[16];
        double Initial_Allowable_Gradient_N2[16];
        double Initial_Allowable_Gradient_He[16];
        double Deco_Gradient_He[16];
        double Deco_Gradient_N2[16];
        int Segment_Number;
        double Run_Time;
        double Segment_Time;
        double Ending_Ambient_Pressure;
        int Mix_Number;
        double Barometric_Pressure;
        BOOL units_fsw;
        double Units_Factor;
        double Helium_Time_Constant[16];
        double Nitrogen_Time_Constant[16];
        double Helium_Pressure[16];
        double Nitrogen_Pressure[16];
        double Initial_Helium_Pressure[16];
        double Initial_Nitrogen_Pressure[16];
        double *Fraction_Helium;
        double *Fraction_Nitrogen;
        double Initial_Critical_Radius_He[16];
        double Initial_Critical_Radius_N2[16];
        double Adjusted_Critical_Radius_He[16];
        double Adjusted_Critical_Radius_N2[16];
        double Max_Crushing_Pressure_He[16];
        double Max_Crushing_Pressure_N2[16];
        double Surface_Phase_Volume_Time[16];
        double Max_Actual_Gradient[16];
        double Amb_Pressure_Onset_of_Imperm[16];
        double Gas_Tension_Onset_of_Imperm[16];
        BOOL Diver_Acclimatized;

        double Last_Direction_Depth;
        double Last_Direction_Time;

        double Start_of_Decompression_Zone;
        BOOL Decompressing;

        decompression_stops *decomp_stops;
        int decomp_stop_index;
        BOOL Real_Time_Decompression;
        double Wait_Time;
} dive_state;

typedef struct{
        double starting_depth;
        int gasmix;
        double rate;
        double step_size;
        double setpoint;
} ascent_summary;

typedef struct{
        int profile_code;
        double starting_depth;
        double ending_depth;
        double rate;
        int gasmix;
        double setpoint;
        double depth;
        double run_time_at_end_of_segment;
        int number_of_ascent_parameter_changes;
        ascent_summary *ascents;
} dive_profile;

typedef struct {
        double fraction_O2;
        double fraction_He;
        double fraction_N2;
} gasmix_summary;

typedef struct {
        char desc[20];
        int num_gas_mixes;
        int repetitive_code;
        gasmix_summary *gasmixes;
        int num_profile_codes;
        dive_profile *dive_profiles;
        double surface_interval_time_minutes;
} single_dive;

typedef struct {
        //input
        single_dive *dives;
        int number_of_dives;

        //altitude
        double Altitude_of_Dive;
        char Diver_Acclimatized_at_Altitude[4];
        double Starting_Acclimatized_Altitude;
        double Ascent_to_Altitude_Hours;
        double Hours_at_Altitude_Before_Dive;

        //settings
        char Units[4]; //"fsw" or "msw"
        BOOL SetPoint_Is_Bar;
        char Altitude_Dive_Algorithm[4]; //"OFF" or "ON"
        double Minimum_Deco_Stop_Time;
        double Critical_Radius_N2_Microns;
        double Critical_Radius_He_Microns;
        char Critical_Volume_Algorithm[4];
        double Crit_Volume_Parameter_Lambda;
        double Gradient_Onset_of_Imperm_Atm;
        double Surface_Tension_Gamma;
        double Skin_Compression_GammaC;
        double Regeneration_Time_Constant;
        double Pressure_Other_Gases_mmHg;

} json_input;

/* Function Prototypes */
void failure(void);
double get_setpoint_or_null(const cJSON *setpoint);
int input_add_dive_profiles(dive_profile *in, cJSON *profile);
int input_add_dive(single_dive *in, cJSON *dive);
void free_dives(json_input *in);
int load_from_json(json_input *in, const char *filename);
void lowercase_string(char *str);
double SCHREINER_EQUATION(double Initial_Inspired_Gas_Pressure, double Rate_Change_Insp_Gas_Pressure, double Interval_Time, double Gas_Time_Constant, double Initial_Gas_Pressure);
double HALDANE_EQUATION(double Initial_Gas_Pressure, double Inspired_Gas_Pressure, double Gas_Time_Constant, double Interval_Time);
int RADIUS_ROOT_FINDER(double A, double B, double C, double Low_Bound, double High_Bound, double *result);
double CALC_BAROMETRIC_PRESSURE(double Altitude, BOOL units_fsw);
double min(double i, double j);
double max(double i, double j);
double CALC_DECO_CEILING(dive_state *dive);
int VPM_ALTITUDE_DIVE_ALGORITHM(json_input *input, dive_state *dive);
void GAS_LOADINGS_CONSTANT_DEPTH(dive_state *dive, double Depth, double Run_Time_End_of_Segment);
void NUCLEAR_REGENERATION(dive_state *dive, double Dive_Time);
double _crushing_pressure_helper(dive_state *dive, double Radius_Onset_of_Imperm_Molecule, double Ending_Ambient_Pressure_Pa, double Amb_Press_Onset_of_Imperm_Pa, double Gas_Tension_Onset_of_Imperm_Pa, double Gradient_Onset_of_Imperm_Pa);
int ONSET_OF_IMPERMEABILITY(dive_state *dive, double Starting_Ambient_Pressure, double Ending_Ambient_Pressure, double Rate, int i);
void CALC_CRUSHING_PRESSURE(dive_state *dive, double Starting_Depth, double Ending_Depth, double Rate);
void GAS_LOADINGS_ASCENT_DESCENT(dive_state *dive, double Starting_Depth, double Ending_Depth, double Rate);
int DECOMPRESSION_STOP(dive_state *dive, double Deco_Stop_Depth, double Step_Size);
double _calculate_deco_gradient(dive_state *dive, double Allowable_Gradient_Molecule, double Amb_Press_First_Stop_Pascals, double Amb_Press_Next_Stop_Pascals);
void BOYLES_LAW_COMPENSATION(dive_state *dive, double First_Stop_Depth, double Deco_Stop_Depth, double Step_Size);
void CALC_MAX_ACTUAL_GRADIENT(dive_state *dive, double Deco_Stop_Depth);
void PROJECTED_ASCENT(dive_state *dive, double Starting_Depth, double Rate, double Step_Size);
void CALC_ASCENT_CEILING(dive_state *dive);
void CALC_SURFACE_PHASE_VOLUME_TIME(dive_state *dive);
void CRITICAL_VOLUME(dive_state *dive, double Deco_Phase_Volume_Time);
int CALC_START_OF_DECO_ZONE(dive_state *dive, double Starting_Depth, double Rate);
void CALC_INITIAL_ALLOWABLE_GRADIENT(dive_state *dive);
void GAS_LOADINGS_SURFACE_INTERVAL(dive_state *dive, double Surface_Interval_Time);
double _new_critical_radius(dive_state *dive, double Max_Actual_Gradient_Pascals, double Adj_Crush_Pressure_Pascals);
void VPM_REPETITIVE_ALGORITHM(dive_state *dive, double Surface_Interval_Time);
int validate_data(json_input *input, dive_state *dive);
int initialize_data(json_input *input, dive_state *dive);
void profile_code_loop(dive_state *dive, single_dive *current_dive);
void deco_stop_loop_block_within_critical_volume_loop(dive_state *dive);
void critical_volume_decision_tree(dive_state *dive);
int critical_volume_loop(dive_state *dive);
void decompression_loop(dive_state *dive, single_dive *current_dive);
int set_gas_mixes(dive_state *dive, single_dive *current_dive);
void repetitive_dive_loop(dive_state *dive, json_input *input);
void free_dive_state(dive_state *dive);
void real_time_dive_loop(dive_state *dive, single_dive *current_dive);
direction current_direction(dive_state *dive, double increment_time);
double find_start_of_decompression_zone(dive_state *dive, dive_profile *current_profile);
int finished_constant_depth_profile(dive_state *dive, dive_profile *current_profile);
void calculate_decompression_stops(dive_state *dive, dive_profile *current_profile);
void critical_volume_decision_tree_to_depth(dive_state *dive, double stop_depth);
int critical_volume_loop_init(dive_state *dive);
void decompress_init(dive_state *dive, dive_profile *current_profile);
#endif
