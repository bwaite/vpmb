#!/usr/bin/env python
import json
import sys

from optparse import OptionParser

NUM_INVALID = 0

def check_range(value, start, end):
    if value >= start and value <= end:
        return True
    else:
        print "%s, is not in the range %s to %s" %(str(value), str(start), str(end))
        globals()['NUM_INVALID'] += 1
        return False
        
def check_tag(tag, expected_object):
    if isinstance(expected_object, list):
        is_type = False
        for _type in expected_object:
            if isinstance(tag, _type):
                is_type = True
                break
        if not is_type:
            print "%s, is not %s" %(tag,str(expected_object))
            globals()['NUM_INVALID'] += 1
            return False
        return True
    
    elif not isinstance(tag, expected_object):
        print "%s, is not %s" %(tag,str(expected_object))
        globals()['NUM_INVALID'] += 1
        return False
    return True

def check_ascent_summary(ascent):
    check_tag(ascent["starting_depth"], [float, int])
    check_tag(ascent["gasmix"], int)
    check_tag(ascent["rate"], [float, int])
    check_range(ascent["rate"], -sys.maxint, 0.0)
    check_tag(ascent["step_size"], [float, int])
    
    if ascent.has_key("setpoint"):
        check_tag(ascent["setpoint"], [float, int])
    
def check_profile(profile):
    if not check_tag(profile["profile_code"], int):
        print "profile_code must be 1, 2, or 99, not %d" %(profile["profile_code"])
        globals()['NUM_INVALID'] += 1
        return False
    
    if profile["profile_code"] == 1:
        check_tag(profile["starting_depth"], [float, int])
        check_tag(profile["ending_depth"], [float, int])
        check_tag(profile["rate"], [float, int])
        check_tag(profile["gasmix"], int)
        if profile.has_key("setpoint"):
            check_tag(profile["setpoint"], [float, int])

    elif profile["profile_code"] == 2:
        check_tag(profile["depth"], [float, int])
        check_tag(profile["run_time_at_end_of_segment"], [float, int])
        check_tag(profile["gasmix"], int)
        if profile.has_key("setpoint"):
            check_tag(profile["setpoint"], [float, int])

    elif profile["profile_code"] == 99:
        check_tag(profile["number_of_ascent_parameter_changes"], int)
        
        if len(profile["ascent_summary"]) != profile["number_of_ascent_parameter_changes"]:
            print "ascent_summary must be equal to the number of elements in the number_of_ascent_parameter_changes array"
            globals()['NUM_INVALID'] += 1

        for ascent in profile["ascent_summary"]:
            check_ascent_summary(ascent)

        
def validate_json(input_file_name):
    valid = True
    f = open(input_file_name)
    j = json.loads(f.read())
    f.close()

    input_settings = j["input"]
    altitude_settings = j["altitude"]
    settings = j["settings"]

    for dive in input_settings:
        check_tag(dive["desc"], unicode)
        check_tag(dive["num_gas_mixes"], int)
        check_tag(dive["gasmix_summary"], list)

        for gasmix in dive["gasmix_summary"]:
            check_tag(gasmix["fraction_O2"], float)
            check_tag(gasmix["fraction_He"], float)
            check_tag(gasmix["fraction_N2"], float)
            sum_of_gases = gasmix["fraction_O2"] + gasmix["fraction_He"] + gasmix["fraction_N2"]
            if sum_of_gases != 1.0:
                print "Gas mixes must add up to 1.0. Currently at %f" % (sum_of_gases)
                globals()['NUM_INVALID'] += 1
            
        if len(dive["gasmix_summary"]) != dive["num_gas_mixes"]:
            print "num_gas_mixes must be equal to the number of elements in the gasmix_summary array"
            globals()['NUM_INVALID'] += 1

        check_tag(dive["profile_codes"],list)
        for profile in dive["profile_codes"]:
            check_profile(profile)
            
        check_tag(dive["repetitive_code"], int)
        if profile.has_key("surface_interval_time_minutes"):
            check_tag(profile["surface_interval_time_minutes"], [float, int])

    # Altitude Settings
    check_tag(altitude_settings["Altitude_of_Dive"], [float, int])
    check_tag(altitude_settings["Diver_Acclimatized_at_Altitude"], unicode)
    check_tag(altitude_settings["Starting_Acclimatized_Altitude"], [float, int])
    check_tag(altitude_settings["Ascent_to_Altitude_Hours"], [float, int])
    check_tag(altitude_settings["Hours_at_Altitude_Before_Dive"], [float, int])

    # Settings
    check_tag(settings["Units"], unicode)
    check_tag(settings["Altitude_Dive_Algorithm"], unicode)
    
    check_tag(settings["Minimum_Deco_Stop_Time"], [float, int])
    check_range(settings["Minimum_Deco_Stop_Time"], 0.0, sys.maxint) # check positive

    check_tag(settings["Critical_Radius_N2_Microns"], [float, int])
    check_range(settings["Critical_Radius_N2_Microns"], 0.2, 1.35)
    
    check_tag(settings["Critical_Radius_He_Microns"], [float, int])
    check_range(settings["Critical_Radius_He_Microns"], 0.2, 1.35)
    
    check_tag(settings["Critical_Volume_Algorithm"], unicode)

    check_tag(settings["Crit_Volume_Parameter_Lambda"],[float, int])
    check_range(settings["Crit_Volume_Parameter_Lambda"], 6500, 8300)
    
    check_tag(settings["Gradient_Onset_of_Imperm_Atm"],[float, int])
    check_range(settings["Gradient_Onset_of_Imperm_Atm"], 5.0, 10.0)

    check_tag(settings["Surface_Tension_Gamma"],[float, int])
    check_range(settings["Surface_Tension_Gamma"], 0.015, 0.065)

    check_tag(settings["Skin_Compression_GammaC"],[float, int])
    check_range(settings["Skin_Compression_GammaC"], 0.160, 0.290)
    
    check_tag(settings["Regeneration_Time_Constant"],[float, int])
    check_range(settings["Regeneration_Time_Constant"],10080, 51840)
    
    check_tag(settings["Pressure_Other_Gases_mmHg"],[float, int])
    
    if NUM_INVALID == 0:
        print "%s is valid" %(input_file_name)
    else:
        print "%s has %d errors" %(input_file_name, NUM_INVALID)

      
if __name__ == '__main__':
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-i",  action="store", dest="input_file_name",
                      default="vpm_decompression_input.json",
                      help="Input file containing dive information")

    (options, args) = parser.parse_args()
    validate_json(options.input_file_name)
