#!/usr/bin/env python
import json
from optparse import OptionParser

def convert_json(input_file_name):
    f = open(input_file_name)
    j = json.loads(f.read())
    f.close()
    
    input_settings = j["input"]
    altitude_settings = j["altitude"]
    settings = j["settings"]
    
    # output altitude
    f = open("ALTITUDE.SET", "w")
    f.write("&Altitude_Dive_Settings\n")
    f.write("Altitude_of_Dive=%f\n"% (altitude_settings["Altitude_of_Dive"]))
    f.write("Diver_Acclimatized_at_Altitude='%s'\n"% (altitude_settings["Diver_Acclimatized_at_Altitude"]))
    f.write("Starting_Acclimatized_Altitude=%f\n"% (altitude_settings["Starting_Acclimatized_Altitude"]))
    f.write("Ascent_to_Altitude_Hours=%f\n"% (altitude_settings["Ascent_to_Altitude_Hours"]))
    f.write("Hours_at_Altitude_Before_Dive=%f\n"% (altitude_settings["Hours_at_Altitude_Before_Dive"]))
    f.write("/")
    
    f.close()
    
    # output settings
    
    f = open("VPMDECO.SET", "w")
    
    f.write("&Program_Settings\n")
    f.write("Units='%s'\n"%(settings["Units"]))
    f.write("Altitude_Dive_Algorithm='%s'\n"%(settings["Altitude_Dive_Algorithm"]))
    f.write("Minimum_Deco_Stop_Time=%f\n"%(settings["Minimum_Deco_Stop_Time"]))
    f.write("Critical_Radius_N2_Microns=%f\n"%(settings["Critical_Radius_N2_Microns"]))
    f.write("Critical_Radius_He_Microns=%f\n"%(settings["Critical_Radius_He_Microns"]))
    f.write("Critical_Volume_Algorithm='%s'\n"%(settings["Critical_Volume_Algorithm"]))
    f.write("Crit_Volume_Parameter_Lambda=%f\n"%(settings["Crit_Volume_Parameter_Lambda"]))
    f.write("Gradient_Onset_of_Imperm_Atm=%f\n"%(settings["Gradient_Onset_of_Imperm_Atm"]))
    f.write("Surface_Tension_Gamma=%f\n"%(settings["Surface_Tension_Gamma"]))
    f.write("Skin_Compression_GammaC=%f\n"%(settings["Skin_Compression_GammaC"]))
    f.write("Regeneration_Time_Constant=%f\n"%(settings["Regeneration_Time_Constant"]))
    f.write("Pressure_Other_Gases_mmHg=%f\n"%(settings["Pressure_Other_Gases_mmHg"]))
    f.write("/")
    
    f.close()
    
    # output dive
    
    f = open("VPMDECO.IN", "w")
    
    for dive in input_settings:
        f.write(dive["desc"] + "\n")
        f.write(str(len(dive["gasmix_summary"])) + "\n")
        for gasmix in dive["gasmix_summary"]:
            f.write(str(gasmix["fraction_O2"]) + ",")
            f.write(str(gasmix["fraction_He"]) + ",")
            f.write(str(gasmix["fraction_N2"]) + "\n")
    
        for profile in dive["profile_codes"]:
            f.write(str(profile["profile_code"]) + "\n")
            if profile["profile_code"] == 1:
                f.write(str(profile["starting_depth"]) + "," + str(profile["ending_depth"])+ "," + str(profile["rate"])+ "," + str(profile["gasmix"])+ "," + str(profile["setpoint"]) + "\n")
            elif profile["profile_code"] == 2:
                f.write(str(profile["depth"])+ "," + str(profile["run_time_at_end_of_segment"])+ "," + str(profile["gasmix"])+ "," + str(profile["setpoint"]) +"\n")
            elif profile["profile_code"] == 99:
                f.write(str(len(profile["ascent_summary"]))+"\n")
                for ascent in profile["ascent_summary"]:
                    f.write(str(ascent["starting_depth"])+ "," + str(ascent["gasmix"])+ "," + str(ascent["rate"])+ "," + str(ascent["step_size"])+ "," + str(ascent["setpoint"]) +"\n")

    
        f.write(str(dive["repetitive_code"])+"\n")
        if dive["repetitive_code"] == 1:
            f.write(str(dive["surface_interval_time_minutes"])+"\n")
    
    
    f.close()


if __name__ == '__main__':
    usage = "usage: %prog [options]"
    parser = OptionParser(usage=usage)

    parser.add_option("-i",  action="store", dest="input_file_name",
                      default="vpm_decompression_input.json",
                      help="Input file containing dive information")

    (options, args) = parser.parse_args()
    convert_json(options.input_file_name)
