#!/usr/bin/env python

import unittest
import vpmb
import json

def is_json_equal(json1, json2):
    """Checks the output json is equal (skips time and other unimportant things)"""
    if len(json1) != len(json2):
        return False
    
    for i, dive in enumerate(json1):
        dive2 = json2[i]

        if dive["dive_profile"] != dive2["dive_profile"]:
            return False
        if dive["gasmix"] != dive2["gasmix"]:
            return False
        if dive["decompression_profile"] != dive2["decompression_profile"]:
            return False

    return True

def load_and_check_file(input_file, expected_results):
    """Load the input file, run the decompression algorithm, and return 'expected'
    and 'results' json files"""
    f = open(input_file)
    j = json.loads(f.read())
    f.close()
    
    program_state = vpmb.DiveState(json_input=j)
    program_state.main()
    results = program_state.get_json()
    
    f = open(expected_results)
    expected = json.loads(f.read())
    f.close()

    return expected, results

def load_and_return_program_state(input_file):
    """Load the input file, create the State object, and return it"""
    f = open(input_file)
    j = json.loads(f.read())
    f.close()
    
    return vpmb.DiveState(json_input=j)
    
class TestDiveResults(unittest.TestCase):
    
    def test_test1(self):
        expected, results = load_and_check_file("test1/vpm_decompression_input.json", "test1/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_msw_simple(self):
        expected, results = load_and_check_file("msw_simple/vpm_decompression_input.json", "msw_simple/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_msw_test(self):
        expected, results = load_and_check_file("msw_test/vpm_decompression_input.json", "msw_test/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_paper_test(self):
        expected, results = load_and_check_file("paper_test/vpm_decompression_input.json", "paper_test/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_small_steps(self):
        expected, results = load_and_check_file("small_steps/vpm_decompression_input.json", "small_steps/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_small_steps_meters(self):
        expected, results = load_and_check_file("small_steps_meters/vpm_decompression_input.json", "small_steps_meters/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_no_deco_stop(self):
        expected, results = load_and_check_file("no_deco_stop/vpm_decompression_input.json", "no_deco_stop/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_gradient_onset_of_imperm_greater(self):
        expected, results = load_and_check_file("gradient_onset_of_imperm_greater/vpm_decompression_input.json", "gradient_onset_of_imperm_greater/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_not_acclimatized(self):
        expected, results = load_and_check_file("not_acclimatized/vpm_decompression_input.json", "not_acclimatized/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_gradient_N2_bubble_formation(self):
        expected, results = load_and_check_file("gradient_N2_bubble_formation/vpm_decompression_input.json", "gradient_N2_bubble_formation/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_altitude_dive1(self):
        expected, results = load_and_check_file("altitude_dive1/vpm_decompression_input.json", "altitude_dive1/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_smooth_dive(self):
        expected, results = load_and_check_file("smooth_dive/vpm_decompression_input.json", "smooth_dive/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_multi_constant_depth(self):
        expected, results = load_and_check_file("multi_constant_depth/vpm_decompression_input.json", "multi_constant_depth/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_multiple_descents_and_ascents(self):
        expected, results = load_and_check_file("multiple_descents_and_ascents/vpm_decompression_input.json", "multiple_descents_and_ascents/expected.json")
        self.assertTrue(is_json_equal(expected, results))

    def test_load_from_file(self):

        program_state = vpmb.DiveState(input_file_name="test1/vpm_decompression_input.json")
        program_state.main()
        results = program_state.get_json()
    
        f = open("test1/expected.json")
        expected = json.loads(f.read())
        f.close()

        self.assertTrue(is_json_equal(expected, results))

class TestAltitudeConfigurations(unittest.TestCase):
    def test_altitude_exception1(self):
        program_state = load_and_return_program_state("altitude_exception1/vpm_decompression_input.json")
        self.assertRaises(vpmb.AltitudeException, program_state.main)
                          
    def test_altitude_exception2(self):
        program_state = load_and_return_program_state("altitude_exception2/vpm_decompression_input.json")
        self.assertRaises(vpmb.AltitudeException, program_state.main)

    def test_altitude_exception3(self):
        program_state = load_and_return_program_state("altitude_exception3/vpm_decompression_input.json")
        self.assertRaises(vpmb.AltitudeException, program_state.main)

    def test_altitude_exception4(self):
        program_state = load_and_return_program_state("altitude_exception4/vpm_decompression_input.json")
        self.assertRaises(vpmb.AltitudeException, program_state.main)

    def test_altitude_exception5(self):
        program_state = load_and_return_program_state("altitude_exception5/vpm_decompression_input.json")
        self.assertRaises(vpmb.AltitudeException, program_state.main)

    def test_altitude_exception6(self):
        program_state = load_and_return_program_state("altitude_exception6/vpm_decompression_input.json")
        self.assertRaises(vpmb.AltitudeException, program_state.main)

class TestGeneralExceptions(unittest.TestCase):
    
    def test_no_input_exception(self):
        self.assertRaises(ValueError, vpmb.DiveState)

    def test_gasmix_exception1(self):
        program_state = load_and_return_program_state("gasmix_exception1/vpm_decompression_input.json")
        self.assertRaises(vpmb.InputFileException, program_state.main)

    def test_gasmix_exception2(self):
        program_state = load_and_return_program_state("gasmix_exception2/vpm_decompression_input.json")
        self.assertRaises(vpmb.InputFileException, program_state.main)

    def test_profile_code_exception1(self):
        program_state = load_and_return_program_state("profile_code_exception1/vpm_decompression_input.json")
        self.assertRaises(vpmb.InputFileException, program_state.main)

    def test_repetative_dive_code_exception1(self):
        program_state = load_and_return_program_state("repetative_dive_code_exception1/vpm_decompression_input.json")
        self.assertRaises(vpmb.InputFileException, program_state.main)


if __name__ == '__main__':
    unittest.main()
