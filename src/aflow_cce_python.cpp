#ifndef _AFLOW_CCE_PYTHON_CPP_
#define _AFLOW_CCE_PYTHON_CPP_
// aflow_cce_python.cpp automatic generated
std::string AFLOW_CCE_PYTHON_PY="\
import json \n\
import os \n\
import subprocess \n\
import warnings \n\
 \n\
VERBOSE = False \n\
 \n\
 \n\
class CCE: \n\
    def __init__(self, aflow_executable='aflow'): \n\
        '''Class implementing CCE python functionality \n\
 \n\
        Methods: \n\
        get_corrections -- get CCE corrections \n\
        Arguments: \n\
        struct_file_path -- path to the structure file \n\
        functionals -- functionals to get corrections for (optional) \n\
        enthalpies_formation_dft -- input DFT formation enthalpies (optional) \n\
        oxidation_numbers -- input oxidation numbers for each atom (optional) \n\
 \n\
        get_oxidation_numbers -- get CCE oxidation numbers \n\
        Arguments: \n\
        struct_file_path -- path to the structure file \n\
 \n\
        get_cation_coordination_numbers -- get CCE cation coordination numbers \n\
        Arguments: \n\
        struct_file_path -- path to the structure file \n\
        ''' \n\
 \n\
        self.aflow_executable = aflow_executable \n\
        # checks whether aflow executable exists \n\
        try: \n\
            subprocess.check_output(self.aflow_executable, shell=True) \n\
        except subprocess.CalledProcessError: \n\
            raise OSError('aflow executable not found: ' + \n\
                          self.aflow_executable) \n\
 \n\
    def run_aflow_command(self, cmd): \n\
        '''checks output of aflow command''' \n\
        try: \n\
            if VERBOSE: \n\
                print('aflow_command(): cmd=' + self.aflow_executable + cmd) \n\
            output = subprocess.check_output(self.aflow_executable + cmd, \n\
                                             shell=True) \n\
            return output \n\
        except subprocess.CalledProcessError: \n\
            raise OSError('An error occurred while executing: ' + \n\
                          self.aflow_executable + cmd) \n\
 \n\
    def add_structure_file_path(self, cmd, struct_file): \n\
        '''adds path to the structure file to the command''' \n\
        if (os.path.exists(struct_file)): \n\
            cmd += ' < ' + struct_file \n\
            return cmd \n\
        else: \n\
            raise OSError(struct_file + ' not found') \n\
 \n\
    def get_json(self, cmd): \n\
        '''sets json as output format, executes aflow command, \n\
        and checks and returns json''' \n\
        # return json output by default \n\
        cmd += ' --print=json' \n\
        # execute command \n\
        output = self.run_aflow_command(cmd) \n\
        # validate json \n\
        try: \n\
            res_json = json.loads(output) \n\
        except json.decoder.JSONDecodeError: \n\
            print(output) \n\
            raise RuntimeError('Output is not valid json.') \n\
        return res_json \n\
 \n\
    def validate_functionals(self, input_functionals): \n\
        '''checks whether for all given functionals corrections are \n\
        available''' \n\
        available_functionals = ['PBE', 'LDA', 'SCAN'] \n\
        invalid_functionals = [n for n in input_functionals if n not in \n\
                               available_functionals] \n\
        valid_functionals = [n for n in input_functionals if n in \n\
                             available_functionals] \n\
        if invalid_functionals: \n\
            warnings.warn('No corrections available for: ' + \n\
                          ','.join(invalid_functionals)) \n\
        if not valid_functionals: \n\
            raise ValueError('No valid functionals provided. ' \n\
                             'Please choose from: ' + \n\
                             ','.join(available_functionals)) \n\
        return valid_functionals \n\
 \n\
    def add_functionals(self, functionals, cmd): \n\
        '''adds functionals to command''' \n\
        functionals_str = ','.join(functionals) \n\
        cmd += ' --functionals=' + functionals_str \n\
        return cmd \n\
 \n\
    def get_corrections(self, struct_file, enthalpies_formation_dft=[], \n\
                        functionals=[], oxidation_numbers=[]): \n\
        '''determines CCE corrections and corrected formation enthalpies if \n\
        precalculated DFT values are given''' \n\
        # initial cce command \n\
        command = ' --get_cce_corrections' \n\
 \n\
        # adding path to structure file \n\
        command = self.add_structure_file_path(command, struct_file) \n\
 \n\
        # handling functionals \n\
        if type(functionals) == str: \n\
            if functionals: \n\
                # convert to list \n\
                functionals = functionals.replace(' ', '') \n\
                functionals = functionals.split(',') \n\
        elif type(functionals) == tuple: \n\
            if len(functionals) > 0: \n\
                # convert to list \n\
                functionals = list(functionals) \n\
        elif type(functionals) == list: \n\
            pass \n\
        else: \n\
            raise TypeError('The functionals argument must be either a list, ' \n\
                            'a tuple, or a string with functionals ' \n\
                            'separated by commas.') \n\
        if len(functionals) > 0: \n\
            functionals = [x.upper() for x in functionals] \n\
            valid_functionals = self.validate_functionals(functionals) \n\
            command = self.add_functionals(valid_functionals, command) \n\
 \n\
        # handling enthalpies_formation_dft \n\
        if type(enthalpies_formation_dft) == str: \n\
            if enthalpies_formation_dft: \n\
                # convert to list \n\
                enthalpies_formation_dft = enthalpies_formation_dft.replace( \n\
                                           ' ', '') \n\
                enthalpies_formation_dft = enthalpies_formation_dft.split(',') \n\
        elif type(enthalpies_formation_dft) == tuple: \n\
            if len(enthalpies_formation_dft) > 0: \n\
                # convert to list \n\
                enthalpies_formation_dft = list(enthalpies_formation_dft) \n\
        elif type(enthalpies_formation_dft) == list: \n\
            pass \n\
        else: \n\
            raise TypeError('The enthalpies argument must be either a list, ' \n\
                            'a tuple, or a string with enthalpies ' \n\
                            'separated by commas.') \n\
        if len(enthalpies_formation_dft) > 0: \n\
            if len(enthalpies_formation_dft) != len(functionals): \n\
                raise ValueError('number of input formation enthalpies ' \n\
                                 'does not equal number of functionals') \n\
            # consider only enthalpies corresponding to valid functionals \n\
            if len(enthalpies_formation_dft) != len(valid_functionals): \n\
                enthalpies_new = [] \n\
                for i in range(len(functionals)): \n\
                    for j in range(len(valid_functionals)): \n\
                        if functionals[i] == valid_functionals[j]: \n\
                            enthalpies_new.append(enthalpies_formation_dft[i]) \n\
                enthalpies_formation_dft = enthalpies_new \n\
            enthalpies_formation = ','.join([str(n) for n in \n\
                                            enthalpies_formation_dft]) \n\
            command += ' --enthalpies_formation_dft=' + enthalpies_formation \n\
 \n\
        # handling oxidation_numbers \n\
        if type(oxidation_numbers) == str: \n\
            if oxidation_numbers: \n\
                # convert to list \n\
                oxidation_numbers = oxidation_numbers.replace(' ', '') \n\
                oxidation_numbers = oxidation_numbers.split(',') \n\
        elif type(oxidation_numbers) == tuple: \n\
            if len(oxidation_numbers) > 0: \n\
                # convert to list \n\
                oxidation_numbers = list(oxidation_numbers) \n\
        elif type(oxidation_numbers) == list: \n\
            pass \n\
        else: \n\
            raise TypeError('The oxidation numbers must be either a list, ' \n\
                            'a tuple, or a string with numbers ' \n\
                            'separated by commas.') \n\
        if len(oxidation_numbers) > 0: \n\
            input_oxidation_numbers = ','.join([str(n) for n in \n\
                                               oxidation_numbers]) \n\
            command += ' --oxidation_numbers=' + input_oxidation_numbers \n\
 \n\
        # get json \n\
        return self.get_json(command) \n\
 \n\
    def get_oxidation_numbers(self, struct_file): \n\
        '''determines oxidation numbers''' \n\
        # initial cce command \n\
        command = ' --get_oxidation_numbers' \n\
 \n\
        # adding path to structure file \n\
        command = self.add_structure_file_path(command, struct_file) \n\
 \n\
        # get json \n\
        return self.get_json(command) \n\
 \n\
    def get_cation_coordination_numbers(self, struct_file): \n\
        '''determines cation coordination numbers''' \n\
        # initial cce command \n\
        command = ' --get_cation_coordination_numbers' \n\
 \n\
        # adding path to structure file \n\
        command = self.add_structure_file_path(command, struct_file) \n\
 \n\
        # get json \n\
        return self.get_json(command) \n\
";
#endif // _AFLOW_CCE_PYTHON_CPP_
