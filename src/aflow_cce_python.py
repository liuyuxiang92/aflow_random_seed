#!/usr/bin/env python

import json
import subprocess
import os

VERBOSE=False

class CCE:
    def __init__(self, aflow_executable = 'aflow'):
        self.aflow_executable = aflow_executable

    def aflow_command(self, cmd):
        try:
            if VERBOSE: print('aflow_command(): cmd = ' + self.aflow_executable + cmd)
            output = subprocess.check_output(
                    self.aflow_executable + cmd,
                    shell=True
                    )
            return output
        except subprocess.CalledProcessError:
            raise AssertionError('aflow executable not found: ' + self.aflow_executable)

    def get_corrections(self, STRUCT_FILE, dft_formation_energies = [], functionals = [], oxidation_numbers = []):
        #initial cce command
        command = ' --cce'

        #adding path to structure file
        if (os.path.exists(STRUCT_FILE)):
          command += '=' + STRUCT_FILE
        else:
          raise OSError('STRUCT_FILE' + ' not found')

        #handling functionals
        AVAILABLE_FUNCTIONALS = ["PBE", "LDA", "SCAN"]

        if len(functionals) == 0:
          command += ""
          input_functionals = []
        elif type(functionals) == list:
          #validate
          invalid_functionals = [n for n in functionals if n not in AVAILABLE_FUNCTIONALS]
          if len(invalid_functionals) > 0:
            raise AssertionError(",".join(invalid_functionals) + " functionals are not available")

          # add input_functionals
          input_functionals = functionals
          input_functionals_str = ",".join(input_functionals)
          command += " --functionals=" + input_functionals_str
        elif type(functionals) == str:
          input_functionals = functionals.replace(" ", "")
          if len(input_functionals) == 0:
            command += ""
          input_functionals = input_functionals.split(",")
          #validate
          invalid_functionals = [n for n in input_functionals if n not in AVAILABLE_FUNCTIONALS]
          if len(invalid_functionals) > 0:
            raise AssertionError(",".join(invalid_functionals) + " functionals are not available")
          # add input_functionals
          input_functionals_str = ",".join(input_functionals)
          command += " --functionals=" + input_functionals_str
        # else:
        #   command += " --functionals=" + ",".join(AVAILABLE_FUNCTIONALS)



        #handling dft_formation_energies
        if type(dft_formation_energies) == list:
          if len(dft_formation_energies) != 0:
            if len(dft_formation_energies) != len(input_functionals):
              raise AssertionError("number of input formation energies does not equal number of functionals")

            input_formation_energies = ",".join([str(n) for n in dft_formation_energies])
            command += " --dft_formation_energies=" + input_formation_energies

        elif type(dft_formation_energies) == str:
          if len(dft_formation_energies) != 0:
            if len(dft_formation_energies.split(',')) != len(input_functionals):
              raise AssertionError("number of input formation energies does not equal number of functionals")

            input_formation_energies = ",".join([n for n in dft_formation_energies.replace(" ", "").split(',')])
            command += " --dft_formation_energies=" + input_formation_energies



        #handling oxidation_numbers
        if type(oxidation_numbers) == list:
          if len(oxidation_numbers) != 0:
            input_oxidation_numbers = ",".join([str(n) for n in oxidation_numbers])
            command += " --oxidation_numbers=" + input_oxidation_numbers
        elif type(oxidation_numbers) == str:
          if len(oxidation_numbers) != 0:
            input_oxidation_numbers = ",".join([n for n in oxidation_numbers.replace(" ","").split(',')])
            command += " --oxidation_numbers=" + input_oxidation_numbers

        #return json output by default
        command += " --print=json"

        #execute command
        output = self.aflow_command(command)

        try:
          res_json = json.loads(output)
        except:
          res_json = output

        return res_json
