import json
import subprocess
import os

VERBOSE=False

class CCE:
    def __init__(self, aflow_executable='aflow'):
        self.aflow_executable=aflow_executable

    def aflow_command(self, cmd):
        "checks output of aflow command"
        try:
          if VERBOSE: print('aflow_command(): cmd=' + self.aflow_executable + cmd)
          output=subprocess.check_output(
                  self.aflow_executable + cmd,
                  shell=True
                  )
          return output
        except subprocess.CalledProcessError:
          raise OSError('aflow executable not found: ' + self.aflow_executable)

    def add_path_to_structure(self, cmd, struct_file):
        "adds path to the structure file to the command"
        if (os.path.exists(struct_file)):
          cmd += ' < ' + struct_file
          return cmd
        else:
          raise OSError(struct_file + ' not found')

    def get_json(self, cmd):
        "sets json as output format, executes aflow command, and checks and returns json"
        #return json output by default
        cmd += " --print=json"
        #execute command
        output=self.aflow_command(cmd)
        #validate json
        try:
          res_json=json.loads(output)
        except:
          print(output)
          raise RuntimeError('The output does not seem to be a valid json.')
        return res_json

    def validate_functionals(self, available_functionals, input_functionals):
        "checks whether for all given functionals corrections are available"
        invalid_functionals=[n for n in input_functionals if n not in available_functionals]
        if invalid_functionals:
          raise ValueError("No corrections available for: " + ",".join(invalid_functionals))

    def add_functionals(self, functionals, cmd):
        "adds functionals to command"
        functionals_str=",".join(functionals)
        cmd += " --functionals=" + functionals_str
        return cmd

    def get_corrections(self, struct_file, enthalpies_formation_dft=[], functionals=[], oxidation_numbers=[]):
        "determines CCE corrections and corrected formation enthalpies if precalculated DFT values are given"
        #initial cce command
        command=' --get_cce_corrections'

        #adding path to structure file
        command=self.add_path_to_structure(command, struct_file)

        #handling functionals
        available_functionals=["PBE", "LDA", "SCAN"]
        if type(functionals) == list and len(functionals) != 0:
          #validate
          self.validate_functionals(available_functionals, functionals)
          # add functionals
          command=self.add_functionals(functionals, command)
        elif type(functionals) == str and functionals:
          functionals=functionals.replace(" ", "")
          functionals=functionals.split(",")
          #validate
          self.validate_functionals(available_functionals, functionals)
          # add functionals
          command=self.add_functionals(functionals, command)

        #handling enthalpies_formation_dft
        if type(enthalpies_formation_dft) == list and len(enthalpies_formation_dft) != 0:
          if len(enthalpies_formation_dft) != len(functionals):
            raise ValueError("number of input formation enthalpies does not equal number of functionals")
          input_enthalpies_formation=",".join([str(n) for n in enthalpies_formation_dft])
          command += " --enthalpies_formation_dft=" + input_enthalpies_formation
        elif type(enthalpies_formation_dft) == str and enthalpies_formation_dft:
          if len(enthalpies_formation_dft.split(',')) != len(functionals):
            raise ValueError("number of input formation enthalpies does not equal number of functionals")
          input_enthalpies_formation=",".join([n for n in enthalpies_formation_dft.replace(" ", "").split(',')])
          command += " --enthalpies_formation_dft=" + input_enthalpies_formation

        #handling oxidation_numbers
        if type(oxidation_numbers) == list and len(oxidation_numbers) != 0:
          input_oxidation_numbers=",".join([str(n) for n in oxidation_numbers])
          command += " --oxidation_numbers=" + input_oxidation_numbers
        elif type(oxidation_numbers) == str and oxidation_numbers:
          input_oxidation_numbers=",".join([n for n in oxidation_numbers.replace(" ","").split(',')])
          command += " --oxidation_numbers=" + input_oxidation_numbers

        #get json
        return self.get_json(command)

    def get_oxidation_numbers(self, struct_file):
        "determines oxidation numbers"
        #initial cce command
        command=' --get_oxidation_numbers'

        #adding path to structure file
        command=self.add_path_to_structure(command, struct_file)

        #get json
        return self.get_json(command)

    def get_cation_coordination_numbers(self, struct_file):
        "determines cation coordination numbers"
        #initial cce command
        command=' --get_cation_coordination_numbers'

        #adding path to structure file
        command=self.add_path_to_structure(command, struct_file)

        #get json
        return self.get_json(command)
