import json
import subprocess
import os

class XtalFinder:

    def __init__(self, aflow_executable='aflow'):
        self.aflow_executable = aflow_executable

    def aflow_command(self, cmd):
        try:
            return subprocess.check_output(
                self.aflow_executable + cmd,
                shell=True
            )
        except subprocess.CalledProcessError:
            print "Error aflow executable not found at: " + self.aflow_executable

    def get_prototype_label(self, input_file, options=None):
        fpath = os.path.realpath(input_file)
        command = ' --prototype'
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json < ' + fpath
        )

        res_json = json.loads(output)
        return res_json

    def compare_materials(self, input_files, options=None):
        command = ' --compare_materials=' + ','.join(input_files)
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )

        res_json = json.loads(output)
        return res_json

    def compare_materials_directory(self, directory, options=None):
        command = ' --compare_materials -D ' + directory
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )

        res_json = json.loads(output)
        return res_json

    def compare_materials_file(self, filename, options=None):
        command = ' --compare_materials -F=' + filename
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )

        res_json = json.loads(output)
        return res_json

    def compare_structures(self, input_files, options=None):
        command = ' --compare_structures=' + ','.join(input_files)
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )

        res_json = json.loads(output)
        return res_json

    def compare_structures_directory(self, directory, options=None):
        command = ' --compare_structures -D ' + directory
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )

        res_json = json.loads(output)
        return res_json

    def compare_structures_file(self, filename, options=None):
        command = ' --compare_structures -F=' + filename
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet'
        )

        res_json = json.loads(output)
        return res_json

    def compare2database(self, input_file, options=None):
        fpath = os.path.realpath(input_file)
        command = ' --compare2database'
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet < ' + fpath
        )

        res_json = json.loads(output)
        return res_json

    def compare2prototypes(self, input_file, options=None):
        fpath = os.path.realpath(input_file)
        command = ' --compare2prototypes'
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --quiet < ' + fpath
        )

        res_json = json.loads(output)
        return res_json

    def get_isopointal_prototypes(self, input_file, options=None):
        fpath = os.path.realpath(input_file)
        command = ' --isopointal_prototype'
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json < ' + fpath
        )

        res_json = json.loads(output)
        return res_json

    def get_unique_atom_decorations(self, input_file, options=None):
        fpath = os.path.realpath(input_file)
        command = ' --unique_atom_decorations'
        output = ''

        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json < ' + fpath
        )

        res_json = json.loads(output)
        return res_json
