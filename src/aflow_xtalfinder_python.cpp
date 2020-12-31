#ifndef _AFLOW_XTALFINDER_PYTHON_CPP_
#define _AFLOW_XTALFINDER_PYTHON_CPP_
// aflow_xtalfinder_python.cpp automatic generated
std::string AFLOW_XTALFINDER_PYTHON_PY="\
import json \n\
import subprocess \n\
import os \n\
 \n\
class XtalFinder: \n\
 \n\
    def __init__(self, aflow_executable='aflow'): \n\
        self.aflow_executable = aflow_executable \n\
 \n\
    def aflow_command(self, cmd): \n\
        try: \n\
            return subprocess.check_output( \n\
                self.aflow_executable + cmd, \n\
                shell=True \n\
            ) \n\
        except subprocess.CalledProcessError: \n\
            print \"Error aflow executable not found at: \" + self.aflow_executable \n\
 \n\
    def get_prototype_label(self, input_file, options=None): \n\
        fpath = os.path.realpath(input_file) \n\
        command = ' --prototype' \n\
        output = '' \n\
 \n\
        if options: \n\
            command += ' ' + options \n\
 \n\
        output = self.aflow_command( \n\
            command + ' --print=json < ' + fpath \n\
        ) \n\
 \n\
        res_json = json.loads(output) \n\
        return res_json \n\
 \n\
    def compare_materials(self, input_files, options=None): \n\
        command = ' --compare_materials=' + ','.join(input_files) \n\
        output = '' \n\
 \n\
        if options: \n\
            command += ' ' + options \n\
 \n\
        output = self.aflow_command( \n\
            command + ' --print=json --screen_only --quiet' \n\
        ) \n\
 \n\
        res_json = json.loads(output) \n\
        return res_json \n\
 \n\
    def compare_materials_directory(self, directory, options=None): \n\
        command = ' --compare_materials -D ' + directory \n\
        output = '' \n\
 \n\
        if options: \n\
            command += ' ' + options \n\
 \n\
        output = self.aflow_command( \n\
            command + ' --print=json --screen_only --quiet' \n\
        ) \n\
 \n\
        res_json = json.loads(output) \n\
        return res_json \n\
 \n\
    def compare_materials_file(self, filename, options=None): \n\
        command = ' --compare_materials -F=' + filename \n\
        output = '' \n\
 \n\
        if options: \n\
            command += ' ' + options \n\
 \n\
        output = self.aflow_command( \n\
            command + ' --print=json --screen_only --quiet' \n\
        ) \n\
 \n\
        res_json = json.loads(output) \n\
        return res_json \n\
 \n\
    def compare_structures(self, input_files, options=None): \n\
        command = ' --compare_structures=' + ','.join(input_files) \n\
        output = '' \n\
 \n\
        if options: \n\
            command += ' ' + options \n\
 \n\
        output = self.aflow_command( \n\
            command + ' --print=json --screen_only --quiet' \n\
        ) \n\
 \n\
        res_json = json.loads(output) \n\
        return res_json \n\
 \n\
    def compare_structures_directory(self, directory, options=None): \n\
        command = ' --compare_structures -D ' + directory \n\
        output = '' \n\
 \n\
        if options: \n\
            command += ' ' + options \n\
 \n\
        output = self.aflow_command( \n\
            command + ' --print=json --screen_only --quiet' \n\
        ) \n\
 \n\
        res_json = json.loads(output) \n\
        return res_json \n\
 \n\
    def compare_structures_file(self, filename, options=None): \n\
        command = ' --compare_structures -F=' + filename \n\
        output = '' \n\
 \n\
        if options: \n\
            command += ' ' + options \n\
 \n\
        output = self.aflow_command( \n\
            command + ' --print=json --screen_only --quiet' \n\
        ) \n\
 \n\
        res_json = json.loads(output) \n\
        return res_json \n\
 \n\
    def compare2database(self, input_file, options=None): \n\
        fpath = os.path.realpath(input_file) \n\
        command = ' --compare2database' \n\
        output = '' \n\
 \n\
        if options: \n\
            command += ' ' + options \n\
 \n\
        output = self.aflow_command( \n\
            command + ' --print=json --screen_only --quiet < ' + fpath \n\
        ) \n\
 \n\
        res_json = json.loads(output) \n\
        return res_json \n\
 \n\
    def compare2prototypes(self, input_file, options=None): \n\
        fpath = os.path.realpath(input_file) \n\
        command = ' --compare2prototypes' \n\
        output = '' \n\
 \n\
        if options: \n\
            command += ' ' + options \n\
 \n\
        output = self.aflow_command( \n\
            command + ' --print=json --screen_only --quiet < ' + fpath \n\
        ) \n\
 \n\
        res_json = json.loads(output) \n\
        return res_json \n\
 \n\
    def get_isopointal_prototypes(self, input_file, options=None): \n\
        fpath = os.path.realpath(input_file) \n\
        command = ' --isopointal_prototype' \n\
        output = '' \n\
 \n\
        if options: \n\
            command += ' ' + options \n\
 \n\
        output = self.aflow_command( \n\
            command + ' --print=json < ' + fpath \n\
        ) \n\
 \n\
        res_json = json.loads(output) \n\
        return res_json \n\
 \n\
    def get_unique_atom_decorations(self, input_file, options=None): \n\
        fpath = os.path.realpath(input_file) \n\
        command = ' --unique_atom_decorations' \n\
        output = '' \n\
 \n\
        if options: \n\
            command += ' ' + options \n\
 \n\
        output = self.aflow_command( \n\
            command + ' --print=json < ' + fpath \n\
        ) \n\
 \n\
        res_json = json.loads(output) \n\
        return res_json \n\
";
#endif // _AFLOW_XTALFINDER_PYTHON_CPP_
