import json 
import subprocess 
import os 
 
class Symmetry: 

    def __init__(self, aflow_executable='aflow'): 
        self.aflow_executable = aflow_executable 


    def aflow_command(self, cmd): 
        try: 
            return subprocess.check_output( 
                self.aflow_executable + cmd, 
                shell=True 
            ) 
        except subprocess.CalledProcessError: 
            raise AssertionError('aflow executable not found: ' + self.aflow_executable)

    def get_filepath(self, filename):
        if os.path.exists(input_file.name):
            return os.path.realpath(filename)
        else:
            raise OSError(filename + ' not found')

    def get_symmetry(self, input_file, tol=None, magmoms=None):
        fpath = get_filepath(input_file.name)
        command = ' --aflowSYM' 
        output = '' 


        if tol: 
            command += '=' + str(tol) 
        if magmoms: 
            command += ' --magmom=' + magmoms 

        output = self.aflow_command( 
            command + ' --print=json --screen_only' + ' < ' + fpath 
        ) 
        res_json = json.loads(output) 
        return res_json 


    def get_edata(self, input_file, tol=None, magmoms=None): 
        fpath = get_filepath(input_file.name)
        command = ' --edata' 
        output = '' 


        if tol: 
            command += '=' + str(tol) 
        if magmoms: 
            command += ' --magmom=' + magmoms 

        output = self.aflow_command( 
            command + ' --print=json' + ' < ' + fpath 
        ) 
        res_json = json.loads(output) 
        return res_json 

    def get_sgdata(self, input_file, tol=None, magmoms=None): 
        fpath = get_filepath(input_file.name)
        command = ' --sgdata' 
        output = '' 


        if tol: 
            command += '=' + str(tol) 
        if magmoms: 
            command += ' --magmom=' + magmoms 

        output = self.aflow_command( 
            command + ' --print=json' + ' < ' + fpath 
        ) 
        res_json = json.loads(output) 
        return res_json 
