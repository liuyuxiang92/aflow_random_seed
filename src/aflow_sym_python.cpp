#ifndef _AFLOW_SYM_PYTHON_CPP_
#define _AFLOW_SYM_PYTHON_CPP_
// aflow_sym_python.cpp automatic generated
std::string AFLOW_SYM_PYTHON_PY="\
import json  \n\
import subprocess  \n\
import os  \n\
  \n\
class Symmetry:  \n\
 \n\
    def __init__(self, aflow_executable='aflow'):  \n\
        self.aflow_executable = aflow_executable  \n\
 \n\
 \n\
    def aflow_command(self, cmd):  \n\
        try:  \n\
            return subprocess.check_output(  \n\
                self.aflow_executable + cmd,  \n\
                shell=True  \n\
            )  \n\
        except subprocess.CalledProcessError:  \n\
            print('Error aflow executable not found at: ' + self.aflow_executable)  \n\
 \n\
 \n\
    def get_symmetry(self, input_file, tol=None, magmoms=None):  \n\
        fpath = os.path.realpath(input_file.name)  \n\
        command = ' --aflowSYM'  \n\
        output = ''  \n\
 \n\
 \n\
        if tol:  \n\
            command += '=' + str(tol)  \n\
        if magmoms:  \n\
            command += ' --magmom=' + magmoms  \n\
 \n\
        output = self.aflow_command(  \n\
            command + ' --print=json --screen_only' + ' < ' + fpath  \n\
        )  \n\
        res_json = json.loads(output)  \n\
        return res_json  \n\
 \n\
 \n\
    def get_edata(self, input_file, tol=None, magmoms=None):  \n\
        fpath = os.path.realpath(input_file.name)  \n\
        command = ' --edata'  \n\
        output = ''  \n\
 \n\
 \n\
        if tol:  \n\
            command += '=' + str(tol)  \n\
        if magmoms:  \n\
            command += ' --magmom=' + magmoms  \n\
 \n\
        output = self.aflow_command(  \n\
            command + ' --print=json' + ' < ' + fpath  \n\
        )  \n\
        res_json = json.loads(output)  \n\
        return res_json  \n\
 \n\
    def get_sgdata(self, input_file, tol=None, magmoms=None):  \n\
        fpath = os.path.realpath(input_file.name)  \n\
        command = ' --sgdata'  \n\
        output = ''  \n\
 \n\
 \n\
        if tol:  \n\
            command += '=' + str(tol)  \n\
        if magmoms:  \n\
            command += ' --magmom=' + magmoms  \n\
 \n\
        output = self.aflow_command(  \n\
            command + ' --print=json' + ' < ' + fpath  \n\
        )  \n\
        res_json = json.loads(output)  \n\
        return res_json  \n\
";
#endif // _AFLOW_SYM_PYTHON_CPP_
