#!/usr/bin/env python

import json
import subprocess
import os

class CHull:

    def __init__(self, aflow_executable='aflow'):
        self.aflow_executable = aflow_executable

    def aflow_command(self, cmd):
        try:
            return subprocess.check_output(
                self.aflow_executable + cmd,
                shell=True
            )
        except subprocess.CalledProcessError:
            print('Error aflow executable not found at: ' + self.aflow_executable)

    def get_hull(self, alloy, options = None):
        command = ' --chull --force'
        if options:
            command += ' ' + options

        output = self.aflow_command(
            command + ' --print=json --screen_only --alloy=' + alloy
        )
        #print('output: ' + output)
        res_json = json.loads(output)
        return res_json

    def get_distance_to_hull(self, alloy, off_hull_point, options = None):
        command = ' --chull --distance_to_hull=' + off_hull_point
        if options:
            command += ' ' + options
        
        output = self.aflow_command(
            command + ' --print=json --screen_only --alloy=' + alloy
        )
        #print('output: ' + output)
        res_json = json.loads(output)
        return res_json
    
    def get_stability_criterion(self, alloy, hull_point, options = None):
        command = ' --chull --stability_criterion=' + hull_point
        if options:
            command += ' ' + options
        
        output = self.aflow_command(
            command + ' --print=json --screen_only --alloy=' + alloy
        )
        #print('output: ' + output)
        res_json = json.loads(output)
        return res_json

    def get_hull_energy(self, alloy, composition, options = None):
        command = ' --chull --hull_energy=' + ','.join([ str(comp) for comp in composition ])
        if options:
            command += ' ' + options
        
        output = self.aflow_command(
            command + ' --print=json --screen_only --alloy=' + alloy
        )
        #print('output: ' + output)
        res_json = json.loads(output)
        return res_json

