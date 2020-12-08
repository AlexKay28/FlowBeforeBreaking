import configparser
import json

config = configparser.ConfigParser()
config.read('config2.ini')

# data for problem object
Dout = float(config['tube params']['Dout'])
t = float(config['tube params']['t'])


p = float(config['info']['p'])
steel_type = str(config['info']['steel_type'])
Nz = float(config['info']['Nz'])
Mx = float(config['info']['Mx'])
My = float(config['info']['My'])


# coeff for Paris' formula
C = float(config['steel']['C'])
m = float(config['steel']['m'])
Rp02_min = float(config['steel']['Rp02_min'])
Rp02_max = float(config['steel']['Rp02_max'])
Rm_min = float(config['steel']['Rm_min'])
Rm_max = float(config['steel']['Rm_max'])
E_module = float(config['steel']['E_module'])
mu = float(config['steel']['mu'])
T = float(config['steel']['T'])

# start point data for deffect
a_0 = 0.2*t
c_0 = 0.5*t
