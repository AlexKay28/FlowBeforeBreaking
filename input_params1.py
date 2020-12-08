import configparser
import json

config = configparser.ConfigParser()
config.read('config1.ini')

# data for problem object
Dout = float(config['tube params']['Dout'])
t = float(config['tube params']['t'])


cycles_info ={'НУЭ': json.loads(config['cycles info']['NUE']),
              'ГИ1': json.loads(config['cycles info']['GI1']),
              'ГИ2': json.loads(config['cycles info']['GI2'])}

# coeff for Paris' formula
C = float(config['coeffs']['C'])
m = float(config['coeffs']['m'])

# start point data for deffect
a_0 = 0.2*t
c_0 = 0.5*t
