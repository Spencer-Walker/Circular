import json
import os
import os.path
import numpy as np
from os import path
import time 
from os.path import exists
main = os.getcwd()

b = 0.8
N = 200
a = b/N
counter = 0
for dir in [main+"/test"]:
  os.chdir(dir)
  os.system("pwd")
  root = dir
  for n in range(N+1):
    time.sleep(0.05)
    F = a*n
    folder_name = "H-" + "{:.3f}".format(F) 
    # here    
    print("mkdir " + folder_name)
    os.system("rm -r " + folder_name)
    os.system("mkdir " + folder_name)
    # here
    with open('input.json', 'r+') as f:
      data = json.load(f)
      data["EPS"]["F"] = F
      data["EPS"]["w"] = 0.25
      data["EPS"]["label"] = "H-" + "{:.3f}".format(F) 
      f.seek(0)        # <--- should reset file position to the beginning.
      json.dump(data, f, indent=4)
      f.truncate()     # remove remaining part

    os.system("cp input.json " +  folder_name +"/")
    os.system("cp submit_sim " + folder_name +"/")
    os.chdir( folder_name)
    os.system("rm *.log")
    if not exists("H-" + "{:.3f}".format(F)+".h5"):
      os.system("sbatch submit_sim")
    os.chdir(root)
print(counter)
