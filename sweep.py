import json
import os
import os.path
import time 
import numpy as np
from os import path

main = os.getcwd()

b = 1.0
a = 0.001
N = 1000

counter = 0
for dir in [main+"/test"]:
  os.chdir(dir)
  os.system("pwd")
  root = dir
  for n in range(N+1):
    time.sleep(1)
    F = a*n
    folder_name = "H-" + "{:.3f}".format(F) 
    print("mkdir " + folder_name)
    os.system("rm -r " + folder_name)
    os.system("mkdir " + folder_name)

    with open('input.json', 'r+') as f:
      data = json.load(f)
      data["EPS"]["F"] = F
      data["EPS"]["label"] = "H-" + "{:.3f}".format(F) 
      f.seek(0)        # <--- should reset file position to the beginning.
      json.dump(data, f, indent=4)
      f.truncate()     # remove remaining part

    os.system("cp input.json " +  folder_name +"/")
    os.system("cp submit_sim " + folder_name +"/")
    os.chdir( folder_name)
    os.system("rm *.log")
    os.system("sbatch submit_sim")
    os.chdir(root)
print(counter)
