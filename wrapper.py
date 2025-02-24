import os
import subprocess
from Bio import SeqIO
#define the directory name dynamicaly
dir_name = "/home/2025/kboreiri/PipelineProject_Kimia_Boreiri/output"
#move into the directory
os.chdir(dir_name)

#varify the current working directory
print("Current Directory:", os.getcwd())



