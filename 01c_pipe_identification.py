import subprocess
import argparse
import os
import os.path
import functions_identification_pipe


parser = argparse.ArgumentParser(description = "Pipeline to identify intratumoral bacteria.")
parser.add_argument("--path_R1", help = "Path to read 1.", required = True)
parser.add_argument("--path_R2", help = "Path to read 2.", required = True)
parser.add_argument("--output", help = "Directory of output, at the end there should be no /.", required = True)
parser.add_argument("--base", help = "Base name of files.", required = True)

args = parser.parse_args()
path_R1 = args.path_R1
path_R2 = args.path_R2
output = args.output
base = args.base


print("[+] Starting Kraken.")
functions_identification_pipe.kraken(path_R1, path_R2, output, base)


print("[+] Starting Bracken.")
functions_identification_pipe.bracken(output, base)