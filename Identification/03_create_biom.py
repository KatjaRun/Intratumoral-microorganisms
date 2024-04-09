import os
import subprocess
import argparse

parser = argparse.ArgumentParser(description = "Pipeline to identify intratumoral bacteria.")
parser.add_argument("--input", help = "Directory of input", required = True)
parser.add_argument("--n_samples", help = "Number of samples being worked on.",  required = True)

args = parser.parse_args()
input = args.input
n_samples = int(args.n_samples)

# Function for creating a BIOM file out of bracken kreport output
def krakenbiom(input):
    input_path = os.listdir(input)
    bracken_reports = []
    output_biom = input + "/Final_nt.biom"

    for files in input_path:
        if files.endswith("_bracken_species.kreport2"):
            path = input + "/" + files
            bracken_reports.append(path)
    print("[+] Creating biom file:")
    print(output_biom)
    print("From ", len(bracken_reports), " files.")
    
    if len(bracken_reports) == n_samples:
        print("Yay, all samples are getting processed.")
    elif len(bracken_reports) < n_samples:
        print("Not all samples are getting processed.")
    elif len(bracken_reports) > n_samples:
        print("More samples than expected, wtf.")
    
    subprocess.run([
        "kraken-biom",
        *bracken_reports,
        "-o", f"{output_biom}",
        "--fmt", "json"
        ])
    
    print("[+] Done!")
    
    
# Running Function
krakenbiom(input)