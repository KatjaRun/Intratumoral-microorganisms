import subprocess
import os.path
import time

# Functions for the identification pipeline

# Kraken2
def kraken(path_R1, path_R2, output, base):
    output_name = base + ".kraken.txt"
    output_path = output + "/" + output_name
    report_name = base + ".kreport2"
    report_path = output + "/" + report_name
    print("File path 1: ", path_R1)
    print("File path 2: ", path_R2)
    print("Output path: ", output_path)
    print("Report path: ", report_path)
    subprocess.call([
        '/home/rungger/.conda/envs/Kraken2/bin/kraken2',
        '--db', '/home/rungger/myScratch/Databases/NT_database/nt-fast',
        '--threads', '20',
        '--output',f'{output_path}',
        '--paired', f'{path_R1}',
        f'{path_R2}',
        '--report', f'{report_path}'])

# Bracken
def bracken(output, base):
    report_name = base + ".kreport2"
    report_path = output + "/" + report_name
    output_name = base + ".bracken"
    output_path = output + "/" + output_name
    print("Report path:", report_path)
    print("Output path: ", output_path)
    subprocess.call([
        "/home/rungger/.conda/envs/Bracken/Bracken-master/bracken",
        "-d", "/home/rungger/myScratch/Databases/NT_database/nt-fast",
        "-i", f"{report_path}",
        "-o", f"{output_path}",
        "-l", "S"])
    
# KrakenBiom for all files
def krakenbiom(output):
    output_path = os.listdir(output)
    bracken_reports = []
    output_biom = output + "/Final.biom"

    for files in output_path:
        if files.endswith("_bracken_genuses.kreport2"):
            bracken_reports.append(os.path.join(output, files))
    print("[+] Creating biom file:")
    print(output_biom)
    
    subprocess.call([
        "kraken-biom",
        *bracken_reports,
        "-o", f"{output_biom}",
        "--fmt", "json"
        ])
