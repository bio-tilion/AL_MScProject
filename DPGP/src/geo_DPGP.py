import GEOparse
from get_ftp  import download_ftp_file


with open("DPGP/data/GEO_acc_list.txt", "r") as f:
    geo_accs = f.read().splitlines()

for geo in geo_accs:
    dir_name = f"DPGP/data/{geo}"

    # Get SOFT file from GEO Database
    gse = GEOparse.get_GEO(geo, destdir=dir_name)

    # Download each supplementary file
    for supp_file in gse.metadata["supplementary_file"]:
        download_ftp_file(supp_file, dir_name)