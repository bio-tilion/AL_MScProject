import os
from urllib.parse import urlparse
from ftplib import FTP


def download_ftp_file(ftp_url, output_folder):
    url = urlparse(ftp_url)
    filename = os.path.basename(url.path)

    with FTP(url.hostname) as ftp:
        # login as anonymous user
        ftp.login()

        os.makedirs(output_folder, exist_ok=True)
        output_file_path = os.path.join(output_folder, filename)

        try:
            with open(output_file_path, 'wb') as f:
                ftp.retrbinary('RETR ' + url.path, f.write)
            print(f"Downloaded {filename} to {output_file_path}")
        except Exception as e:
            print(f"Failed to download {filename}. Error: {str(e)}")
