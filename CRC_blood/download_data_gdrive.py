# Downloads the data directly from google drive
import gdown
import os
import argparse

def download_data(dir_url):
    # download the whole directory from dir_url
    gdown.download_folder(dir_url,  quiet=False, use_cookies=False)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--link', type=str)
    args = parser.parse_args()
    download_data(args.link)

if __name__ == "__main__":
    main()
