#!/usr/bin/env python3
import argparse
import os


def main():
    parser = argparse.ArgumentParser(description="Create overall index file that mirrors the Alliance reports S3 "
                                                 "bucket")
    parser.add_argument("-i", "--input-dir", metavar="input_dir", dest="input_dir", type=str,
                        default="./", help="working directory where input files are located")
    args = parser.parse_args()

    print("<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
    print("<ListBucketResult xmlns=\"http://s3.amazonaws.com/doc/2006-03-01/\">")
    print("\t<Name>wb-db-reports</Name>")
    print("\t<Prefix></Prefix>")
    print("\t<Marker></Marker>")
    print("\t<MaxKeys>1000</MaxKeys>")
    print("\t<IsTruncated>false</IsTruncated>")

    for release_version in os.listdir(args.input_dir):
        for release_type in os.listdir(release_version):
            for release_date in os.listdir(release_type):
                for json_file_path in os.listdir(release_date):
                    if json_file_path.endswith(".json"):
                        print("\t<Contents>")
                        print("\t\t<Key>gene-descriptions/" + release_version + "/" + release_type + "/" +
                              release_date + "/" + json_file_path + "</Key>")
                        print("\t</Contents>")

    print("</ListBucketResult>")


if __name__ == '__main__':
    main()
