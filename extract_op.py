import os
import pysam
import argparse
import sys

def extrack_denovo(input_vcf, min_support=1):


    output_vcf = open(os.path.join(input_vcf.replace(".vcf", ".denovo_s{}.vcf".format(min_support))), "w")

    input_vcf = pysam.VariantFile(input_vcf)

    output_vcf.write(str(input_vcf.header))

    for record in input_vcf:

        if "POLY" in str(record):
            continue

        if record.info["SUPPORT"] < min_support:
            continue

        if "_NewCOMP_NewCOMP" in str(record):
            output_vcf.write(str(record))


    input_vcf.close()
    output_vcf.close()


def extrack_somatic(input_vcf, min_support):

    output_vcf = open(os.path.join(input_vcf.replace(".vcf", ".somatic_s{}.vcf".format(min_support))), "w")

    input_vcf = pysam.VariantFile(input_vcf)

    output_vcf.write(str(input_vcf.header))

    for record in input_vcf:

        if "POLY" in str(record):
            continue

        if record.info["SUPPORT"] < min_support:
            continue

        if "_NewCOMP" in str(record):
            output_vcf.write(str(record))

    input_vcf.close()
    output_vcf.close()


def extrack_support(input_vcf, min_support):

    output_vcf = open(os.path.join(input_vcf.replace(".vcf", ".extract_s{}.vcf".format(min_support))), "w")

    input_vcf = pysam.VariantFile(input_vcf)

    output_vcf.write(str(input_vcf.header))

    for record in input_vcf:

        if record.info["SUPPORT"] < min_support:
            continue

        output_vcf.write(str(record))

    input_vcf.close()
    output_vcf.close()


def parse_arguments(arguments=sys.argv[1:]):
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description="Extract candidate denovo/somatic SVs from SVision-pro's output VCF")

    inout_params = parser.add_argument_group("Input/Output parameters")
    inout_params.add_argument('--input_vcf', dest="input_vcf", type=os.path.abspath, required=True, help='Absolute path to the input VCF')
    inout_params.add_argument('--min_supp', dest="min_supp", type=int, default=1, help='Minimum support read number required for SV calling (default: %(default)s)')
    inout_params.add_argument('--extract', dest='extract', required=True, choices=['denovo', 'somatic', "support"],  help='Extract calls by denovo, somatic or supports')

    return parser.parse_args(arguments)


if __name__ == '__main__':

    options = parse_arguments()

    extract_mode = options.extract

    if extract_mode == "denovo":
        extrack_denovo(options.input_vcf, options.min_supp)

    elif extract_mode == "somatic":
        extrack_somatic(options.input_vcf, options.min_supp)

    elif extract_mode == "support":
        extrack_support(options.input_vcf, options.min_supp)

    else:

        print("No such extract mode")
        exit(-1)

