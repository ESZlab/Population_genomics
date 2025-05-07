import pandas as pd
import argparse
import os

def convert_gen_to_vcf(input_file, output_file, nchr, lenchr):
    gen_file = pd.read_csv(input_file, sep='\t', index_col=False)
    gen_file.columns = [
        'CHROM', 'POS', 'REF', 'ALT',
        'NCJ1', 'NCJ2', 'NCJ3', 'NCJ4', 'NCJ5', 'NCJ6',
        'USCA1', 'USCA2', 'USCA3', 'USCA4', 'USCA5', 'USCA6',
        'AZM1', 'AZM2', 'AZM3', 'AZM4', 'AZM5', 'AZM6', 'AZM7', 'AZM8', 'AZM9',
        'AZJ1', 'AZJ2', 'AZJ3', 'AZJ4', 'AZJ5', 'AZJ6', 'AZJ7', 'AZJ8', 'AZJ9',
        'ITTC1', 'ITTC2', 'ITTC3', 'ITTC4', 'ITTC5', 'ITTC6'
    ]

    # Insert required VCF columns
    gen_file.insert(2, 'ID', '.')
    gen_file.insert(5, 'QUAL', '.')
    gen_file.insert(6, 'FILTER', '.')
    gen_file.insert(7, 'INFO', '.')
    gen_file.insert(8, 'FORMAT', 'GT')

    with open(output_file, 'w') as vcf_file:
        # VCF headers
        vcf_file.write("##fileformat=VCFv4.2\n")
        vcf_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

        for i in range(nchr):
            vcf_file.write(f"##contig=<ID={i+1},length={lenchr},assembly=SIMULATION.fa>\n")

        # Write column headers
        vcf_file.write('#' + '\t'.join(gen_file.columns) + '\n')

        # Replace genotype codes
        replace_values = {0: '0/0', 1: '0/1', 2: '1/1'}
        for column in gen_file.columns:
            if column not in ['CHROM', 'POS', 'ID', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'REF', 'ALT']:
                gen_file[column] = gen_file[column].replace(replace_values)

        # Write data rows
        for _, row in gen_file.iterrows():
            row_str = '\t'.join(map(str, row.values))
            vcf_file.write(row_str + '\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert GEN file to VCF.')
    parser.add_argument('--input', required=True, help='Path to input .gen file')
    parser.add_argument('--output', required=True, help='Path to output .vcf file')
    parser.add_argument('--nchr', type=int, default=10000000, help='Number of chromosomes to list in contigs')
    parser.add_argument('--lenchr', type=int, default=100, help='Length of each chromosome')

    args = parser.parse_args()
    convert_gen_to_vcf(args.input, args.output, args.nchr, args.lenchr)
