# Given sumstats file, generate a file that can be used by subsequent scripts to
#    generate MsCAVIAR input files.
import argparse, math, sys


def parseargs():    # handle user arguments
	parser = argparse.ArgumentParser(description="Given sumstats file, generate" +
		" a formatted file used by subsequent scripts to generate MsCAVIAR input files.")
	parser.add_argument('--infile', required=True, help = 'Summary statistics file.')
	parser.add_argument('--outfile', required=True, help = 'Output file name.')
	parser.add_argument('--chromosome', required=True,
		help = 'Name of column with chromosomes for each SNP.')
	parser.add_argument('--bp', '--pos', required=True,
		help = 'Name of column with position (BP) for each SNP.')
	parser.add_argument('--snp_id', '--snp_name', required=True,
		help = 'Name of column with names or IDs for each SNP.')
	parser.add_argument('--ref_allele', required=True,
		help = 'Name of column with reference allele for each SNP.')
	parser.add_argument('--alt_allele', required=True,
		help = 'Name of column with alternate allele for each SNP.')
	parser.add_argument('--zscore', default='NONE',
		help = 'Name of column with z-scores for each SNP. Optional.')
	parser.add_argument('--beta', '--effect_size', default='NONE',
		help = 'Name of column with beta (effect size) for each SNP. Optional.')
	parser.add_argument('--odds_ratio', default='NONE',
		help = 'Name of column with odds ratio for each SNP. Optional.')
	parser.add_argument('--se', '--stderr', '--standard_error', default='NONE',
		help = 'Name of column with standard error for each SNP. Optional.')
	parser.add_argument('--reformat_snp_ids', action='store_true',
		help = 'Standardize SNP IDs between studies by setting them to chr:pos.')
	parser.add_argument('--delimiter', default='\t',
		help = 'Column delimiter in --infile. Default is tab.')
	args = parser.parse_args()
	return args


# given the column names and the file header, extract the column numbers for
#    the fields for the output file, in order.
def get_cols_from_header(args, header):
	ordered_cols = []
	# guaranteed fields: chromosome, bp, snp_id, ref_allele, alt_allele
	ordered_cols.append(header.index(args.chromosome))
	ordered_cols.append(header.index(args.bp))
	ordered_cols.append(header.index(args.snp_id))
	ordered_cols.append(header.index(args.ref_allele))
	ordered_cols.append(header.index(args.alt_allele))

	# user can provide either zscore or beta/OR and stderr
	if args.zscore != 'NONE':
		ordered_cols.append(header.index(args.zscore))
	else:
		if args.beta != 'NONE':
			ordered_cols.append(header.index(args.beta))
		else:
			ordered_cols.append(header.index(args.odds_ratio))
		ordered_cols.append(header.index(args.se))
	return ordered_cols


# if zscore not given, calculate it from either beta or odds ratio and stderr
def calc_zscore(measure, value, stderr):
	float_stderr = float(stderr)
	#if float_stderr == 0.0:  # prevent divide by zero
	#	# we set it to the minimum nonzero value they could have specified
	#	#    based on the number of digits they gave for the stderr
	#	digits_given = len(stderr.split('.')[1]) + 10
	#	float_stderr = 10 ** (-1 * digits_given)
	if measure == 'beta':  # beta is the measure: zscore = beta / stderr
		return str(float(value) / float_stderr)
	else:  # odds ratio is the measure: zscore = log(odds_ratio) / stderr
		return str(math.log(float(value)) / float_stderr)


def main():
	args = parseargs()
	if args.zscore == 'NONE':
		if (args.beta == 'NONE' and args.odds_ratio == 'NONE') or args.se == 'NONE':
			sys.exit('Error: must specify one of the following: ' +
				'--zscore; --beta and --se; or --odds_ratio and --se.')

	with(open(args.infile, 'r')) as infile:
		with(open(args.outfile, 'w')) as outfile:
			header = infile.readline().strip().split()  # infile header
			ordered_cols = get_cols_from_header(args, header)
			outfile.write('chr pos rsid A0 A1 Zscore\n')  # outfile header
			for line in infile:
				splits = line.strip().split(args.delimiter)
				output_fields = [splits[i] for i in ordered_cols]
				if args.reformat_snp_ids:
					# attempt to standardize chromosome name to the style "chr1"
					#    and snp ID to style "chromosome:position"
					chrom, pos = output_fields[:2]
					chrom = 'chr' + chrom.strip('chromosomeCHROMOSOME')
					output_fields[0] = chrom
					output_fields[2] = chrom + ':' + pos

				if args.zscore == 'NONE':
					# if zscore not provided, calc from beta or odds ratio and stderr
					if args.beta != 'NONE':
						zscore = calc_zscore('beta', output_fields[-2], output_fields[-1])
					else:
						zscore = calc_zscore('odds_ratio', output_fields[-2], output_fields[-1])
					output_fields = output_fields[:-2]
				else:
					zscore = output_fields[-1]
					output_fields = output_fields[:-1]
				for i in output_fields:
					outfile.write(i + ' ')  # write columns according to outfile header order
				outfile.write(zscore + '\n')


if __name__ == '__main__':
	main()
#
