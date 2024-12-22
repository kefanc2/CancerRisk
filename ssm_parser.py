import sys

def parse_line_maf2ssm(line, idx, gene_idx, ref_idx=40, alt_idx=39):
    fields = line.strip().split("\t")
    ref = fields[ref_idx]
    alt = fields[alt_idx]
    gene = fields[gene_idx]
    result = 's' + str(idx) + '\t' + gene + '\t' + str(ref) + '\t' + str(alt) + '\t'
    result += '0.999\t0.499\n'
    return result

def parse_ssm(inpath, outpath):
    mafin = open(inpath, "r")
    header = 'id\tgene\ta\td\tmu_r\tmu_v\n'
    fout = open(outpath, "w")
    fout.write(header)
    # skip header
    for line in mafin:
        if line.startswith("#"):
            continue
        else:
            break
    maf_header = line.strip().split("\t")
    gene_idx = maf_header.index("Hugo_Symbol")
    idx = 0
    for line in mafin:
        fout.write(parse_line_maf2ssm(line, idx, gene_idx))
        idx += 1

if __name__ == "__main__":
    """
    usage: python ssm_parser.py input output
    turns a maf file into a ssm file
    """
    args = sys.argv[1:]
    inpath = args[0]
    outpath = args[1]
    parse_ssm(inpath, outpath)