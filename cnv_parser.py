import sys

def parse_line_maf2cnv(line, idx):
    fields = line.strip().split("\t")
    result = 'c' + str(idx) + '\t'
    return result

def make_line_cnv(idx, a, d, ssms=[], physical_cnvs=[]):
    result = 'c' + str(idx) + '\t' + str(a) + '\t' + str(d) + '\t'
    for ssm in ssms:
        result += ssm + ';'
    return result

def parse_line_cnv(fields, idx, ssm_dict):
    result = 'c' + str(idx) + '\t'
    a = fields[-1]
    d = fields[-2]
    ssms = ssm_dict[fields[0]]
    result += str(a) + '\t' + str(d) + '\t'
    for ssm in ssms:
        result += ssm + ';'
    result += '\t'
    return result



def parse_cnv(cnv_path, ssm_path, outpath = './cnv_data.txt'):
    ssm_in = open(ssm_path, "r")
    
    # skip header
    ssm_in.readline()
    ssm_dict = {}
    for line in ssm_in:
        fields = line.strip().split("\t")[1]
        gene_name = fields[1]
        ssm_id = fields[0]
        if gene_name not in ssm_dict:
            ssm_dict[gene_name] = []
        ssm_dict[gene_name].append(ssm_id)

    # start processing cnv data
    cnv_in = open(cnv_path, "r")
    cnv_header = line.strip().split("\t")
    line = cnv_in.readline()
    fout = open(outpath, "w")
    header = 'cnv\ta\td\tssms\tphysical_cnvs\n'
    fout.write(header)
    idx = 0
    for line in cnv_in:
        fields = line.strip().split("\t")
        if len(fields) == 8 and fields[-1] != fields[-2]:
            print(fields)
            if fields[0] in ssm_dict:
                fout.write(parse_line_cnv(fields, idx, ssm_dict))
                idx += 1

if __name__ == "__main__":
    args = sys.argv[1:]
    sample_barcode = args[0]
    ssm_path = 'result/'+sample_barcode+'.txt'
    cnv_path = 'data/'+sample_barcode+'/'+args[1]
    outpath = 'result/'+args[2]
    parse_cnv(cnv_path, ssm_path, outpath)