from Bio import SeqIO, SeqUtils
import drep



def get_composite_gc(lens, GCs):
    return reduce(
        lambda x,y: x+y, 
        map(
            lambda x,y: x*y, 
            lens, 
            GCs
        ), 
        0
    ) / sum(lens)

def get_bin_stats(bin_fp):
    stats = {}
    lens = []
    GCs = []
    for seq in SeqIO.parse(bin_fp, 'fasta'):
        lens.append(len(seq))
        GCs.append(SeqUtils.GC(seq.seq))
    
    stats['GC'] = get_composite_gc(lens, GCs)
    stats['length'] = sum(lens)
    stats['N50'] = drep.d_filter.calc_n50(bin_fp)

    return stats








