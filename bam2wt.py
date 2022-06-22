# read bed with seq from stdout
# use wavelets to localize CpG sites in both time and frequency (i.e., the order and density) domains

import pywt # mamba install pywavelets
import numpy as np
import sys, argparse

def truncate() -> int:
    """Read arguments from command line

    Returns:
        int: the number of bases to keep, the wavelet level
    """
    parser = argparse.ArgumentParser()

    parser.add_argument('bases', type=int, help='- how many bases to keep in 5\' to 3\' direction')
    parser.add_argument('level', type=int, help='- what is the wavelet decomposition level, i.e., the level 2 for 1000 bases will report 250 points')

    args = parser.parse_args()

    return args.bases, args.level

def get_mother_WT(data: np.ndarray, restored = False, level = 1) -> np.ndarray:
    """Wavelet decomposition

    Args:
        data (np.ndarray): array of time series values
        restored (bool, optional): return the "smoothed" = True or "compressed" = False. Defaults to False.
        level (int, optional): Wavelets decomposition level ("compression"). Defaults to 1.

    Returns:
        np.ndarray: result of wavelet decomposition
    """
    wavelet_family = 'haar'
    mode = 'smooth'
    coeff_all = pywt.wavedec(
        data, 
        wavelet = wavelet_family, 
        mode=mode, 
        level=level)
    cA, cD = coeff_all[0], coeff_all[1:]
    
    if restored:
        n = data.shape[0]
        omp0 = pywt.upcoef('a', cA, wavelet_family, level=level)[:n]
        return omp0 # restored mother wavelet
    else:
        return cA # compressed mother wavelet

def seq2boolean(seq: str) -> np.ndarray:
    """Convert a sequence to its boolean representation
    replace CG (case sensitive) by "1", other bases by "0"

    Args:
        seq (str): DNA sequence

    Returns:
        np.ndarray: result
    """
    seq = seq.replace('CG', '10').lower() # replace CG

    for _ in 'atgcn': # all other bases
        seq = seq.replace(_, '0')

    return np.asarray(list(seq)).astype(int)

def read_stdin() -> str:
    for line in sys.stdin:
        yield line.rstrip()
        
# main

bases, level = truncate()

read_stdin_generator = read_stdin()       

line = next(read_stdin_generator, False)

while line:

    # get tab separated string like this:
    # chr19	54579701	54589400	UID-1	0	-	ctc...

    chrom, start, end, name, quality, strand, seq = line.split('\t')

    if len(seq) < bases:
        seq = seq.rjust(bases, '0')

    seq = seq[:bases]

    boolean_seq = seq2boolean(seq)
    
    wavelet = get_mother_WT(boolean_seq, False, level)

    # return the region details and the seq wavelet
    print(
        f'{name}::{chrom}:{start}-{end}({strand})',
        *wavelet,
        sep='\t',
        file=sys.stdout)

    line = next(read_stdin_generator, False)

sys.exit(0)
