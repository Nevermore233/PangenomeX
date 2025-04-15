import os
import time
import warnings
from tqdm import trange
import pandas as pd
import argparse
from utils import *
warnings.filterwarnings('ignore')

parser = argparse.ArgumentParser()
parser.add_argument('-config', type=str, help="Path to the '.config' file ", required=True)
args = parser.parse_args()


def init_dict(seq):
    initialized_dict = {}
    for key, value in seq.items():
        initialized_dict[key] = 'R' * len(value)
    return initialized_dict


def fill_seq_from_bam(sequence, bam_file):
    alignments = read_bam_file(bam_file)

    for alignment in alignments:
        rname = alignment['RNAME']
        pos = alignment['POS']
        seq = alignment['SEQ']

        if rname in sequence:
            start_idx = pos - 1
            end_idx = start_idx + len(seq)

            if not isinstance(sequence[rname], np.ndarray):
                sequence[rname] = np.array(list(sequence[rname]))

            sequence[rname][start_idx:end_idx] = list(seq)

    for rname in sequence:
        sequence[rname] = ''.join(sequence[rname])

    return sequence


def replace_with_ref(sequence, filled_sequence):
    replaced_sequence = filled_sequence.copy()
    for key, value in filled_sequence.items():
        original_value = sequence[key]
        replaced_value = ''
        for i in range(len(value)):
            if value[i] == 'R':
                replaced_value += original_value[i]
            else:
                replaced_value += value[i]
        replaced_sequence[key] = replaced_value
    return replaced_sequence


def write_seq_to_file(sequence, filename):
    with open(filename, 'w') as file:
        file.write(f'>{file}\n{sequence}\n')


def main():
    config_file = args.config
    paths = read_config(config_file)

    # Input
    # Trained files
    train_file_list = paths['train_file_list']
    train_df = pd.read_csv(train_file_list, index_col=0)

    # Tested files
    test_file_list = paths['test_file_list']
    test_df = pd.read_csv(test_file_list, index_col=0)

    ref_file = paths['ref_file']

    sequence = read_fasta_file(ref_file)
    initialized_sequence = init_dict(sequence)

    for i in trange(train_df.shape[0]):
        time.sleep(0.01)
        mapping = f"sample_{i}"
        train_filename = train_df.loc[train_df['mapping'] == mapping]['file_name'].values[0]
        print(f'process {train_filename} ..................')

        # output: seq.fas
        file_name = os.path.splitext(train_filename)[0] + '.fas'
        filled_sequence = fill_seq_from_bam(initialized_sequence, train_filename)
        replaced_sequence = replace_with_ref(sequence, filled_sequence)
        concatenated_sequence = ''.join(replaced_sequence.values())
        write_seq_to_file(concatenated_sequence, file_name)

    for i in trange(test_df.shape[0]):
        time.sleep(0.01)
        mapping = f"sample_{i}"
        test_filename = test_df.loc[test_df['mapping'] == mapping]['file_name'].values[0]
        print(f'process {test_filename} ..................')

        # output: seq.fas
        file_name = os.path.splitext(test_filename)[0] + '.fas'
        filled_sequence = fill_seq_from_bam(initialized_sequence, test_filename)
        replaced_sequence = replace_with_ref(sequence, filled_sequence)
        concatenated_sequence = ''.join(replaced_sequence.values())
        write_seq_to_file(concatenated_sequence, file_name)


if __name__ == '__main__':
    st = time.time()
    main()
    et = time.time()
    rt = et - st
    print(f"Finish! runtime: {rt}sec")
