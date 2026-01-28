#! /usr/bin/env python3

import sys
import random  

def sliding_window_cut(input_file, output_file, window_size, step_size):
    if input_file.endswith('.gz'):
        import gzip
        opener = gzip.open
    else:
        opener = open
        
    with opener(input_file, 'rt') as f_in, open(output_file, 'w') as f_out:
        a = 0
        header, seq, sep, qual = '', '', '', ''
        for line in f_in:
            line = line.rstrip('\n')  
            if a % 4 == 0:  
                header = line
                if not header.startswith('@'):
                    header = '@' + header.lstrip('@')  
            elif a % 4 == 1:  
                seq = line
            elif a % 4 == 2:  
                sep = line  
            elif a % 4 == 3:  
                qual = line
                
                if len(seq) < window_size:
                    a += 1
                    continue  
                
                if len(seq) != len(qual):
                    a += 1
                    continue  
                
                read_id = header.split()[0].lstrip('@')
                seq_len = len(seq)
                current_direction = random.choice(['left', 'right']) ## 随机从左右切
                if current_direction == 'left':
                    win_starts = range(0, seq_len - window_size + 1, step_size)
                elif current_direction == 'right':
                    win_starts = range(seq_len - window_size, -1, -step_size)

                
                for start in win_starts:
                    end = start + window_size
                    sub_seq = seq[start:end]
                    sub_qual = qual[start:end]                                        
                    new_header = f"@{read_id}:window_{start}-{end}"
                    f_out.write(f"{new_header}\n{sub_seq}\n{sep}\n{sub_qual}\n")
            a += 1

if __name__ == "__main__":
    input_file, output_file = sys.argv[1], sys.argv[2]
    window_size, step_size = int(sys.argv[3]), int(sys.argv[4])
    sliding_window_cut(input_file, output_file, window_size, step_size)