#! /usr/bin/env python3

## the script is for cutting reads in BAM file
from pysam import AlignmentFile, view, sort
from random import randint
import sys, os, re
from collections import OrderedDict

def arg(argv):

    import argparse
    parser = argparse.ArgumentParser(description="Analysis using Python", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-i', '--input', help="Alignment File with samtools indexed", type=str, required=True)
    parser.add_argument('--prefix', help="Prefix for output file", type=str)
    parser.add_argument('--length', help="Lengths for cut reads. If multiple lengths were to be supplied, then separate them with ','", type=str)
    parser.add_argument('--direct', help="Direction for cutting reads [str, {'left' (default), 'random', 'right'}]", type=str, default='left')
    parser.add_argument('--paired', help="Single-end or paired-end [Bool]. Note: For paired-end mode, the script needs further modification", action='store_true', default=False)
   
    args = parser.parse_args() 
    sys.stderr.write('args: '+str(args)+'\n')
    main(args.input, args.prefix, length_list=[int(x) for x in args.length.split(',')], direction=args.direct, paired=args.paired)

def main(inFile, out_prefix, length_list=None, direction='left', paired=False):

    bam_in = AlignmentFile(inFile)
    header = view('-H', inFile)
    for length in length_list:
        # sys.stderr.write('CURRENT: LENGTH='+str(length)+'\n')
        out_file = open(out_prefix+'.LENGTH_'+str(length)+'.sam', 'w')
        out_file.write(header); query_dict = OrderedDict()
        for bam in bam_in.fetch():
            query_id = bam.query_name
            #sys.stderr.write('\nquery_id: '+query_id+'\n')
            flag_info = bam.flag
            #sys.stderr.write('flag_info: '+str(flag_info)+'\n')
            ref_id = bam_in.get_reference_name(bam.reference_id)
            if direction == 'random': 
                if len(bam.query_sequence)-length < 0: continue
                rand_seed = randint(0, len(bam.query_sequence)-length)
            else: rand_seed = 0
            #sys.stderr.write('rand_seed: %s\nORIGINAL startPos: %s\n' % (rand_seed, bam.reference_start+1))
            CIGAR_info, ref_pos = CIGAR_startPos_modification(bam.cigartuples, bam.reference_start+1, bam.query_sequence, length=length, direction=direction, rand_seed=rand_seed)
            map_q = bam.mapping_quality
            #sys.stderr.write('FINAL CIGAR_info: %s; ref_pos: %s\n' % (CIGAR_info, ref_pos))
            if paired:
                mate_id = bam.next_reference_name
                mate_pos = bam.next_reference_start
                t_length = bam.template_length
            else: 
                mate_id = "*"; mate_pos = 0; t_length = 0
            seq = seq_modification(bam.query_sequence, length=length, direction=direction, rand_seed=rand_seed)
            quality = qual_modification(''.join([chr(x+33) for x in bam.query_qualities]), length=length, direction=direction, rand_seed=rand_seed)
            tags = '\t'.join(['%s:%s:%s' % (x[0],'i',x[1]) for x in bam.get_tags(with_value_type=True)])
            if not paired:
                out_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (query_id, flag_info, ref_id, ref_pos, map_q, CIGAR_info, mate_id, mate_pos, t_length, seq, quality, tags)) 
            else:
                query_dict.setdefault(query_id, []).append([flag_info, ref_id, ref_pos, map_q, CIGAR_info, mate_id, mate_pos, t_length, seq, quality, tags])

        if paired:
            #sys.stderr.write('ALL information has loaded.\n')
            #sys.stderr.write('length: '+str(length)+'\n')
            for query_id in query_dict:
                if len(query_dict[query_id]) > 2:
                    #sys.stderr.write('Current query_id: '+str(query_id)+'\n')
                    paired_1 = [_q for _q in query_dict[query_id] if bin(_q[0])[2:][::-1][6] == '1'] # first in pair
                    paired_2 = [_q for _q in query_dict[query_id] if bin(_q[0])[2:][::-1][6] == '0'] # second in pair
                    if len(paired_1) == 0 and len(paired_2) != 0:
                        for _p2 in paired_2:
                            _p2[6] += 1
                            out_file.write('%s\t%s\n' % (query_id, '\t'.join([str(x) for x in _p2])))
                        continue
                    elif len(paired_2) == 0 and len(paired_1) != 0:
                        for _p1 in paired_1:
                            _p1[6] += 1
                            out_file.write('%s\t%s\n' % (query_id, '\t'.join([str(x) for x in _p1])))
                        continue
                    for _p1 in paired_1:
                        _p2_list = [_p for _p in paired_2 if _p[1] == _p1[5]]
                        _p2 = sorted([(abs(_p[2]-_p1[6]), _p) for _p in _p2_list], key=lambda x:x[0])[0][1]
                        _p1[6] = _p2[2]; _p2[6] = _p1[2]
                        _p1_right = _p1[2]+CIGAR_length(_p1[4])-1; _p2_right = _p2[2]+CIGAR_length(_p2[4])-1
                        sites = [_p1[2], _p2[2], _p1_right, _p2_right]
                        if _p1[2] < _p2[2]: _p1[7] = max(sites) - min(sites) + 1
                        else: _p1[7] = (-1) * (max(sites) - min(sites) + 1)
                        if _p2[2] < _p1[2]: _p2[7] = max(sites) - min(sites) + 1
                        else: _p2[7] = (-1)*(max(sites) - min(sites) + 1)
                        out_file.write('%s\t%s\n' % (query_id, '\t'.join([str(x) for x in _p1])))
                        out_file.write('%s\t%s\n' % (query_id, '\t'.join([str(x) for x in _p2])))
                    #sys.stderr.write('QUERY_ID %s has more than two mapping results.\n' % (query_id))
                    continue
                elif len(query_dict[query_id]) == 1:
                    query_dict[query_id][0][6] += 1
                    out_file.write('%s\t%s\n' % (query_id, '\t'.join([str(x) for x in query_dict[query_id][0]])))
                    #sys.stderr.write('QUERY_ID %s has no paired reads\n' % (query_id))
                    continue
                query_dict[query_id][0][6] = query_dict[query_id][1][2]
                _p1_right = query_dict[query_id][0][2]+CIGAR_length(query_dict[query_id][0][4])-1
                _p2_right = query_dict[query_id][1][2]+CIGAR_length(query_dict[query_id][1][4])-1
                sites = [query_dict[query_id][0][2], query_dict[query_id][1][2], _p1_right, _p2_right]
                if query_dict[query_id][0][2] < query_dict[query_id][1][2]: query_dict[query_id][0][7] = max(sites) - min(sites) + 1
                else:  query_dict[query_id][0][7] = (-1)*(max(sites) - min(sites) + 1)
                query_dict[query_id][1][6] = query_dict[query_id][0][2]
                if query_dict[query_id][1][2] < query_dict[query_id][0][2]: query_dict[query_id][1][7] = max(sites) - min(sites) + 1
                else: query_dict[query_id][1][7] = (-1)*(max(sites) - min(sites) + 1)
                for _q in query_dict[query_id]:
                    out_file.write('%s\t%s\n' % (query_id, '\t'.join([str(x) for x in _q])))
    out_file.close()

def CIGAR_length(CIGAR):

    CIGAR_tuple = [(int(x[:-1]), x[-1]) for x in re.findall(r'\d*[MNIDSH]', CIGAR)]
    LENGTH = sum([X[0] for X in CIGAR_tuple if X[1] in ['M', 'N', 'D']])
    return(LENGTH)

def CIGAR_startPos_modification(CIGAR_info, startPos, seq, length=0, direction='left', rand_seed=0):
    # [(0, 30), ...], (operation, length), operation order: match, ins, del, ref_skip, soft_clip, hard_clip, pad, equal, diff, back

    total_length = 0
    operation_dict = {0:'M', 1:'I', 2:'D', 3:'N', 4:'S', 5:'H', 6:'P', 7:'=', 8:'X', 9:'B', 10:'NM'}
    out_CIGAR = ''

    if direction == 'left':
        for _CIGAR in CIGAR_info:
            _opr = _CIGAR[0]; _len = _CIGAR[1]
            if _opr in [0,1,4]: # match, ins, soft_clip
                if total_length + _len >= length:
                    tmp_len = length-total_length
                    out_CIGAR += str(tmp_len)+operation_dict[_opr]
                    break
                else: 
                    total_length += _len
                    out_CIGAR += str(_len)+operation_dict[_opr]
            elif _opr in [2,3]: # del, ref_skip
                out_CIGAR += str(_len)+operation_dict[_opr]
            elif _opr == 5: continue # hard_clip 
            else: return('UNKNOWN_CIGAR', startPos)
        return(out_CIGAR, startPos)

    elif direction == 'right':
        for _CIGAR in CIGAR_info[::-1]:
            _opr = _CIGAR[0]; _len = _CIGAR[1]
            if _opr in [0,1,4]: # match, ins, soft_clip
                if total_length + _len >= length:
                    tmp_len = length-total_length
                    out_CIGAR = str(tmp_len)+operation_dict[_opr] + out_CIGAR
                    break
                else: 
                    total_length += _len
                    out_CIGAR = str(_len)+operation_dict[_opr] + out_CIGAR
            elif _opr in [2,3,5]: # del, ref_skip, hard_clip
                out_CIGAR = str(_len)+operation_dict[_opr] + out_CIGAR
            else: return('UNKNOWN_CIGAR', startPos)
        out_CIGAR_list = [[x[:-1], x[-1]] for x in re.findall(r'\d*[MNIDSH]', out_CIGAR)]
        if out_CIGAR_list[0][-1] == 'I': out_CIGAR_list[0][-1] = 'S'
        out_CIGAR = ''.join([''.join(x) for x in out_CIGAR_list])

        tmp_len = 0; s_len = 0  # to cound the startPos
        # sys.stderr.write('startPos: '+str(startPos)+'\n')
        for _CIGAR in CIGAR_info:
            _opr = _CIGAR[0]; _len = _CIGAR[1]
            #sys.stderr.write('_opr: %s; _len: %s\n' % (_opr, _len))
            #sys.stderr.write('tmp_len: %s; startPos: %s; s_len: %s\n' % (tmp_len, startPos, s_len))
            if _opr in [0]:
                # sys.stderr.write('len_seq: %s; length: %s; tmp_len: %s\n' % (len(seq), length, tmp_len))
                if tmp_len + _len < len(seq) - length+1:
                    startPos += _len; tmp_len += _len
                else:
                    startPos += (len(seq)-length)-tmp_len#; tmp_len += _len
                    break
            elif _opr in [4]: 
                if tmp_len+_len < len(seq)-length+1: tmp_len += _len; s_len += _len
                else: break
            elif _opr in [1]: 
                if tmp_len + _len < len(seq) - length+1: tmp_len += _len
                else: 
                    break
            elif _opr in [2,3,5]: startPos += _len
        #sys.stderr.write('startPos: '+str(startPos)+'\n')
        return(out_CIGAR, startPos)

    else: # direction == 'random'
        moved_len = 0; total_length = 0; start = False
        for _CIGAR in CIGAR_info:
            _opr = _CIGAR[0]; _len = _CIGAR[1]
            #sys.stderr.write('original CIGAR: %s\nmoved_len: %s; startPos: %s; total_length: %s\n' % (str(_len)+operation_dict[_opr], moved_len, startPos, total_length))
            if _opr in [0]: # match
                if moved_len + _len <= rand_seed and not start:
                    moved_len += _len; startPos += _len
                else:
                    if not start: # the site where start loading CIGAR information
                        if _len - (rand_seed - moved_len) > length:
                            out_CIGAR += str(length) + operation_dict[_opr]
                            startPos += rand_seed - moved_len
                            break
                        else:
                            out_CIGAR += str(_len - (rand_seed - moved_len)) + operation_dict[_opr]
                            total_length += _len - (rand_seed - moved_len)
                            startPos += rand_seed - moved_len
                            start = True; continue
                        
                    else: # the CIGAR region is for continuing the loading
                        if total_length + _len <= length:
                            out_CIGAR += str(_len)+operation_dict[_opr]
                            total_length += _len
                        else:
                            out_CIGAR += str(length-total_length)+operation_dict[_opr]
                            break
            elif _opr in [4]: # soft_clip
                if moved_len + _len <= rand_seed and not start:
                    moved_len += _len
                else:
                    if not start: # the site where start loading CIGAR information
                        if _len - (rand_seed - moved_len) > length:
                            out_CIGAR += str(length) + operation_dict[_opr]
                            break
                        else:
                            out_CIGAR += str(_len - (rand_seed - moved_len)) + operation_dict[_opr]
                            total_length += _len - (rand_seed - moved_len)
                        # if _len - (rand_seed - moved_len) <= length:
                        #     out_CIGAR += str(_len - (rand_seed - moved_len)) + operation_dict[_opr]
                        #     total_length += _len - (rand_seed - moved_len)
                        # else:
                        #     out_CIGAR += str(length - (rand_seed - moved_len)) + operation_dict[_opr]
                        start = True; continue
                    else: # the CIGAR region is for continuing the loading
                        if total_length + _len <= length:
                            out_CIGAR += str(_len)+operation_dict[_opr]
                            total_length += _len
                        else:
                            out_CIGAR += str(length-total_length)+operation_dict[_opr]
                            break
            elif _opr in [1]:
                if moved_len + _len <= rand_seed and not start:
                    moved_len += _len# ; startPos += _len
                else:
                    if not start: # the site where start loading CIGAR information
                        if _len - (rand_seed - moved_len) < length:
                            out_CIGAR += str(_len - (rand_seed - moved_len)) + 'S'
                            total_length += _len - (rand_seed - moved_len)
                        else: 
                            out_CIGAR += str(length - (rand_seed - moved_len)) + 'S'
                        startPos += rand_seed - moved_len
                        start = True; continue
                    else: # the CIGAR region is for continuing the loading
                        if total_length + _len < length:
                            out_CIGAR += str(_len)+'I'
                            total_length += _len
                        else: 
                            out_CIGAR += str(length-total_length)+'S'
                            break
            elif _opr in [2,3]: # del, ref_skip
                if not start: startPos += _len
                else: out_CIGAR += str(_len)+operation_dict[_opr]
            elif _opr == 5: continue # hard_clip 
            else: return('UNKNOWN_CIGAR', startPos)
        return(out_CIGAR, startPos)

def seq_modification(seq, length=0, direction='left', rand_seed=0):

    if direction == 'left':
        return(seq[:length])
    elif direction == 'right':
        return(seq[::-1][:length][::-1])
    else: # direction == 'rand':
        return(seq[rand_seed:rand_seed+length])

def qual_modification(qual, length=0, direction='left', rand_seed=0):

    if direction == 'left':
        return(qual[:length])
    elif direction == 'right':
        return(qual[::-1][:length][::-1])
    else: # direction == 'random':
        return(qual[rand_seed:rand_seed+length])

if __name__ == '__main__':
    arg(sys.argv)
