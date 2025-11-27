#!/usr/bin/env python3

# ./filter_reads.py -f ____ -o _____ -m ______ -u _____ -s _____ -t _____ --ht _____ -q _____ --hq _____


# Load Packages 
import subprocess # allows python to use yq in the activated conda envs TNSeq_env (upstream of running this .py script)
import argparse 
import gzip 

# MANAGING INPUTS 
# defining input args 
def get_args():
    parser = argparse.ArgumentParser(description="QC script to remove low-quality reads from FASTQ file.")
    parser.add_argument("-f", required=True, help="Absolute path to UNZIPPED INPUT .fastq (pre-filter)") 
    parser.add_argument("-o", required=True, help="Absolute path to OUTPUT MATCHED READS .fastq (post-filter — GENOMIC DNA ONLY, TRANSPONS REMOVED — for Alignment)")
    parser.add_argument("-m", required=True, help="Absolute path to OUTPUT MATCHED READS .fastq (post-filter — TRANSPOSONS KEPT — for FASTQC REPORT)")
    parser.add_argument("-u", required=True, help="Absolute path to OUTPUT UNMATCHED READS .fastq (reads that were filtered out)")
    parser.add_argument("-s", required=True, help="Absolute path to OUTPUT Pipeline Summary .tsv (file will store calculated pipeline performance numbers to be used in the Final QC Report)")
    parser.add_argument("-t", required=True, help="Transposon Border Sequence (str)")
    parser.add_argument("-ht", required=True, help="Hamming Distance Threshold (# of mismatched bases permitted in TBS to be kept)")
    parser.add_argument("-q", required=True, help="Quality Score Threshold (Phred Score Threshold #)")
    parser.add_argument("-hq", required=True, help="Hamming Distance Threshold (# of low-quality bases permitted in sequence to be kept)")
    return parser.parse_args()

# setting input args--> variables 
args=get_args()
input_fastq_file=args.f
output_fastq_file=args.o
matched_fastq_file=args.m
unmatched_fastq_file=args.u
pipeline_summary_file=args.s
transpon_border_seq=args.t
ham_dist_threshold_tbs=int(args.ht)
qs_threshold=int(args.q)
ham_dist_threshold_qs=int(args.hq)

# DEFINING FUNCTIONS 
def convert_phred(letter: str) -> int:
    '''Converts a single character into a phred score'''
    return ord(letter) -33

def hamdist_qs(qs_line: str,ham_dist_threshold: int, qs_threshold: int):
    '''
    Assesses whether or not current qs_line meets qs-threshold depending on hamming_distance (# of pos.s permitted to fall below qs-threshold)
    Input: qs_line (str)[normally from R2 or R3 for demultiplex.], ham_dist threshold, qs_threshold 
    Output: a boolean (T/F) saying whether or not given qs_line passed the quality-score threshold 
    #IMPORTANT: REQUIRES bioinfo.convert_phred ( a fxn that converts a qs-line from ASCII-->actual phred scores)
    '''
    #Converting QS line (ASCII-->Phred (assumes +33 encoding))
    conv_qs=[]
    for letter in qs_line:
        phred_score=convert_phred(letter)
        conv_qs.append(phred_score)
    #Calc. Ham_Dist (a counter for # of bases w/ a qs that fell BELOW threshold)
        #Init. Ham_Dist 
    ham_dist=int(0)
        #cac. ham_dist based on qs threshold 
    for score in conv_qs:
        if score <=qs_threshold:
            ham_dist+=1
        #Filtering based on threshold
    if ham_dist >= ham_dist_threshold:
        return False
    elif ham_dist < ham_dist_threshold:
        return True

def hamdist_tbs(seq_line: str,ham_dist_threshold: int,tranposon_border_seq: str):
    '''
    Checks if transposon barcode is in the sequence, with leniency for mismatches
    Input: seq_line (str), ham_dist_threshold (max # of mismatched bases permitted to still be considered a "matched read"), tranposon_border_seq (LAST ~12 bases of transposon border seq.)
    Output: 
        - if tranposon seq. detected, will return the 0-based POS where the detected tranposon-border seq. ENDS in the seq. 
        - if tranposon seq. NOT detected, will return 'None' 
    '''
    # Looping over substrings the length of the tranposon_border_seq (sliding a k-mer window across the seq.) --> looking for a match 
        # Calc. k of kmer-window = length of TBS
    tbs_len=len(tranposon_border_seq)
        # Loop through each kmer-window in seq.
    for i in range(len(seq_line)-tbs_len+1): # +1 to make it 1-based
        tbs_window = seq_line[i:i+tbs_len]
        # Calc. Ham_Dist (# of mismatched bases) 
            # Init. Ham_Dist 
        ham_dist=int(0)
            # FOR every aligned pos. of a given barcode_window --> calc. ham_dist (# of mismatched bases)
                # zip(a,b) --> returns a tuple for EVERY aligned position b/w strings a & b
        for a, b in zip(tbs_window, tranposon_border_seq):
            #IF the tuple's bases mismatch --> Increment ham_dist +1 for given barcode_window 
            if a != b:
                ham_dist+=1
        # If hamming distance passes threshold (doesn't exceed permitted # of mismatches), return the position (0-based) where transposon barcode ends 
        if ham_dist <= ham_dist_threshold:
            return i + tbs_len #(0-based cutoff pos.)
    # If none of the kmer windows pass ham_dist_threshold, return None (False) 
    return None 

# 0. Set up 
    # a. Initialize Record-Variables 
record_holder = ["", "", "", ""] # empty LIST to temporarily hold processing records ( 4 lines @ a time) during each loop 
    # b. Init. Counter Variables 
counter = 0 # tracks which line we're on (within current record) — aka the index in record_holder list (which is 0-based)
unmatched_reads = 0 
matched_reads = 0 

# 1. Open input FASTQ File (reading), output filtered FASTQ File [transposon removed] (writing), output unmatched FASTQ File (writing), & output matched FASTQ File [transposons kept](writing)
# a. Open FASTQ file in data folder 
with gzip.open (input_fastq_file,"rt") as input_fastq, \
    open(output_fastq_file, "w") as output_fastq, \
    open(matched_fastq_file, "w") as matched_fastq, \
    open(unmatched_fastq_file, "w") as unmatched_fastq:
# 2. Build / Isolate record — readings lines until we've built a record 
    # a. Iterate through EACH LINE of FASTQ file (to avoid loading everything into mem.)
    for line in input_fastq:
        # remove leading & trailing charcters in string (line)
        line = line.strip() 
        # append line to record_holder[counter] — counter # specifies which element within the list to append to 
        record_holder[counter] = line 
        counter+=1 # increment record-line-counter +1 
        # ATP full record built 
        if counter ==4: 
            # set each item in list into variables that can be appended into a string later
            header = record_holder[0]
            seq_line = record_holder[1]
            plus_line = record_holder[2]
            qs_line = record_holder[3]
            # At this point we've built the entire record (as a list), and line variable is holding onto the \
            # header line of next record (will be added to the reset record at end of this if). 
#3. Determine if current record will be filtered out or kept + IF KEPT, REMOVE TBS / Isolate ONLY Genomic Seq. 
    # CALC. TBS MATCH HAMDIST. — Is this a valid transposon-insertion read? 
            # a. Determine if TBS detected in seq. + if so, store TBS end POS # 
            tbs_pos = hamdist_tbs(seq_line,ham_dist_threshold_tbs,transpon_border_seq)
    # CALC. QS HAMDIST. + DETERMINE IF RECORD NEEDS TO BE FILTERED OUT 
            
            # FILTERED RECORDS --> UNMATCHED FASTQ 
            # a. Filter out records that do not meet either of the hamdist. filters 
            if tbs_pos == None or hamdist_qs(qs_line,ham_dist_threshold_qs,qs_threshold)==False:
                unmatched_fastq.write("\n".join(record_holder)+"\n") # filtered out record --> list is written to unmatched FASTQ file as a joined cont. string, ending w/ & delimited by "\n"s
                unmatched_reads+=1 # increment unmatched reads tracker
            
            # KEPT RECORDS --> MATCHED FASTQ 
            # For records that passed filters.... 
            else: 
                # b. Write original record w/ transposons KEPT to output FASTQ --> FASTQC 
                matched_fastq.write("\n".join(record_holder)+"\n")

                # c. Isolate only genomic DNA in seq. line  
                genomic_seq = seq_line[tbs_pos:] #slice string from tbs_pos to end 
                # d. Isolate only genomic DNA's quality scores in qs line 
                trimmed_qs_line = qs_line[tbs_pos:] #slice line from tbs_pos to end 
                # e. Replace seq + QS lines in record_holder list w/ trimmed versions 
                record_holder[1] = genomic_seq
                record_holder[3] = trimmed_qs_line
                # f. Write updated record w/ tranposons REMOVED (genomic DNA only) to output FASTQ --> Alignment 
                output_fastq.write("\n".join(record_holder)+"\n")
                matched_reads+=1 # increment matched reads tracker
                
            # Reset Record & Start the New Record (add header-line) + Reset counter to 1 
            record_holder = ["", "", "", ""]
            counter=0 

# Closing Files 
input_fastq.close()
output_fastq.close()
unmatched_fastq.close()
matched_fastq.close()

# CALC. & REPORT % MATCHED READS --> FINAL REPORT STAT 
    # Calc. Values 
percent_matched_reads= (matched_reads/(matched_reads+unmatched_reads))*100 
percent_reads_removed= (unmatched_reads/(matched_reads+unmatched_reads))*100 
    # Write into File 
with open (pipeline_summary_file,"w") as psf:
    psf.write("metric\tvalue\n")
    psf.write(f"percent_matched_reads\t{percent_matched_reads}\n")
    psf.write(f"percent_reads_removed\t{percent_reads_removed}\n")
psf.close()



        


