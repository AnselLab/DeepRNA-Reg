#### Author: Harshaan Sekhon (sskhon2014@berkeley.edu)


import argparse
import sys
import numpy as np
from datascience import *
from scipy.integrate import simps
from scipy import stats
import pandas as pd
from scipy.signal import savgol_filter
from scipy.stats import percentileofscore

import time
import multiprocess
from multiprocess import Pool
import tqdm
import subprocess
from io import BytesIO
    
class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def DEEP_CLIP(loci):

    from keras.models import load_model
    
    pbar = tqdm.tqdm(total=len(loci))
    
    

    def get_depth_data(track_files,track_names,chrom,start,stop,strand,track_type):
        def view_region(track_file,strand,region):
            return subprocess.Popen(("samtools", "view",
                                     strand_to_flag[use_strand],
                                     "-b",track_file,
                                     region),stderr=subprocess.PIPE,stdout=subprocess.PIPE)
        mydepths = pd.DataFrame([0]*(stop-start+1),index=range(start,stop+1),columns=["depth"])
        depth_list = pd.DataFrame(0,index=range(start,stop),columns=track_names)
        strandinvert = {"+":"-","-":"+"}
        strand_to_flag = {"+":"-F 0x10",
            "-":"-f 0x10"}
        for n,track_file in enumerate(track_files):
            use_strand=strand
            region = chrom + ":" + str(start) + "-" + str(stop)
            if track_type[n] == "as":
                use_strand = strandinvert[strand]
            # Get sequences from a given region (in binary bam format still)
            ps =view_region(track_file,strand_to_flag[use_strand],region)
            sout,err = ps.communicate() # get stdout, stderr
            ## CHECK TO MAKE SURE THE REFERENCE GENOME CHROMOSOME IS FINE.
            if len(err)>0: # is there anytihn in stder?
                if b"specifies an unknown reference name" in err:
                    # SWITCH REFERENCE
                    temp_chrom = chrom.replace("chr","")
                    region = temp_chrom + ":" + str(start) + "-" + str(stop)
                    ps =view_region(track_file,strand_to_flag[use_strand],region)
                    sout,err = ps.communicate()
            if len(err)>0:
                raise NameError("Unknown samtools error. Ran: samtools view %s -b %s %s | samtools depth - " % (strand_to_flag[use_strand],track_file,region))
            # Run samtools depth on the sequences retrieved
            ps2 = subprocess.Popen(("samtools", "depth","-"),stdin=subprocess.PIPE,stdout=subprocess.PIPE)
            output,err = ps2.communicate(input=sout)
            sample_depths = pd.read_table(BytesIO(output),names=["chrom","depth"],index_col=1)
            if len(sample_depths.index)>0:
                mydepths.depth = sample_depths.depth
                depth_list[track_names[n]] = sample_depths.depth
                depth_list = depth_list.fillna(value=0)
        return depth_list
    def ranges(nums):
        nums = sorted(set(nums))
        gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s+1 < e]
        edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
        summ = list(zip(edges, edges))
        start = [summ[x][0] for x in np.arange(len(summ))]
        stop = [summ[x][1] for x in np.arange(len(summ))]

        return Table().with_columns("start", start, "stop", stop, "range", [i-j for i,j in zip(stop,start)])

    def range_generator(summary_tbl):

        cond_1 = ranges(summary_tbl.where("Greater Binding in Condition 1", are.equal_to(True)).column(0))
        cond_2 = ranges(summary_tbl.where("Greater Binding in Condition 2", are.equal_to(True)).column(0))
        return cond_1.where("range", are.above(3)), cond_2.where("range", are.above(3)), summary_tbl

    ### We are throwing away any regions called that have a lenth lower than 4 nucleotides

    
    def MA_normalization(read_data):
        M = np.log(read_data.column(1) + .999999) - np.log(read_data.column(2) + .999999)
        A = (np.log(read_data.column(1) + .999999) + np.log(read_data.column(2) + .999999))
        slope, intercept, r_value, p_value, std_err =  stats.linregress(A,M)
        adjusted_M = M - ((slope*A)+intercept)

        return adjusted_M


    model = load_model(model_name)


    def get_pred(arr):
        a = ((arr- np.mean(arr))/np.std(arr))
        b = ((arr- np.mean(arr))/np.std(arr))
        seq = [[x,y] for x,y in zip(a,b)]
        seq = np.array(seq)
        X, y = seq[:, 0], seq[:, 1]
        X = X.reshape((len(X), 1, 1))
        return X, y


    def auc_generator(condition, summ_tbl):
        depths = []
        startZ = condition.column("start")
        stopZ = condition.column("stop")
        for i in np.arange(condition.num_rows):
            start = startZ[i]
            stop = stopZ[i]
            depth = simps(summ_tbl.where(0, are.between(start, stop+1)).column("Adjusted_Fitted_M"))
            if depth <0:
                depth = -1*depth # We multiply by -1 because we don't care wether integrated area is positive or negative
                                # We only care for the intensity (magnitude) of the integral
            depths.append(depth)
        return np.asarray(depths)
    
    cond1_starts, cond1_stops, cond1_strands, cond1_chroms, cond1_geneIDs, cond1_depths = [], [], [], [], [], []
    cond2_starts, cond2_stops, cond2_strands, cond2_chroms, cond2_geneIDs, gene_runs, cond2_depths = [], [], [], [], [], [], []


    for locus in loci:
        pbar.update(1)
        chrom = locus[0]
        start = int(locus[1])
        stop = int(locus[2])
        geneid = locus[3]
        strand = locus[5]
        region = chrom + ":" + str(start) + "-" + str(stop)
        
        depths = get_depth_data(make_array(f1, f2),make_array(f1, f2),chrom,start,stop,strand,["s", "s"])
        f1_depths = depths[f1].tolist()
        f2_depths = depths[f2].tolist()

        reading_data = Table().with_columns("Unnamed: 0",depths.index.tolist(), f1,f1_depths, f2, f2_depths)

        
        if sum(reading_data.column(1)) > 0 and sum(reading_data.column(2)) > 0 and reading_data.num_rows >21:
            difference = reading_data.column(2) - reading_data.column(1)
            if sum(difference  == difference[0]) != len(difference): #There cannot be a constant difference in read depth
                if sum(np.isnan(MA_normalization(reading_data))) != len(reading_data.column(1)):
                    predictor = savgol_filter(MA_normalization(reading_data), 21,3)
                    reading_data = reading_data.with_column("m_val", predictor )#.select(0,3)
                    vals = []
                    X,y = get_pred(predictor)
                    yhat = model.predict([X], verbose=0)

                    if sum(np.isnan(yhat))[0] != len(yhat):
                        # If we get an array of all NANs, ignore it.
                        yhat = savgol_filter([yhat[z][0] for z in np.arange(len(yhat))], 21,1)
                        for b in np.arange(len(yhat)):
                            val = yhat[b]
                            min_pos = np.argmin([abs(1.-val), abs(val), abs(-1-val)])
                            if min_pos ==0:
                                val = 1
                            if min_pos ==1:
                                val = 0
                            if min_pos ==2:
                                val = -1
                            vals.append(val)
                        reading_data = reading_data.with_columns("vals", vals)
                        ##PLOT DATA TEST
                        ##reading_data.plot(0,width=15, height=7)
                        reading_data = reading_data.with_column("Greater Binding in Condition 1",reading_data.column("vals")==1 )
                        reading_data = reading_data.with_column("Greater Binding in Condition 2",reading_data.column("vals")==-1 )
                        tbl = reading_data.relabeled("Unnamed: 0", "pos").relabeled("m_val", "Adjusted_Fitted_M").relabeled("vals", "lstm predictions")
                        cond1, cond2, summ_tbl = range_generator(tbl)
                        cond1 = cond1.with_column("strand", [strand for i in cond1.column(0)]).with_column("chrom", [chrom for i in cond1.column(0)]).with_column("geneid", [geneid for z in cond1.column(0)])
                        cond2 = cond2.with_column("strand", [strand for i in cond2.column(0)]).with_column("chrom", [chrom for i in cond2.column(0)]).with_column("geneid", [geneid for z in cond2.column(0)])

                        [cond1_geneIDs.append(i) for i in cond1.column("geneid")]
                        [cond1_chroms.append(i) for i in cond1.column("chrom")]
                        [cond1_starts.append(i) for i in [int(y) for y in cond1.column("start")]]
                        [cond1_stops.append(i) for i in [int(y) for y in cond1.column("stop")]]
                        [cond1_strands.append(i) for i in cond1.column("strand")]
                        [cond1_depths.append(i) for i in auc_generator(cond1, summ_tbl)]

                        [cond2_geneIDs.append(i) for i in cond2.column("geneid")]
                        [cond2_chroms.append(i) for i in cond2.column("chrom")]
                        [cond2_starts.append(i ) for i in [int(y) for y in cond2.column("start")]]
                        [cond2_stops.append(i) for i in [int(y) for y in cond2.column("stop")]]
                        [cond2_strands.append(i) for i in cond2.column("strand")]
                        [cond2_depths.append(i) for i in auc_generator(cond2, summ_tbl)]


    condition_1 = Table().with_columns("geneid", cond1_geneIDs, "chrom",cond1_chroms, "start",[int(y) for y in cond1_starts], "stop", [int(y) for y in cond1_stops],"strands",cond1_strands, "AUC Differential Binding",cond1_depths )
    condition_2 = Table().with_columns("geneid", cond2_geneIDs, "chrom",cond2_chroms, "start",[int(y) for y in cond2_starts], "stop", [int(y) for y in cond2_stops],"strands",cond2_strands, "AUC Differential Binding", cond2_depths)
    
    return [condition_1.to_df(),condition_2.to_df()]

def q_C_list(cond_2,DBE_col):
    arr = cond_2.column(DBE_col)

    # pre-sort array
    arr_sorted =  sorted(arr)

    # calculate percentiles using scipy func percentileofscore on each array element
    s = pd.Series(arr)
    percentiles = s.apply(lambda x: percentileofscore(arr_sorted, x))
    def quality_score(pct):
        if pct <= 25:
            return 1
        if pct <= 50:
            return 2
        if pct <=70:
            return 3
        if pct <= 80:
            return 4
        if pct <= 85:
            return 5
        if pct <= 90:
            return 6
        if pct <= 95:
            return 7
        if pct <= 97.5:
            return 8
        if pct <= 99:
            return 9
        if pct <= 100:
            return 10
    cond_2 = cond_2.with_columns("DBE Percentile", percentiles, "DBE Score", [quality_score(p) for p in percentiles])
    return cond_2

def concat_files(TABLEZ):
    starts, stops, strands, chroms, geneIDs, depths = [], [], [], [], [], []
    for tbl in TABLEZ:
        [starts.append(i) for i in tbl.column("start")]
        [stops.append(i) for i in tbl.column("stop")]
        [strands.append(i) for i in tbl.column("strands")]
        [chroms.append(i) for i in tbl.column("chrom")]
        [geneIDs.append(i) for i in tbl.column("geneid")]
        [depths.append(i) for i in tbl.column("AUC Differential Binding")]
    
        
    tbl = Table().with_columns("chrom", chroms,"start", starts, "stop", stops,"geneid",geneIDs,"Differential Binding Enrichment(DBE)",depths, "strand", strands )
    tbl = q_C_list(tbl, 'Differential Binding Enrichment(DBE)')
    return tbl
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description='Welcome to DeepRNAreg, a powerful tool for differential CLIP-Seq analysis.')
    parser.add_argument("--bam1", help="This is the absolute path to the first BAM file.")
    parser.add_argument("--bam2", help="This is the absolute path to the second BAM file.")
    parser.add_argument("--bed", help="This is the absolute path to a BED format file containing the genomic loci upon which DeepRNAreg should be executed.")


    args = parser.parse_args()
    bed_filename = args.bed
    f1 = args.bam1
    f2 = args.bam2
    if None in [bed_filename, f1, f2]:
        #sys.exit(f"{bcolors.WARNING}  PLace Warning Here{bcolors.UNDERLINE}")
        sys.exit(f"\033[31;1;4mError: One of the required files is not supplied.\033[0m")

        #sys.exit(''+ '\033[93m')
        
    model_name = 'DeepRNAreg.h5'

    logo = '\n\nWelcome to...\n\n\n##########      ##########    ##########  ##########              ##########  ####        ####   ############\n############    ##########    ##########  ####   ####             ##########  ####        ####   ####     ####\n##############  ####          ####        ####    ####            ####        ####        ####   ####      ####\n####      ####  ####          ####        ####    ####            ####        ####        ####   ####      ####\n####      ####  ##########    ##########  ####   ####             ####        ####        ####   ####     ####\n####      ####  ##########    ##########  ##########              ####        ####        ####   ############\n####      ####  ####          ####        ####                    ####        ####        ####   ####\n##############  ####          ####        ####                    ####        ####        ####   ####\n############    ##########    ##########  ####                    ##########  ##########  ####   ####\n##########      ##########    ##########  ####                    ##########  ##########  ####   ####\n\n\n'
        
    #print(f"{bcolors.HEADER} Welcome...{bcolors.ENDC}")
    print(logo)
    print(f"{bcolors.HEADER}BED-File: " + bed_filename)
    print('First BAM File: '+f1)
    print('Second BAM File: '+ f2, end='\n\n\n')
    
    #
    print(f"{bcolors.ENDC}Running Analysis...")
    
    
    loci = np.loadtxt(bed_filename, usecols=range(6),dtype=str)
    

    print((f"{bcolors.OKBLUE}Number of CPU Cores Detected: #num{bcolors.ENDC}").replace('#num', str(multiprocess.cpu_count())))
    start_time = time.time()
    with Pool(processes=multiprocess.cpu_count(),maxtasksperchild=20) as pool:
        results = pool.map(DEEP_CLIP, np.array_split(loci, multiprocess.cpu_count()) )
    stop_time = time.time()


    pool.close()
    pool.join()

    print("--- Runtime: %s seconds ---" % (stop_time - start_time))
    results = [[Table.from_df(j) for j in i] for i in results]
    cond1 =concat_files([i[0] for i in results])
    cond2 =concat_files([i[1] for i in results])
    cond1.to_csv('enriched_cond1.csv')
    cond2.to_csv('enriched_cond2.csv')
    print(f"{bcolors.HEADER}{bcolors.UNDERLINE}\n\n\n\n\n\nDeepRNAreg ran successfully! Output saved to enriched_cond1.csv and enriched_cond2.csv.\n\nGoodbye! \n\n{bcolors.ENDC}")
