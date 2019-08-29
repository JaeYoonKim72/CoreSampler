import os, sys, collections, time, copy
import itertools
import numpy as np

def ver_check():
    version_major = sys.version_info.major
    if version_major == 3: 
        itertools.imap = lambda *args, **kwargs: list(map(*args, **kwargs))
        print("Operating python version %s"%(version_major))
    else:
        itertools.imap = itertools.imap
        print("Operating python version %s"%(version_major))
        print("The genocore program only works with Python version 3.")
        sys.exit()

def now():
    now_time = time.strftime("%Y-%m-%d %p %H:%M:%S", time.localtime())
    return now_time

def dt(time_dt):
    if time_dt >= 60:
        time_dt = time_dt / 60
        time_int = "Min"
    else:
        time_dt = time_dt
        time_int = "Sec"
    run_time = str(round(time_dt,2)) 
    return run_time, time_int

def calc_maf(l):
    n0, n1, n2 = l.count('0'), l.count('1'), l.count('2')
    n = n0 + n1 + n2
    p = float(((2*n0)+n1)/(2*n))
    q = 1 - p 
    maf1 = min(p,q)
    mono = "TRUE" if maf1 <=0 else "FALSE"
    loh = "TRUE" if n1 <=0 else "FALSE"
    calc_maf_line = [n, n0, n1, n2, p, q, maf1, mono, loh]
    return calc_maf_line
    
def geno_count(geno_line):
    geno_set = np.array(list(set(geno_line)),dtype=float)
    geno_set_nan = len(geno_set[~np.isnan(geno_set)])
    return geno_set_nan

def overlap_score(snp_array):
    snp_array[np.isnan(snp_array)] = 99.0
    genotype, genotype_n = np.unique(snp_array, return_counts=True)
    genotype_dic = dict(zip(genotype, genotype_n))
    genotype_dic[99.0] = 0.0
    genotype_count = np.fromiter((genotype_dic[x] for x in snp_array), dtype=float)
    return genotype_count

def calc_dy(idxy, prenumxy, coverage1xy):
    if idxy == 0:
        dyy = coverage1xy
    elif idxy == 1 and prenumxy ==0:
        dyy = coverage1xy
    elif idxy == 1 and prenumxy !=0:
        dyy = np.float(coverage1xy) - np.float(coverage_table[1][2])
    elif idxy > 1 and prenumxy !=0:
        dyy = np.float(coverage1xy) - np.float(coverage_table[idxy][2])
    else:
        dyy = np.float(coverage1xy) - np.float(coverage_table[idxy-1][2])
    return dyy

def cover(idxx, prenumx, samplex, var_numx, ide_numx):
    if (idxx == 0) or (idxx == 1 and len(coverage_table) == 0):
        coverage_table.append(["Iteration","Sample_name", "Coverage", "Difference"])
        print(','.join(["Iteration","Sample_name", "Coverage", "Difference"]), file=Temp_file)
    np.seterr(divide='ignore', invalid='ignore')
    coverage1x = str(np.nanmean(ide_numx/var_numx*100))
    dyx = calc_dy(idxx, prenumx, coverage1x)
    coverage1_print = [str(idxx), samplex, coverage1x, str(dyx)]
    print(','.join(coverage1_print), file=Temp_file)
    coverage_table.append(coverage1_print)
    return coverage1_print

def CoreSampler(infile, preset, sample_n, maf, outname):
    print("Python Version Check...")
    ver_check()
    
    infile = infile
    pfile = preset
    nsamples = np.int(sample_n)
    output = outname
    maf = np.float(maf)

    ##variable setting
    var_num=[]
    ide_num=[]
    data_set = []
    coverage = []
    rnames = []
    cnames = []
    rm_maf = []
    nr = 0
    nc = 0
    result = []
    result_idx = []
    coreset = []
    prenum = 0
    global coverage_table
    coverage_table=[]
    Temp_file_name = output + '_TMP.csv'
    Cover_file_name = output + "_Coverage.csv"
    Coreset_file_name = output + "_CoreSample.csv"
    Maf_file_name = output + "_Removed_markers.csv"
    Core_sample_name = output + "_CoreSample.list.txt"
    Core_file = open(Core_sample_name, 'w')
    global Temp_file
    Temp_file = open(Temp_file_name, 'a')
    Cover_file = open(Cover_file_name, 'w')
    Coreset_file = open(Coreset_file_name, 'w')

    print("Data loading...", now())
    np_s_time = time.time()        
    for line2 in itertools.imap(lambda x:x.rstrip('\n').rstrip('\r').split(','), open(infile)):
        if line2[0] == "":
            cnames = line2[1:]
            continue
        cm = calc_maf(line2[1:])
        if maf != 0.0:
            if cm[6] <= maf:
                maf_line = [line2[0], str(cm[0]), str(cm[4]), str(cm[5])]
                rm_maf.append(maf_line)
                continue
        rnames.append(line2[0])
        preprocessed_line = [x if x in ["0", "1", "2"] else np.nan if x.isdigit() else np.nan for x in line2[1:]]
        var_num_t = np.array(list(set(preprocessed_line)), dtype=float)
        var_num_t = np.array(var_num_t)[~np.isnan(var_num_t)]
        var_num_t = var_num.append(len(var_num_t))
        data_set.append(preprocessed_line)
    else:
        nr = len(rnames)
        nc = len(cnames)
        var_num = np.array(var_num, dtype=float)
        data_set = np.array(data_set, dtype=float)
        counts = copy.deepcopy(data_set)
        counts_cnames = np.array(copy.deepcopy(cnames))
        np_e_time = time.time()
        np_dt = dt(np_e_time - np_s_time)
        if len(rm_maf) != 0:
            Maf_file = open(Maf_file_name, 'w')
            rm_maf = [["SNP name", "Total_Allele_counts", "Ref_Allele_freq", "Alt_Allele_frea"]] + rm_maf
            for marker in rm_maf:
                print(','.join(marker), file=Maf_file)
            else:
                Maf_file.close()
        print("Loading complete!")
        print("Loading time", np_dt[0], np_dt[1])
        print("Allocated memory size: " + str(round((sys.getsizeof(data_set)/1024),2)) + " Kbytes")
   
    if pfile != None:
        preset = [x for x in itertools.imap(lambda x:x.rstrip('\n').rstrip('\r'),open(pfile))] if pfile != "NN" else "NN"
        prenum = len(preset)
        select_counts_idx = [cnames.index(i) for i in preset]
        result.extend(preset)
        result_idx.extend(select_counts_idx)
        ide_num = np.apply_along_axis(geno_count, 1, data_set[:,result_idx])
        counts = np.delete(counts, select_counts_idx, 1)
        counts_cnames = np.delete(counts_cnames, select_counts_idx)
        coverage = cover(0, prenum, "preset", var_num, ide_num)
        print("Preset coverage is", str(round(np.float(coverage[2]),6)), "%")

    t_s_time = time.time()
    mpe=np.arange(1,nsamples-prenum+1).tolist()
    for idx in mpe:
        print("Selecting", idx, "-th sample", now())
        i_s_time = time.time()

        not_na_counts=np.apply_along_axis(lambda x:len(x)-np.sum(np.isnan(x)),0,counts)
        candidate_counts_idx = np.where(not_na_counts == np.max(not_na_counts))[0]
        candidate_counts_names = counts_cnames[candidate_counts_idx]
        step0=np.apply_along_axis(lambda x: overlap_score(x),1,counts)
        counts[counts == 99] = np.nan
        step01 = step0[:,candidate_counts_idx]
        overlap = np.apply_along_axis(np.mean, 0, step01)
    
        select = np.where(overlap == np.max(overlap))[0]
        select_name = candidate_counts_names[select]
        if len(select) != 1:
            minnum = np.apply_along_axis(lambda x: np.min(var_num[np.where(~np.isnan(x))]),0,counts[:,candidate_counts_idx[select]])
            minselidx = np.where(minnum == min(minnum))[0]
            select = select[minselidx]
            select = select[0]
            select_name = candidate_counts_names[select]
        select_counts_idx = candidate_counts_idx[select]
        if type(select_name) == np.str_:
            select_name = [select_name]
        select_dataset_idx = [cnames.index(x) for x in select_name]
        result.append(select_name[0])
        result_idx.append(select_dataset_idx[0])
        rm_idx = select_counts_idx
        ide_num = np.apply_along_axis(geno_count, 1, data_set[:,result_idx])
        coverage = cover(idx, prenum, select_name[0], var_num, ide_num)
        print("Selected sample is", select_name[0])   
        print("Coverage and Difference:", str(round(np.float(coverage[2]),6)), "%", str(round(np.float(coverage[3]),6)), "%")

        counts_cnames_total_index = np.arange(len(counts_cnames))
        counts_cnames_resi_del_rm = np.setdiff1d(counts_cnames_total_index, rm_idx)
        counts_mask_resi_del_rm = np.apply_along_axis(lambda x:np.isin(x[counts_cnames_resi_del_rm], x[rm_idx]),1,counts)
        counts = np.delete(counts, rm_idx, 1)
        counts[counts_mask_resi_del_rm] = np.nan
        counts_cnames = np.delete(counts_cnames, rm_idx)

        i_e_time = time.time()
        np_dt = dt(i_e_time - i_s_time)
        print("Running tme time", np_dt[0], np_dt[1])
        print('-'*50)
    else:
        t_e_time = time.time()
        np_dt = dt(t_e_time - t_s_time)
        print("Total running time", np_dt[0], np_dt[1])
        coreset = data_set[:,result_idx]
        coreset = np.array(np.array(coreset, dtype='i4'), dtype=str)
        coreset[coreset == '-2147483648'] = 'NA'
        coreset[coreset == '-9223372036854775808'] = 'NA'
        rnames = np.array(rnames).reshape((nr,1))
        coreset = np.hstack((rnames,coreset))
        print("Making Selected sample list file...")
        for sname in result:
            print(sname, file=Core_file)
        else:
            print('Complete')
            Core_file.close()
        result = np.array([""]+result)
        coreset = np.vstack((result,coreset))
        Temp_file.close()
        print("Making coverage tabel csv file...")
        for cover_line in coverage_table:
            print(','.join(cover_line), file=Cover_file)
        else:
            print("Complete")
            Cover_file.close()
        print("Making Selected samples csv file...")
        for core_line in coreset:
            print(','.join(core_line), file=Coreset_file)
        else:
            print("Complete")
            Coreset_file.close()  
        return "Done"
