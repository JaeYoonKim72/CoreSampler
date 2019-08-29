import sys, time, gzip
import numpy as np

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

def VCF_to_CSV(Infile, Outfile, phased, gziped):
    print("Converting VCF to CSV ...", now())
    s_time = time.time()    
    if phased == 'N':
        pn = '/'
    elif phased == 'Y':
        pn = '|'
    else:
        print("Invalid phased command")
        sys.exit()
    infile = Infile
#    outfile_name = infile.rstrip('.vcf') +'.csv'
    outfile_name = Outfile
    outfile = open(outfile_name, 'w')

    if gziped == "N":
        openfile = open(infile)
    elif gziped == "Y":
        openfile = gzip.open(infile, 'rb')
    else:
        print("Invalid gziped command")
        sys.exit()

    for line in map(lambda x:x, openfile):
        if gziped == "N":
            line = line.rstrip('\n').rstrip('\r').split('\t')
        else:
            line = line.decode().rstrip('\n').rstrip('\r').split('\t')
        if line[0].startswith('##'):
            continue
        if line[0].startswith('#'):
            data_name = line[9:]
            new_head = [''] + data_name
            print(','.join(new_head), file=outfile)
            continue
        CHR = line[0]
        POS = line[1]
        RS = line[2]
        REF = line[3]
        ALT = line[4]
        NEW_POS = CHR + '_' + POS +'_' + REF + ALT + '_' + RS

        new_data = [] 
        data = line[9:]

        for x in data:
            xx = x.split(':')[0].split(pn)
            if xx[0] == ".":
                new_x = int(-1)
            else:
                new_x = int(xx[0]) + int(xx[1])
            new_data.append(str(new_x))
        new_line = [NEW_POS] + new_data
        print(','.join(new_line), file=outfile)
    else:
        e_time = time.time()    
        np_dt = dt(e_time - s_time)
#        print("Coverting complete!")
        print("Running time", np_dt[0], np_dt[1])
        outfile.close()
        return "Done"
        
def VCF_select(VCF, Sample_list, outname, gziped):
    Sample_list_file = open(Sample_list)
    Sample_list = [x.rstrip('\n').rstrip('\r').split('\t')[0] for x in Sample_list_file]
    print("Selecting samples from the VCF file...", now()) 
    Total_idx = []
    outfile = open(outname, 'w')

    if gziped == "N":
        openfile = open(VCF)
    elif gziped == "Y":
        openfile = gzip.open(VCF, 'rb')
    else:
        print("Invalid gziped command")
        sys.exit()

    s_time = time.time()
    for line in map(lambda x:x, openfile):
        if gziped == "N":
            line = line.rstrip('\n').rstrip('\r').split('\t')
        else:
            line = line.decode().rstrip('\n').rstrip('\r').split('\t')
        if line[0].startswith('##'):
            print('\t'.join(line), file=outfile)
            continue
        if line[0].startswith('#'):
            Info_idx = np.arange(9)
            sample_idx = np.array([line.index(x) for x in Sample_list], dtype=int)
            Total_idx = np.hstack([Info_idx, sample_idx])
            line = np.array(line)
            Sel_line = line[Total_idx]
            print('\t'.join(Sel_line), file=outfile)
            continue
        line = np.array(line)
        Sel_line = line[Total_idx]
        print('\t'.join(Sel_line), file=outfile)
    else:
        outfile.close()
        e_time = time.time()
        np_dt = dt(e_time - s_time)
        print("Running time", np_dt[0], np_dt[1])
        return "Done"















