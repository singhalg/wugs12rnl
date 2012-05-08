#-------------------------------------------------------------------------------
# Name:        module1
# Purpose:
#
# Author:      gsinghal
#
# Created:     07/05/2012
# Copyright:   (c) gsinghal 2012
# Licence:     <your licence>
#-------------------------------------------------------------------------------


import sys, re, os
import psutil as PS
from subprocess import Popen, PIPE, STDOUT
import pickle
import time



def fastqDownload(sampleName):
    fn = 'gsinghal'+sampleName+'.log'
    log = open(fn, 'w')

    file_loc_NCBI = "ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/"
    file_loc_EBI = "ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/"

    fh = open('YRI_lowCoverage.csv', 'rU')

    file_PID_dict = {}

    gunzipProcs = []
    switch = True
    seq_data = fh.readlines()
    fh.close()
    psList = {}

    for line in seq_data:
        flds = line.split(',')
        if (flds[9].strip() == sampleName) and (flds[20].strip()== '0')  and (flds[12].strip()=='ILLUMINA'):

            if switch:

                file_loc = file_loc_NCBI + flds[0].strip()
            else:
                file_loc = file_loc_EBI + flds[0].strip()

            if switch:
                switch=False
            else:
                switch=True

            cmd = " wget --retry-connrefused -b "+ file_loc
            print cmd
            log.write(cmd)
            log.write('\n')
            job = Popen(cmd,  shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)

            output = job.stdout.read()

            s = output.find('pid') + 4
            e = output.find('.')
            PID = output[s:e]

            fq_file_name = ret_FQfile_name(file_loc)
            file_PID_dict[PID] = fq_file_name


            proc = PS.Process(int(PID))

            psList[PID] = proc

            time.sleep(10)
    #print file_PID_dict
    processes = file_PID_dict.keys()
    for each in processes:
        p = psList[each]
        p.wait()
        fileName = file_PID_dict[each]
        cmd = 'gunzip  ' + fileName
        job = Popen(cmd,  shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
        print cmd
        log.write(cmd)
        log.write('\n')
        gunzipProcs.append(job)
        time.sleep(10)

    for aproc in gunzipProcs:
        aproc.wait()

    log.write('File download complete')
    log.write('gunzip extraction complete')
    os.system(' bash mail_finish_download ')
    log.close()

def bashScript(nodes, ppn, wt, working_dir, reference, samplename):

    fn = 'gsinghal'+sampleName+'.log'
    log = open(fn, 'a')


    jobname = 'gsinghal_' + samplename + '_'+reference
    bsName1 = jobname+ '_paired'
    bsName2 = jobname+ '_se'

    fhout1 = open(bsName1, 'w')

    options = "#!/bin/bash" + '\n'+ "#PBS -l nodes="+nodes+":ppn="+ppn+",walltime="+wt+ ',mem=10gb \n' + "#PBS -N "+jobname + '\n' + "#PBS -d " + working_dir + '\n' + '#PBS -m abe '+'\n' +'#PBS -q dque_smp ' +'\n'

    fhout1.write(options)
    fhout2.write(options)
    single_end = ''

    fh = open('YRI_lowCoverage.csv', 'rU')
    seq_data = fh.readlines()
    fh.close()

    unpaired = ''
    mate1 = ''
    mate2 = ''


    for line in seq_data:
        flds = line.split(',')
        if (flds[9].strip() == samplename) and (flds[20].strip()=='0') and (flds[12].strip()=='ILLUMINA'):
            file_loc = flds[0].strip()
            fq_file_name = ret_FQfile_name2(file_loc)
            if flds[18].strip() == 'PAIRED':
                print fq_file_name
                if re.search('_1', fq_file_name):
                    mate1 += fq_file_name + ','
                elif re.search('_2', fq_file_name):
                    mate2 += fq_file_name + ','
                else :
                    unpaired += fq_file_name + ','
            elif  flds[18].strip() == 'SINGLE':
                single_end += fq_file_name + ','


    proc = str(int(nodes)*int(ppn))
# writing bowtie command for paired mates
    outFile1 = samplename+'_'+reference+'_paired.sam'
    stdoutFile = samplename+'_'+reference+'_paired_bowtie2.txt  '

    command = 'bowtie2 -p ' + proc +'  --sensitive -x  ' + reference+ '  -1 ' + mate1[:-1] + ' -2 ' + mate2[:-1] +  ' -S '+outFile1 +' &>  '+ stdoutFile

    fhout1.write(command)
    fhout1.write('\n')
    fhout1.write('bash mail_bowtie_paired \n')
    log.write('Writing bowtie commands to the bash script')
    log.write(command)
    log.write('\n')

# writing bowtie command for single ended reads
    outFile2 = samplename+'_'+reference+'_se.sam'
    stdoutFile = samplename+'_'+reference+'_se_bowtie2.txt  '

    command = 'bowtie2 -p ' + proc +'  --sensitive -x  ' + reference+ '  -U ' + single_end[:-1] +  ' -S '+outFile2 +' &>  '+ stdoutFile

    fhout2.write(command)
    fhout2.write('\n')

    log.write(command)
    log.write('\n')

    fhout2.write('bash mail_bowtie_se \n')



# samtools view -bS  NA18487_YRI.sam > NA18487_YRI.bam \n '
    command = 'samtools view -bS  ' + outFile1 + ' > ' + outFile1[:-4] + '.bam'
    fhout1.write(command)
    fhout1.write('\n')
    fhout1.write('bash mail_samtools_paired \n')
    log.write('Writing samtools commands to the bash script')
    log.write(command)
    log.write('\n')


    command = 'samtools view -bS  ' + outFile2 + ' > ' + outFile2[:-4] + '.bam'
    fhout2.write(command)
    fhout2.write('\n')
    fhout2.write('bash mail_samtools_se \n')
    log.write(command)
    log.write('\n')
    fhout1.write('sleep 10')
    fhout2.write('sleep 10')
    log.close()
    fhout1.close()
    fhout2.close()
    return bsName1, bsName2


def runBS(bsName):
    cmd = "qsub " + bsName

    os.system(cmd)


def ret_FQfile_name(file_loc):
    start = file_loc.find('read')
    if start < 0:
        print 'COULD NOT DOWNLOAD ', file_loc
    else:
        s = start+ 5

    return file_loc[s:].strip()

def control(sampleName):

    fastqDownload(sampleName)
    #bashScript(nodes, ppn, wt, working_dir, reference, samplename)
    bsName1, bsName2 = bashScript('1', '48', '12:00:00' , '/scratch/gsinghal/YRI_LC/' , 'YRIref_index', sampleName)
##    runBS(bsName1)
##    runBS(bsName2)
    bsName1, bsName2 = bashScript('1', '48', '12:00:00' , '/scratch/gsinghal/YRI_LC/' , 'hg19', sampleName)
##    runBS(bsName1)
##    runBS(bsName2)

def main():
    control('NA18498')


if __name__ == '__main__':
    main()
