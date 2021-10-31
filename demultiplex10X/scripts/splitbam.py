import os
import csv
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)

#adapted from scope https://github.com/rujinwang/SCOPE/blob/master/inst/docs/split_script.py

data = snakemake.params.barcodes
logging.info(data)
tcount=0
bcount=0
with open(snakemake.input.bam) as file:
    for line in file:
      tcount+=1
      line=line.strip()
      tags = line.split()[11:]
      tags_sort = sorted(tags)
      for tag in tags_sort:
          tag_split=tag.split(':')
          if 'CB' in tag_split:
            barcode=tag_split[-1]
            if barcode in data:
              bcount+=1
              logging.info('Barcode present in bam file ({})'.format(barcode))
              logging.info(str(tcount)+" "+str(bcount))
              directory = "results/align/"
              outfile = directory+"/"+barcode+".sam"
              try:
                f=open(outfile,"a+")
                print >> f,line
                f.close()
              except:
                os.mkdir(directory)
                f=open(outfile,"a+")
                print >> f,line
                f.close()
              break

pct=float(bcount)/float(tcount)
print(pct)