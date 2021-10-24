import pandas as pd
import os
import numpy as np
from sklearn.ensemble import RandomForestClassifier

def train_classifier(filename):
    training_data = read_from_h5(filename, '/training_data')

    labels = training_data["label"]

    del training_data["label"]

    clf = RandomForestClassifier(n_estimators=500)

    model = clf.fit(training_data, labels)

    features = training_data.columns.values
    model.feature_names_ = features

    return model

def extract_insert_metrics(file):
    ''' Extract median and mean insert size '''

    # picardtools insertmetrics completes with code 0 and doesn't generate metrics file
    # if inputs don't have sufficient read count
    if not os.path.isfile(file):
        return 0, 0, 0

    mfile = open(file)
    
    targetlines = []
    
    line = mfile.readline()
    
    while line != '':
        if line.startswith('## METRICS CLASS'):
            targetlines.append(mfile.readline().strip().split('\t'))
            targetlines.append(mfile.readline().strip().split('\t'))
            break
        line = mfile.readline()

    mfile.close()

    header, data = targetlines

    header = [v.lower() for v in header]
    header = {v:i for i,v in enumerate(header)}

    median_ins_size = data[header['median_insert_size']]
    mean_ins_size = data[header['mean_insert_size']]
    std_dev_ins_size = data[header['standard_deviation']]

    return median_ins_size, mean_ins_size, std_dev_ins_size

def extract_flagstat_metrics(file):
    """
    extract from flagstat
    """

    df = pd.read_table(file,
                        sep=r'\s\+\s0\s',
                        header=None,
                        names=['value', 'type'],
                        engine='python')

    tot_reads = df[df['type']=='in total (QC-passed reads + QC-failed reads)']['value']
    tot_mpd_reads = df[(df['type'].str.contains('mapped') == True ) & ( df['type'].str.contains('mate mapped') == False)]
    tot_dup_reads = df[df['type']=='duplicates']['value']
    tot_prop_paired = df[df['type'].str.contains('properly paired') ]

    assert len(tot_reads) == 1
    assert len(tot_mpd_reads) == 1
    assert len(tot_dup_reads) == 1
    assert len(tot_prop_paired) == 1

    tot_reads = tot_reads.iloc[0]
    tot_mpd_reads = tot_mpd_reads['value'].iloc[0]
    tot_dup_reads = tot_dup_reads.iloc[0]
    tot_prop_paired = tot_prop_paired['value'].iloc[0]


    return tot_reads, tot_mpd_reads, tot_dup_reads, tot_prop_paired

def format_data(metrics,
              colnames):

    data = []
    for colname in colnames:
        coldata = metrics[colname]

        if colname == 'scaled_halfiness':
            # haploid poison adds inf, replace with big number since 0 is considered good
            # and we want to score to decrease
            coldata = coldata.replace(np.inf, 1e10)
        data.append(coldata)

    data = pd.concat(data, axis=1)

    data = data.replace(-np.inf, np.nan)
    data = data.fillna(0)

    return data

def read_from_h5(filename, tablename):
    with pd.HDFStore(filename) as h5store:
        data = h5store[tablename]
    return data

metrics = pd.read_csv(snakemake.input.metrics)
flgstat = extract_flagstat_metrics(snakemake.input.flgstat)
insrt = extract_insert_metrics(snakemake.input.insrt)

metrics["total_reads"] = flgstat[0]
metrics["total_mapped_reads"] = flgstat[1]
metrics["total_duplicate_reads"] = flgstat[2]
metrics["total_properly_paired"] = flgstat[3]
metrics["percent_duplicate_reads"] = metrics["total_duplicate_reads"] / metrics["total_mapped_reads"]
metrics["median_insert_size"] = insrt[0]
metrics["mean_insert_size"] = insrt[1]
metrics["standard_deviation_insert_size"] = insrt[2]

print(snakemake.config)

#add quality score
model = train_classifier(snakemake.config['classifier_training_data'])
data = format_data(metrics, model.feature_names_)
predictions = model.predict_proba(data)
metrics["quality"] = predictions[0][1]

metrics.to_csv(snakemake.output.metrics)