
import os, io, sys, pymzml, argparse, time, math

from pyteomics import fasta, parser, mass

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from scipy.stats import zscore, ttest_ind, f_oneway

pd.options.mode.chained_assignment = None
IAA_MOD_MASS = mass.calculate_mass(formula='CH2CONH')

argparser = argparse.ArgumentParser(description = 'Find signals in MS1 data matching tryptic peptides')

# file IO
argparser.add_argument('fasta',
                    type = str,
                    help = 'FASTA sequence file for target protein')
argparser.add_argument('mzml',
                    type = str,
                    help = 'Directory containing one or more mzML files to be analysed')

argparser.add_argument('--outDirectory',
                    type = str,
                    help = 'Directory to which output file should be written. Defaults to mzml directory if not specified')
argparser.add_argument('--outPrefix',
                    type = str,
                    help = 'File name prefix given to generated outputs. Defaults to a timestamp if not specified')
argparser.add_argument('--groups',
                    nargs = '+',
                    help = 'Treatment group labels used for determination of statistical significance. These labels are compared to the input mzML file\
                       names to assign individual files to treatment groups. If not given, all files are treated independently and statistical testing not applied. \
                       Group labels must not contin spaces.')

# MS matching
argparser.add_argument('--ppm',
                    type = float,
                    default = 5,
                    help = 'MS1 ppm tolerance for expected monoisotopic masses. Window is +/- ppm')
argparser.add_argument('--minFrac',
                    type = float,
                    default = 0.3,
                    help = 'Minimum fraction (0-1) of MS1 spectra in which a given target ion must be present for further consideration in the analysis')

# target generation properties`
argparser.add_argument('--missed',
                    type = int,
                    default = 1,
                    help = 'Number of missed enzymatic cleavages to allow for in silico protein digestion')
argparser.add_argument('--minLength',
                    type = int,
                    default = 5,
                    help = 'Minimum peptide length for consideration in analysis')
argparser.add_argument('--maxCharge',
                    type = int,
                    default = 4,
                    help = 'Maximum peptide charge to be considered in analysis')
argparser.add_argument('--minCharge',
                    type = int,
                    default = 1,
                    help = 'Minimum peptide charge to be considered in analysis')

# clustering
argparser.add_argument('--numClusters',
                    type = int,
                    default = 5,
                    help = 'Number of clusters to produce in HCA analysis')
argparser.add_argument('--noZtransform',
                   action = 'store_true',
                   help = 'Do not perform z-transform on peptide intensity values')

# data
argparser.add_argument('--normaliseIntensities',
                   action = 'store_true',
                   help = 'Normalise target intensity values to spectral maxima')
argparser.add_argument('--noLog',
                   action = 'store_true',
                   help = 'Do not log-transform target intensity values.')
argparser.add_argument('--pvalThreshold',
                    type = float,
                    default = 1,
                    help = 'Negative Log 10 p-value threshold required of peptides for inclusion')

# reporting
argparser.add_argument('--reportAll',
                    action = 'store_true',
                    help = 'If given, all charge states for all peptides will be listed separately in the outputs. \
                       Otherwise, only the most abundant charge state of each peptide in each group will be reported',
                    )

options = argparser.parse_args()

class Target(object):

    def __init__(self, peptide, charge, mods_mass):

        self.charge = charge
        self.peptide = peptide

        # calculate neutral mass of peptide
        self.target = mass.calculate_mass(
                sequence = self.peptide, ion_type = 'M', charge = charge) + float(mods_mass) / float(charge)

        # calculate upper and lower m/z limits
        self.targetLL = self.target - self.target / 1000000 * options.ppm
        self.targetHL = self.target + self.target / 1000000 * options.ppm

        self.targetIntensityDIct = {}
        self.targetScanCounter = {}
        return

    def setIntensity(self, file, value):
        try:
            self.targetIntensityDIct[file].append( value )
        except KeyError:
            self.targetIntensityDIct[file] = [ value ]
        return

    def incrementTargetScanCounter(self, file):
        try:
            self.targetScanCounter[file] += 1
        except:
            self.targetScanCounter[file] = 1
        return

    def getAverageIntensity(self, file):
        values = [float(v) for v in self.targetIntensityDIct[file]]
        try:
            return float(sum(values)) / float(len(values))
        except ZeroDivisionError:
            return 0

class Peptide(object):

    def __init__(self, peptide):
        self.peptide = peptide
        self.peptideLength = len(peptide)

        self.peptideMzDict = {}

        # add IAA Cys residues?
        self.neutral_unmodified = mass.calculate_mass(
                sequence = self.peptide, ion_type = 'M', charge = 0)

        self.neutral_modification_masses = self.peptide.lower().count('c') * IAA_MOD_MASS
        self.total_neutral_mass = self.neutral_unmodified + self.neutral_modification_masses

        self.targetList = []

        for z in range(options.minCharge, options.maxCharge + 1):
            self.targetList.append(
                Target(
                    self.peptide, z, self.neutral_modification_masses
                )
            )

        return

    def findChargeStatesToReport(self, groupMap, groups):

        valuesToReturn = []

        quantSum = 0
        for group in groups:
            files = groupMap[group]

            # need to get the intensities for each target for each file in this group
            groupDataList = []
            for t in self.targetList:
                groupData = {}

                # save target for later
                groupData['target'] = t
                groupData['charge'] = t.charge

                # average intensity excludes 0 values for now - data is fairly bouncy - do we want this?
                groupData['abundances'] = [t.getAverageIntensity(f) for f in files]
                groupData['numNonZero'] = len([_ for _ in groupData['abundances'] if _ != 0])

                # not all targets have non-zero intensity
                try:
                    groupData['avgAbundance'] = float(sum(groupData['abundances'])) #/ float(len(groupData['abundances']))
                except ZeroDivisionError:
                    groupData['avgAbundance'] = 0

                groupDataList.append(groupData)

            # sort targets based on average abundance - highest first
            groupDataList = sorted(groupDataList, key = lambda i: i['avgAbundance'], reverse = True)

            # first value is max avg intensity
            mostIntenseTarget = groupDataList[0]

            quantVals = []
            for f in files:
                # quite a bit of bouncy-ness to the data
                # --- require that any peak be present in > some fraction of all scans
                pointValues = [pv for pv in mostIntenseTarget['target'].targetIntensityDIct[f] if pv != 0]
                if len(pointValues) < len(mostIntenseTarget['target'].targetIntensityDIct[f]) * options.minFrac:
                    quantVals.append(0)
                else:
                    # not all file have same # of scans - average
                    quantVals.append(mostIntenseTarget['target'].getAverageIntensity(f))

            if self.peptide == 'ADLAKYICDNQDTISSK':
                print(group)
                print(len(pointValues), options.minFrac)
                print(groupDataList)

            quantSum += sum(quantVals)
            # create list of values to print
            valuesToReturn.extend(quantVals)
            valuesToReturn.append(mostIntenseTarget['target'].charge)
            valuesToReturn.append(mostIntenseTarget['target'].target)
            valuesToReturn.append(mostIntenseTarget['target'].targetLL)
            valuesToReturn.append(mostIntenseTarget['target'].targetHL)
        return valuesToReturn, quantSum

def digetsProteinFromFASTA():
    sequenceIter = fasta.read(source = options.fasta)
    uniquePeptides = set()
    for s in sequenceIter:
        newPeptides = parser.cleave(s.sequence, 'trypsin', missed_cleavages = options.missed, min_length = options.minLength)
        uniquePeptides.update(newPeptides)

    uniquePeptides = list(uniquePeptides)
    return [Peptide(x) for x in uniquePeptides]

def findPeptideIntensities():

    # digest sequence and return list of Peptide objects for each peptide
    peptides = digetsProteinFromFASTA()

    files = os.listdir(options.mzml)
    dataFiles = sorted([x for x in files if '.mzml' in x.lower()])

    for dataFile in dataFiles:
        print('\tProcessing %s' %dataFile)

        mzml_file = os.path.join(options.mzml, dataFile)
        run = pymzml.run.Reader(mzml_file)

        for n, spec in enumerate(run):
            mzs = spec.mz
            ints = spec.i

            scanMaxIntensity = ints.max()

            for p in peptides:
                for t in p.targetList:
                    mask = np.where (
                        (mzs > t.targetLL) & (mzs < t.targetHL))

                    if len(ints[mask]) == 0:
                        t.setIntensity(dataFile, 0)
                    else:
                        # max values align with fresstyle ion counts better
                        # makes result easier to check
                        # is there any downside here?
                        dataPoint = ints[mask].max()

                        if options.normaliseIntensities:
                            dataPoint = ints[mask].max() / scanMaxIntensity * 100
                        if not options.noLog:
                            if dataPoint >= 1:
                                dataPoint = math.log10(dataPoint)
                            else:
                                dataPoint = 0

                        t.setIntensity(dataFile, dataPoint)

                    t.incrementTargetScanCounter(dataFile)

    of1 = io.StringIO()

    if options.reportAll:
        toPrint = ['peptide', 'neutral_sequence_mass', 'neutral_modifications_mass', 'neutral_peptide_mass', 'length', 'charge','mz', 'mz_lower', 'mz_upper']
        toPrint += dataFiles
        of1.write('%s\n' %('\t'.join([str(x) for x in toPrint])))
        for p in peptides:
            for t in p.targetList:
                toPrint = [
                    p.peptide, p.neutral_unmodified, p.neutral_modification_masses, p.total_neutral_mass,
                    p.peptideLength, t.charge, t.target, t.targetLL, t.targetHL
                ]
                quantVals = []
                for dataFile in dataFiles:
                    # quite a bit of bouncy-ness to the data
                    # --- require that any peak be present in > some fraction of all scans

                    pointValues = [pv for pv in t.targetIntensityDIct[dataFile] if pv != 0]
                    if len(pointValues) < len(t.targetIntensityDIct[dataFile]) * options.minFrac:
                        quantVals.append(0)
                    else:
                        # not all file have same # of scans - average
                        quantVals.append(t.getAverageIntensity(dataFile))

                # don't write any peptides/charge states that are not observed
                # in at least one sample
                if sum(quantVals) == 0: continue

                toPrint += quantVals
                of1.write('%s\n' %('\t'.join([str(x) for x in toPrint])))
    else:
        # headers will look different when charge data are summarised to peptides
        # -- need to report charge used for each group
        groupMap = makeGroupMap(dataFiles)
        groups = list(groupMap.keys())

        toPrint = ['peptide', 'neutral_sequence_mass', 'neutral_modifications_mass', 'neutral_peptide_mass', 'length']

        for group in groups:
            toPrint.extend(groupMap[group])
            toPrint.extend( [
                '%s charge' %group, '%s mz' %group, '%s mz_lower' %group, '%s mz_upper' %group
            ])

        of1.write('%s\n' %('\t'.join([str(x) for x in toPrint])))

        for p in peptides:
            toPrint = [ p.peptide, p.neutral_unmodified, p.neutral_modification_masses, p.total_neutral_mass, p.peptideLength ]
            valuesToWrite, quantSum  = p.findChargeStatesToReport(groupMap, groups)

            if quantSum == 0: continue

            toPrint.extend(valuesToWrite)
            of1.write('%s\n' %('\t'.join([str(x) for x in toPrint])))

    # return pointer to beginning of file
    of1.seek(0)

    # create dataframe from reulsts
    df = pd.read_csv(of1, delimiter = '\t')

    if options.reportAll:
        # keys denote a given peptide in a given charge state
        # these need to be unique to display properly in following plots
        df['key'] = df['peptide'] + '_+' + df['charge'].astype(str)
        df = df.set_index('key')
    else:
        # not present if charge states are summarised
        # default back to peptide sequence
        df['key'] = df['peptide']
        df = df.set_index('key')

    print('\tDone...')
    return df, dataFiles, peptides

def makeGroupMap(dataFiles):
    groupMap = {}

    if options.groups:
        for g in options.groups:
            for dataFile in dataFiles:
                if g in dataFile:
                    try:
                        groupMap[g].append(dataFile)
                    except:
                        groupMap[g] = [dataFile]
    else:
        # no groupings specified
        # treat each sample as its own group
        for dataFile in dataFiles:
            groupMap[dataFile] = [dataFile]

    return groupMap

def doClustering(df, dataFiles, resultsDir):

    # create subset from quantification columns
    quantDF, cbLabel, quantCols = getQuantDataFrame(df)

    # for Ali - 2 March 21
    quantDF = quantDF.fillna(0)

    # assign samples to groups if specified
    if options.groups:

        groupMap = makeGroupMap(dataFiles)

        print('\n\tSample group assignment')
        print('\t-----------------------')
        for k,v in groupMap.items():
            for sample in v:
                print('\tGroup: %s, Sample: %s' %(k, sample))

        significantRows = []
        for index, row in quantDF.iterrows():
            groupedRowData = [row[entries].to_list() for entries in groupMap.values()]
            f, pval = f_oneway(*groupedRowData)
            nlogp = math.log10(pval) * -1

            if nlogp > options.pvalThreshold:
                significantRows.append(index)

        print('\n\t%s targets before statistical filtering' %len(quantDF))
        quantDF = quantDF[quantDF.index.isin(significantRows)]
        df = df[df.index.isin(significantRows)]
        print('\t%s targets after statistical filtering' %len(quantDF))
    else:
        print('\n\t%s targets identified' %len(quantDF))

    cluster = AgglomerativeClustering(n_clusters = options.numClusters)
    cluster.fit_predict(quantDF)

    groupLabels = list(set(cluster.labels_))

    # create color palette for the user-specified number of clusters
    # These are the colors that appear alongside the dendrograms
    lut = dict(zip(groupLabels, sns.color_palette("hls", options.numClusters)))

    # map these colors to the relevant site entries
    df['cluster'] = cluster.labels_
    clusters = df["cluster"]

    row_colors = clusters.map(lut)

    plot = sns.clustermap(
        quantDF, cmap = 'coolwarm', cbar_kws={'label': cbLabel}, row_colors = row_colors)

    # write figure to file
    outFigure = os.path.join(resultsDir, options.outPrefix + 'clustermap')

    savefig(plot, plot, outFigure)

    # write data frame to excel file
    outExcelFile = os.path.join(resultsDir, options.outPrefix + 'datatable.xlsx')

    df = df.sort_values('cluster')
    writer = pd.ExcelWriter(outExcelFile, engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Sheet1')

    workbook  = writer.book
    worksheet = writer.sheets['Sheet1']

    # create hex strings for colours to use in shading output table cells
    colorDict = {}
    for k,v in lut.items():
        values = [int(float(_*255)) for _ in list(v)]
        colorDict[k] = '#%02x%02x%02x' % tuple(values)

    formatBook = {}
    for g in groupLabels:
        formatBook[g] = workbook.add_format({'bg_color': colorDict[g]})

    # need to get index of cluster column to find which excel column should be shaded
    clusterColIndex = list(df).index('cluster')

    start_row = 1
    start_col = clusterColIndex + 1 # excel is 1-indexed
    end_row = len(df)
    end_col = start_col

    for key, formati in formatBook.items():
        worksheet.conditional_format(start_row, start_col, end_row, end_col,
            {'type':     'cell',
            'criteria': '=',
            'value':    key,
            'format':   formati})

    writer.save()
    print('\tDone...')
    return df, lut

def doAnalysis(clusterMatrix, lut, resultsDir):

    quantDF, cbLabel, quantCols = getQuantDataFrame(clusterMatrix)

    targetColumns = ['cluster', 'length', 'charge','mz', 'neutral_peptide_mass']

    subClusterMatrix = clusterMatrix[targetColumns]
    densityPlot = sns.pairplot( data=subClusterMatrix, hue = 'cluster', palette = lut)
    savefig(densityPlot, densityPlot, os.path.join(resultsDir, options.outPrefix + "pairs.png"))

    sns.set_style("darkgrid")
    fig, ax = plt.subplots()
    catPlot = sns.countplot(x="charge",  data=clusterMatrix, hue = 'cluster', palette = lut, ax = ax)
    ax.legend(loc='upper right')
    savefig(catPlot, fig, os.path.join(resultsDir, options.outPrefix + "cluster_charge_histogram.png"))
    print('\tDone...')
    return

def getQuantDataFrame(df):
    quantCols = [x for x in list(df) if '.mzml' in x.lower()]

    # create subset from quantification columns
    quantDF = df[quantCols]
    print(quantDF)

    if not options.noZtransform:
        quantDF_T = quantDF.T
        quantDF_T_Z = quantDF_T.apply(zscore)
        quantDF_Z = quantDF_T_Z.T
        quantDF = quantDF_Z
        cbLabel = 'Z-score'
    else:
        cbLabel = 'Raw Intensity'

    return quantDF, cbLabel, quantCols

def savefig(plot, fig, path, resolution = 400):

    try:
        plot.savefig(path, dpi = resolution)
    except:
        try:
            plot = plot.get_figure()
            plot.savefig(path, dpi = resolution)
        except:
            plot.figure.savefig(path, dpi = resolution)

    try:
        plt.close(fig)
    except:
        try:
            plt.close(plot)
        except:
            # special method for seaborn clustermap
            # the mpl ax object is hiding in an ax_heatmap attribute
            # https://stackoverflow.com/questions/32868423/plot-on-top-of-seaborn-clustermap
            try:
                plt.close(plot.ax_heatmap.get_figure())
            except:
                # hope for the best
                return

def appendToDict(d, key, value):
    try:
        d[key].append(value)
    except:
        d[key] = [value]
    return

def doPCAs(peptides, dataFiles, resultsDir):

    pcaData = {}

    for dataFile in dataFiles:

        for i in range(len(peptides[0].targetList[0].targetIntensityDIct[dataFile])):
            appendToDict(pcaData, 'dataFile', dataFile)

        for p in peptides:
            for t in p.targetList:
                for v in t.targetIntensityDIct[dataFile]:
                    key = '%s_%s'%(p.peptide, t.charge)
                    appendToDict(pcaData, key, v)

    pcaDataFrame = pd.DataFrame(pcaData, columns = list(pcaData.keys()))

    Y = pcaDataFrame.pop('dataFile')
    X = pcaDataFrame

    # returns first N_compoents principal components
    pca = PCA(n_components=2)

    # trailing .transform() needed to return actual values
    # otherwise returns object
    # gives array of num features (samples) rows by num components cols
    X_r = pca.fit(X).transform(X)

    # to shade points in plot by saple group
    # need to subset plots
    # get truth mask for which rows in X-r array correspond to each group
    fig, ax = plt.subplots()

    groups = set(list(Y.to_list()))
    for group in groups:
        mask = Y.isin([group])
        pcaPlot = ax.scatter(X_r[mask,0], X_r[mask,1], label = group)

    ax.legend()
    ax.set_xlabel('PC1 (%.1f%%)' %(pca.explained_variance_ratio_[0]*100))
    ax.set_ylabel('PC2 (%.1f%%)' %(pca.explained_variance_ratio_[1]*100))

    plt.savefig(os.path.join(resultsDir, options.outPrefix + 'pca.png'))
    print('\tDone...')
    return

def main():

    # if not output path specified - write to same directory as mzML files
    if not options.outDirectory:
        options.outDirectory = options.mzml

    if options.outPrefix:
        options.outPrefix += '_'
    else:
        options.outPrefix = ''

    # create results directory
    resultsDir = os.path.join(options.outDirectory, 'resutls' )
    try:
        os.mkdir( resultsDir )
    except FileExistsError:
        pass

    print('Processing sample data')
    print('======================')
    # read peptides and return dataframe containing peptides and intensity values
    peptideIntensityMatrix, dataFiles, peptides = findPeptideIntensities()

    print('\n\nRunning PCA')
    print('===========')
    # pca plots
    doPCAs(peptides, dataFiles, resultsDir)

    print('\n\nClustering Feature Data')
    print('=======================')
    # do clustering
    clusterMatrix, lut = doClustering(peptideIntensityMatrix, dataFiles, resultsDir)

    if options.reportAll:
        print('\n\nRunning Summary Analysis')
        print('========================')
        #  analyse
        doAnalysis(clusterMatrix, lut, resultsDir)

    print('\n\nAll Done!')

if __name__ == '__main__':
    main()
