
import os, io, sys, pymzml, argparse, time

from pyteomics import fasta, parser, mass

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from sklearn.cluster import AgglomerativeClustering
from scipy.stats import zscore

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
        if not value:
            try:
                assert file in self.targetIntensityDIct.keys()
                return
            except:
                self.targetIntensityDIct[file] = []
                return
        try:
            self.targetIntensityDIct[file].append( value )
        except KeyError:
            self.targetIntensityDIct[file] = [ value ]
        return

    def updateTargetScanCounter(self, file):
        try:
            self.targetScanCounter[file] += 1
        except KeyError:
            self.targetScanCounter[file] = 1
        return

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
        print('Processing %s' %dataFile)

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
                        t.setIntensity(dataFile, None)
                    else:
                        # max values align with fresstyle ion counts better
                        # makes result easier to check
                        # is there any downside here?

                        if options.normaliseIntensities:
                            t.setIntensity(dataFile, ints[mask].max() / scanMaxIntensity * 100)
                        else:
                            t.setIntensity(dataFile, ints[mask].max())

                    t.updateTargetScanCounter(dataFile)

    of1 = io.StringIO()

    toPrint = ['peptide', 'neutral_sequence_mass', 'neutral_modifications_mass', 'neutral_peptide_mass', 'length', 'charge','mz', 'mz_lower', 'mz_upper']
    toPrint += dataFiles
    of1.write('%s\n' %('\t'.join([str(x) for x in toPrint])))

    for p in peptides:
        for t in p.targetList:
            toPrint = [p.peptide, p.neutral_unmodified, p.neutral_modification_masses, p.total_neutral_mass, p.peptideLength, t.charge, t.target, t.targetLL, t.targetHL]
            quantVals = []
            for dataFile in dataFiles:
                # quite a bit of bouncy-ness to the data
                # --- require that any peak be present in > some fraction of all scans
                if len(t.targetIntensityDIct[dataFile]) < t.targetScanCounter[dataFile] * options.minFrac:
                    quantVals.append(0)
                else:
                    # not all file have same # of scans - average
                    quantVals.append(int(sum(t.targetIntensityDIct[dataFile])/len(t.targetIntensityDIct[dataFile])))

            # don't write any peptides/charge states that are not observed
            # in at least one sample
            if sum(quantVals) < 1: continue

            toPrint += quantVals
            of1.write('%s\n' %('\t'.join([str(x) for x in toPrint])))

    # return pointer to beginning of file
    of1.seek(0)

    # create dataframe from reulsts
    df = pd.read_csv(of1, delimiter = '\t')

    # keys denote a given peptide in a given charge state
    # these need to be unique to display properly in following plots
    df['key'] = df['peptide'] + '_+' + df['charge'].astype(str)
    df = df.set_index('key')

    return df


def doClustering(df):

    # create subset from quantification columns
    quantDF, cbLabel, quantCols = getQuantDataFrame(df)

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
    outFigure = os.path.join(
        options.outDirectory, options.outPrefix + '_clustermap')

    savefig(plot, plot, outFigure)

    # write data frame to excel file
    outExcelFile = os.path.join(options.outDirectory, options.outPrefix + '_datatable.xlsx')

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
    return df, lut

def doAnalysis(clusterMatrix,lut):

    quantDF, cbLabel, quantCols = getQuantDataFrame(clusterMatrix)

    # density Plots
    densityDirectory = os.path.join(options.outDirectory, 'density' )
    try:
        os.mkdir( os.path.join(options.outDirectory, 'density' ))
    except FileExistsError:
        pass

    targetColumns = ['cluster', 'length', 'charge','mz', 'neutral_peptide_mass']

    subClusterMatrix = clusterMatrix[targetColumns]
    densityPlot = sns.pairplot( data=subClusterMatrix, hue = 'cluster', palette = lut)
    savefig(densityPlot, densityPlot, os.path.join(densityDirectory, "pairs.png"))

    sns.set_style("darkgrid")
    fig, ax = plt.subplots()
    catPlot = sns.countplot(x="charge",  data=clusterMatrix, hue = 'cluster', palette = lut, ax = ax)
    ax.legend(loc='upper right')
    savefig(catPlot, fig, os.path.join(densityDirectory, "cluster_charge_histogram.png"))

def getQuantDataFrame(df):
    quantCols = [x for x in list(df) if '.mzml' in x.lower()]

    # create subset from quantification columns
    quantDF = df[quantCols]

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

def main():

    # if not output path specified - write to same directory as mzML files
    if not options.outDirectory:
        options.outDirectory = options.mzml

    # assign datestamp to output files if no prefix given
    if not options.outPrefix:
        options.outPrefix = str(int(time.time()))

    # read peptides and return dataframe containing peptides and intensity values
    peptideIntensityMatrix = findPeptideIntensities()

    # do clustering
    clusterMatrix, lut = doClustering(peptideIntensityMatrix)

    #  analyse
    doAnalysis(clusterMatrix, lut)

if __name__ == '__main__':
    main()
