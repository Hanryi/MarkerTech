"""
Created on April 22, 2021

@Author
Han Rui (154702913@qq.com)

@Usage
> python3 Minima.py <i_csv> <o_xlsx.dir> <o_png.dir>
# <i_csv>: table containing frequency info (Num, Frequency, UMI)
#          1. Num: number of sample which has a record on the site
#          2. Frequency: appearance probability of the mut base provided by user
#          3. UMI: number of reads with unique UMI marker
# <o_xlsx.dir>: directory for saving genotype results
# <o_png.dir>: directory for saving KDE plots

@Function
Genotyping for each sample by minima dots in the KDE curve
"""

import os
import sys

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV
from sklearn.model_selection import LeaveOneOut

i_csv = sys.argv[1]
o_dir_xlsx = sys.argv[2]
o_dir_png = sys.argv[3]

# Read as a row vector
freqTab = pd.read_csv(i_csv, header=0, names=["Num", "Frequency", "UMI"])
num = ['Num'] + [i for i in freqTab.Num]
frequency = ['Frequency'] + [i for i in freqTab.Frequency]

freqSeq = freqTab.Frequency
freqArr = np.array(freqSeq)
step = freqArr.ptp()
# LOO of Cross-validation for best bandwidth of KDE
grid = GridSearchCV(
    estimator=KernelDensity(kernel='gaussian'),
    param_grid={'bandwidth': 10 ** np.linspace(-1.2, -0.3, 100)},
    cv=LeaveOneOut(),
)
grid.fit(freqArr.reshape(-1, 1))
bandwidth_best = grid.best_params_["bandwidth"]
bandwidth_best *= step

# Build suitable model by using best bandwidth
model = KernelDensity(bandwidth=bandwidth_best, kernel='gaussian')
model.fit(freqArr.reshape(-1, 1))
# X inputs and Y outputs
freq_range = np.linspace(freqSeq.min() - 1, freqSeq.max() + 1, 1000)
freq_log_prob = model.score_samples(freq_range.reshape(-1, 1))
freq_prob = np.exp(freq_log_prob)

# Find minima dots
xMini = []
yMini = []
for i in range(1, len(freq_range) - 1):
    if not 0 <= freq_range[i] <= 1:
        continue
    if freq_prob[i - 1] > freq_prob[i] and freq_prob[i + 1] > freq_prob[i]:
        xMini.append(freq_range[i])
        yMini.append(freq_prob[i])

# Plot
csvName = os.path.split(i_csv)[1].split('.')[0]
plt.title(csvName)
plt.xlim(0, 1)
sns.set()
sns.despine(top=True, right=False, offset=5)
# Rug
sns.distplot(
    freqArr, kde=False,
    rug=True,
    rug_kws={
        'color': '#2978b5', 'alpha': 0.5,
    },
    hist=False, norm_hist=True,
    hist_kws={
        'color': '#8BB0A6', 'alpha': 0.25,
    },
)
# KDE curve
plt.plot(freq_range, freq_prob, color='#2978b5', label='KDE')
plt.ylim(ymin=0)
plt.fill(freq_range, freq_prob, color='#fbe0c4', alpha=0.5)
# Minima dots
plt.scatter(x=xMini, y=yMini, s=16, c='#0061a8', label='minima')
plt.legend()
# Text
for i in range(len(xMini)):
    plt.text(
        xMini[i], yMini[i],
        str("%.3f" % (xMini[i])) + "\n",
        ha="center", va="bottom"
    )
# plt.show()

# Png
plt.savefig(
    os.path.join(
        o_dir_png, csvName + ".png"
    ),
    format='png',
    dpi=600,
)
plt.close()

# Genotype
genotype = ["Genotype"]
xZone = [0] + xMini + [1]
for i in freqSeq:
    n = 0
    for j in range(len(xZone) - 1):
        if xZone[j] <= i <= xZone[j + 1]:
            genotype.append(n)
        n += 1
genotypeMat = [num, frequency, genotype]
genotypeFrame = pd.DataFrame(genotypeMat).T
# Xlsx
genoTab = os.path.join(o_dir_xlsx, csvName + ".xlsx")
sheetName = 'Genotype'
frameWriter = pd.ExcelWriter(genoTab, engine='xlsxwriter')
genotypeFrame.to_excel(
    frameWriter, sheet_name=sheetName, header=False, index=False
)

""" Customize style.
"""
genotypeBook = frameWriter.book
genotypeSheet = frameWriter.sheets[sheetName]
headerStyle = genotypeBook.add_format({
    'valign': 'vcenter',
    'align': 'center',
    'color': '#2f5b66',
    # 'fg_color': '#F4B084',
    'font_name': 'Times New Roman',
    'bold':  True,
    'text_wrap': True,
})
defaultCellStyle = genotypeBook.add_format({
    'valign': 'vcenter',
    'align': 'center',
})
# Set row style include header row
genotypeSheet.set_row(0, 28, cell_format=headerStyle)
for row in range(1, genotypeFrame.shape[0]):
    genotypeSheet.set_row(row, 16)
# Set column style
genotypeSheet.set_column('A:C', cell_format=defaultCellStyle)
# Most chars in each column
charLen = [0 for _ in range(genotypeFrame.shape[1])]
for column in range(genotypeFrame.shape[1]):
    for sample in range(1, genotypeFrame.shape[0]):
        if len(str(genotypeFrame[column][sample])) + 2 > charLen[column]:
            charLen[column] = len(str(genotypeFrame[column][sample])) + 2
# Set column width (pixel)
genotypeSheet.set_column('A:A', max(6.64, charLen[0]))
genotypeSheet.set_column('B:B', max(10.64, charLen[1]))
genotypeSheet.set_column('C:C', max(9.64, charLen[2]))
# Save as file
frameWriter.save()
