"""
Created on March 13, 2022

@Author
Han Rui (154702913@qq.com)

@Usage
> python3 MethHist.py <meth.csv> <o_FencePlot.d>
# <meth.csv>: a sorted table (by chr and pos), its header containing
#             (Chr, Pos, Context, C, C+T, Active)
#             1. Chr: chromosome name
#             2. Pos: methylation position
#             3. Context: methylation type
#             4. C: unique cytosine
#             5. C+T: sum of unique thymine and cytosine
#             6. Active: active fraction
# <o_FencePlot.d>: a directory using to save fence plot figure
#                  (Fence Plot: named by @Author)

@Function
Convert a large amount of position and methylation active information in the
input table into easy-to-understand picture.

@Note
Fence Plot: A single fence plot showing as one subplot in the output figure.
            Its name comes from its horizontal discrete axes, because the space
            between the axes is similar to the space between the planks.
"""

import os
import sys
import matplotlib.pyplot as plt

import seaborn as sns


class MethData:
    """ Structured methylation data.

    This class mainly build structured data by the method named context_len().
    Moreover, when the object is instantiated, all the name and limit interval
    of chromosomes are added to class attribute.
    """
    lim_pref = "limit_"
    chromosome = []

    def __init__(self, filepath: str):
        self.path = filepath
        self.name = os.path.basename(filepath).split('.')[0]

        # Limits of upper and lower position in chromosome
        region_tab = open(self.path)
        next(region_tab)
        for line in region_tab:
            line_vec = line.split(',')
            if not hasattr(self, self.lim_pref + line_vec[0]):
                self.chromosome.append(line_vec[0])
                setattr(
                    self, self.lim_pref + line_vec[0],
                    [int(line_vec[1])] * 2
                )
            lim = getattr(self, self.lim_pref + line_vec[0])
            if int(line_vec[1]) > lim[1]:
                setattr(
                    self, self.lim_pref + line_vec[0],
                    [lim[0], int(line_vec[1])]
                )

    def context_ten(self):
        """ Convert the input file to a data tensor.

        :return: tensor
            format: {<context_1>: <chr_1>[(<pos>, <active>), ()], [],}
            e.g.  : {
                        'CG': [
                                [(20, 0.8), ...],
                                ...
                              ],
                        ...
                    }
        """
        region_tab = open(self.path)
        next(region_tab)
        meth_tensor = {}
        for line in region_tab:
            line_vec = line.split(',')
            if line_vec[2] not in meth_tensor.keys():
                meth_tensor[line_vec[2]] = [[] for _ in self.chromosome]
            meth_tensor[line_vec[2]][
                self.chromosome.index(line_vec[0])
            ].append(
                (int(line_vec[1]), eval(line_vec[5]))
            )
        return meth_tensor


class FencePlot(MethData):
    """ Fence plot.

    Preprocess tensor data and map limitation position to true interval when
    this class is instantiated. And the method named active_hist() can plot a
    single fence plot with a definite axes in the figure.
    """
    inr_pref = "interval_"
    dpi = 600
    single_void = 25
    sub_axis_level = -0.025
    pillar_width = 2
    # Font size
    endpoint_size = 10
    chr_name_size = 15
    act_name_size = 15
    # Font family
    base_font = 'FangSong'
    demo_font = 'Times New Roman'

    def __init__(self, filepath: str):
        # Allow Chinese and negative sign
        plt.rcParams['font.sans-serif'] = [self.base_font]
        plt.rcParams['axes.unicode_minus'] = False

        def step(start, length):
            return [start, start + length]

        super(FencePlot, self).__init__(filepath)
        # Adapt limitation to interval position
        length_seq = map(
            lambda p: p[-1] - p[0],
            (getattr(self, self.lim_pref + ch) for ch in self.chromosome)
        )
        for item in range(len(self.chromosome)):
            # True coordinate interval
            if item == 0:
                setattr(
                    self, self.inr_pref + self.chromosome[0],
                    step(self.single_void, next(length_seq))
                )
                continue
            setattr(
                self, self.inr_pref + self.chromosome[item],
                step(
                    getattr(
                        self, self.inr_pref + self.chromosome[item - 1]
                    )[1] + 2 * self.single_void, next(length_seq)
                )
            )
        self.x_max_value = getattr(
            self, self.inr_pref + self.chromosome[-1]
        )[1] + self.single_void
        self.meth_data = self.context_ten()

    def fence_color(self):
        """ Generate Palette.

        :return: (list of RGB tuple)
        """
        color_pal = sns.color_palette(
            palette="Set1", n_colors=len(self.chromosome), desat=1
        )
        return color_pal

    def active_hist(self, axes, sub_fence):
        """ Plot a fence subplot.

        :param axes: (Object)
        :param sub_fence: (list of tuple)

        :return: None
        """
        # Size of axis
        plt.axis(
            [0, self.x_max_value, -0.1, 1]
        )

        # Remove tick of x-axis and set step of y-tick
        axes.get_xaxis().set_visible(False)
        axes.yaxis.set_major_locator(
            plt.MultipleLocator(0.25)
        )
        # Remove frame and plot a line to replace y-axis
        for key, spine in axes.spines.items():
            spine.set_color(None)
        plt.plot([0, 0], [0, 1], linewidth=1, color='black')
        # Plot discrete chromosome axis
        for ch in self.chromosome:
            plt.plot(
                getattr(self, self.inr_pref + ch), [self.sub_axis_level] * 2,
                color=self.fence_color()[self.chromosome.index(ch)],
                linestyle="--", linewidth=0.5, alpha=0.75,
            )
            plt.scatter(
                getattr(self, self.inr_pref + ch), [self.sub_axis_level] * 2,
                color=self.fence_color()[self.chromosome.index(ch)], s=2
            )
            plt.text(
                getattr(self, self.inr_pref + ch)[0], -0.05,
                getattr(self, self.inr_pref + ch)[0], size=self.endpoint_size,
                horizontalalignment='center', verticalalignment='top',
                color=self.fence_color()[self.chromosome.index(ch)],
                bbox=dict(facecolor="white", alpha=0)
            )
            plt.text(
                getattr(self, self.inr_pref + ch)[1], -0.05,
                getattr(self, self.inr_pref + ch)[1], size=self.endpoint_size,
                horizontalalignment='center', verticalalignment='top',
                color=self.fence_color()[self.chromosome.index(ch)],
                bbox=dict(facecolor="white", alpha=0)
            )

        # Fence plot
        for ch in self.chromosome:
            for p in sub_fence[self.chromosome.index(ch)]:
                plt.plot(
                    [
                        getattr(self, self.inr_pref + ch)[0] + (
                                p[0] - getattr(self, self.lim_pref + ch)[0]
                        )
                    ] * 2, [0, p[1]],
                    linewidth=self.pillar_width, alpha=0.75,
                    color=self.fence_color()[self.chromosome.index(ch)]
                )
            # Chromosome name annotation
            plt.text(
                sum(getattr(self, self.inr_pref + ch)) / 2, -0.1, ch,
                horizontalalignment='center', verticalalignment='top',
                bbox=dict(facecolor="white", alpha=0), size=self.chr_name_size
            )
        return None


if __name__ == "__main__":
    fenceObj = FencePlot(sys.argv[1])
    context = list(fenceObj.meth_data.keys())

    # Figure size and subplot spacing
    fig = plt.figure(
        figsize=(fenceObj.x_max_value / 1000 * 16, 4 * len(context))
    )
    fig.subplots_adjust(hspace=0.25, wspace=0)
    # Plot fence plots
    for n in range(1, len(context) + 1):
        ax = fig.add_subplot(len(context), 1, n)
        ax.set_ylabel(
            context[n - 1] + ' Active',
            fontdict=dict(
                fontsize=fenceObj.act_name_size,
                family=fenceObj.demo_font,
                weight='bold'
            )
        )
        fenceObj.active_hist(ax, fenceObj.meth_data[context[n - 1]])

    plt.savefig(
        os.path.join(sys.argv[2], fenceObj.name + ".png"),
        dpi=fenceObj.dpi
    )
