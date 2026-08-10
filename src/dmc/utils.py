import sys
from time import strftime
import numpy as np
from scipy import stats
import logging
import matplotlib.pyplot as plt
from matplotlib import rc
from scipy import optimize

# use latex formatting for plots
rc('text', usetex=True)


def plot_known_predicted_ages(known_ages,
                              predicted_ages,
                              outfile,
                              title="EpigeneticPacemaker",
                              xlab='Chronological Age',
                              ylab='EPM Age',
                              font_size=16):
    # define optimization function
    def func(x, a, b, c):
        return a * np.asarray(x)**0.5 + c
    # fit trend line
    popt, pcov = \
        optimize.curve_fit(func, [1 + x for x in known_ages], predicted_ages)
    # get r squared
    rsquared = r2(predicted_ages, func([1 + x for x in known_ages], *popt))
    # format plot label
    plot_label = \
        f'$f(x)={popt[0]:.2f}x^{{1/2}} {popt[2]:.2f}, R^{{2}}={rsquared:.2f}$'
    # initialize plt plot
    fig, ax = plt.subplots(figsize=(12, 12))
    # plot trend line
    ax.plot(sorted(known_ages),
            func(sorted([1 + x for x in known_ages]),
            *popt), 'r--', label=plot_label)
    # scatter plot
    ax.scatter(known_ages, predicted_ages, marker='o', alpha=0.8, color='k')
    ax.set_title(title, fontsize=font_size + 2)
    ax.set_xlabel(xlab, fontsize=font_size)
    ax.set_ylabel(ylab, fontsize=font_size)
    ax.tick_params(axis='both', which='major', labelsize=font_size)
    ax.legend(fontsize=font_size)
    plt.savefig(outfile)



def r2(x, y):
    # return r squared
    return stats.pearsonr(x, y)[0]**2


def pearson_correlation(meth_matrix, phenotype):
    """
    calculate pearson correlation coefficient between rows of input matrix
    and phenotype
    """
    # calculate mean for each row and phenotype mean
    matrix_means = np.mean(meth_matrix, axis=1)
    phenotype_mean = np.mean(phenotype)
    transformed_matrix = meth_matrix - matrix_means.reshape([-1, 1])
    transformed_phenotype = phenotype - phenotype_mean
    covariance = np.sum(transformed_matrix * transformed_phenotype, axis=1)
    variance_meth = np.sqrt(np.sum(transformed_matrix ** 2, axis=1))
    variance_phenotype = np.sqrt(np.sum(transformed_phenotype ** 2))
    return covariance / (variance_meth * variance_phenotype)


def printlog(mesg):
    '''print progress message'''
    mesg = "@ " + strftime("%Y-%m-%d %H:%M:%S") + ": " + mesg
    print(mesg, file=sys.stderr)


def traof(x, adult_age=20):
    x = (x+1)/(1 + adult_age)
    if x <= 1:
        y = np.log(x)
    else:
        y = x - 1
    return y


def anti_traof(x, adult_age=20):
    if x < 0:
        y = (1 + adult_age)*np.exp(x) - 1
    else:
        y = (1 + adult_age)*x + adult_age
    return y


def plot_coef(infile, outfile, rfile):
    """
    Visualize CpGs' coefficients in the model. Missing CpGs are highlighted.

    Parameters
    ----------
    infile : str
        the name of input file (*.coef.tsv)
    outfile : str
        The name of output figure file. The suffix of this file must be \
        '.pdf' or '.png'.
    rfile : str
        The name of R script file. This file generates the above figure file.

    Returns
    -------
    None.

    """
    ROUT = open(rfile, 'a')
    print("\n\n###### Rscript to generate coefficients of clock CpGs ######\n\n", file=ROUT)
    print("d = read.table(file='%s', sep='\\t', header=T)" % infile, file=ROUT)
    print("d_rank = d[order(d$Coef),]", file=ROUT)
    print("rownames(d_rank) = 1:length(d_rank$Coef)", file=ROUT)
    print("index = as.numeric(rownames(d_rank))", file=ROUT)
    print("", file=ROUT)

    print("y_lab = 'Coefficient'", file=ROUT)
    print("x_lab = 'Rank of CpGs'", file=ROUT)
    print("", file=ROUT)

    print("included_char = 1", file=ROUT)
    print("included_col = 'grey'", file=ROUT)
    print("included_cex = 1", file=ROUT)
    print("included_cpgs = length(which(d_rank$Found=='True'))", file=ROUT)
    print("included_legend = paste('Found CpG (', included_cpgs, ')', sep='')",
          file=ROUT)
    print("", file=ROUT)

    print("missed_char = 13", file=ROUT)
    print("missed_col = 'red'", file=ROUT)
    print("missed_cex = 1", file=ROUT)
    print("missed_cpgs = length(which(d_rank$Found=='False'))", file=ROUT)
    print("missed_legend = paste('Missed CpG (', missed_cpgs, ')', sep='')",
          file=ROUT)
    print("", file=ROUT)

    print("lx = 1", file=ROUT)
    print("ly = max(d_rank$Coef)", file=ROUT)
    print("", file=ROUT)

    if outfile.endswith('.pdf'):
        print("pdf(file='%s', width=8, height=8)" % outfile, file=ROUT)
    else:
        print("png(file='%s', width=1080, height=1080)" % outfile, file=ROUT)
    print("plot(index, d_rank$Coef, ylab = y_lab, xlab = x_lab, col = \
         included_col, pch = included_char, cex = included_cex)", file=ROUT)
    print("points(index[d_rank$Found=='False'], \
          d_rank$Coef[d_rank$Found=='False'], col = missed_col, \
          pch = missed_char, cex = missed_cex)", file=ROUT)
    print("abline(h = 0,lty = 'dashed')", file=ROUT)
    print("legend(lx, ly,legend=c(included_legend, missed_legend), \
          col=c(included_col, missed_col), pch = c(included_char, \
          missed_char), pt.cex = c(included_cex, missed_cex))", file=ROUT)
    print("dev.off()", file=ROUT)
    ROUT.close()


def plot_corr(cage, dage, outfile, rfile):
    """
    Generate corrleation plot between the chronological age and epigenetic age.

    Parameters
    ----------
    cage : list
        List of chronological ages.
    dage : list
        List of predicted DNAm ages.
    outfile : str
        The name of output figure file. The suffix of this file must be \
        '.pdf' or '.png'.
    rfile : str
        The name of R script file. This file generates the above figure file.

    Returns
    -------
    None.
    """
    r, pval = stats.pearsonr(cage, dage)
    r = str(round(r, 3))
    pval = "{:.2e}".format(pval)
    ROUT = open(rfile, 'a')
    print("\n\n###### Rscript to generate scatter plot ######\n\n", file=ROUT)
    print("c_age <- c(%s)" % ','.join([str(i) for i in cage]), file=ROUT)
    print("d_age <- c(%s)" % ','.join([str(i) for i in dage]), file=ROUT)
    print("x_lab = 'Chronological age'", file=ROUT)
    print("y_lab = 'DNAm age'", file=ROUT)
    print("", file=ROUT)

    if outfile.endswith('.pdf'):
        print("pdf(file='%s', width=8, height=8)" % outfile, file=ROUT)
    else:
        print("png(file='%s', width=1080, height=1080)" % outfile, file=ROUT)

    print("plot(c_age, d_age, xlab=x_lab, ylab=y_lab, col='blue', sub=\"r = %s, P = %s\")" % (r, pval), file=ROUT)
    print("text(min(c_age) + 15, max(d_age), labels=c(\"r = %s, P = %s\"))" % (r, pval), file=ROUT)
    print("abline(a=0, b=1, lty='dashed', col='red', lwd=0.75)", file=ROUT)
    print("abline(glm(d_age~c_age), lwd=1, col='blue')", file=ROUT)
    #print("abline(v=seq(0,100,5),h=seq(0,100,5), lty='dashed', lwd=0.5, col='grey')", file=ROUT)
    print("dev.off()", file=ROUT)
    ROUT.close()


def config_log(switch, logfile=None):
    """
    Configureing the logging module.

    Parameters
    ----------
    switch : bool
        Debugging switch.
    Returns
    -------
    None.

    """
    if switch is True:
        if logfile is None:
            logging.basicConfig(
                format="%(asctime)s [%(levelname)s]  %(message)s",
                datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
        else:
            logging.basicConfig(
                filename=logfile,
                format="%(asctime)s [%(levelname)s]  %(message)s",
                datefmt='%Y-%m-%d %I:%M:%S', level=logging.DEBUG)
    else:
        if logfile is None:
            logging.basicConfig(
                format="%(asctime)s [%(levelname)s]  %(message)s",
                datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)
        else:
            logging.basicConfig(
                filename=logfile,
                format="%(asctime)s [%(levelname)s]  %(message)s",
                datefmt='%Y-%m-%d %I:%M:%S', level=logging.INFO)
