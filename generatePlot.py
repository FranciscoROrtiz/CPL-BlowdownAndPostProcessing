import matplotlib.pyplot as plt


def generatePlot(dataset, xValues, Info):
    """
    Plots
    Args:
        dataset (numpy array):
        xValues (numpy array):
        Info (list):
    """
    Input = Info[0]
    Output = Info[1]
    UnitConversion = Info[2]

    plt.plot(xValues, dataset*UnitConversion)
    plt.xlabel(f"{Input}")
    plt.ylabel(f"{Output}")
    plt.title(f"{Output} vs. {Input}")
    plt.grid()
    # show the figure
    plt.show()
