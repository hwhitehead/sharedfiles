def set_plot_defaults():
    """
    assign default plot settings
    """
    ## FIGURE
    #plt.rcParams["text.usetex"] = "True"

    ## FONT
    #plt.rcParams['font.serif']=['cm']
    #plt.rcParams['font.family']='serif'
    # plt.rcParams['font.serif']=['cm']
    plt.rcParams['text.latex.preamble']=r'\usepackage{amsmath}'

    plt.rcParams['font.size']=18
    plt.rcParams['xtick.labelsize']=15
    plt.rcParams['ytick.labelsize']=15
    plt.rcParams['legend.fontsize']=14
    plt.rcParams['axes.titlesize']=16
    plt.rcParams['axes.labelsize']=16
    plt.rcParams['axes.linewidth']=2
    plt.rcParams["lines.linewidth"] = 2.2
    ## TICKS
    plt.rcParams['xtick.top']='True'
    plt.rcParams['xtick.bottom']='True'
    plt.rcParams['xtick.minor.visible']='True'
    plt.rcParams['xtick.direction']='out'
    plt.rcParams['ytick.left']='True'
    plt.rcParams['ytick.right']='True'
    plt.rcParams['ytick.minor.visible']='True'
    plt.rcParams['ytick.direction']='out'
    plt.rcParams['xtick.major.width']=1.5
    plt.rcParams['xtick.minor.width']=1
    plt.rcParams['xtick.major.size']=4
    plt.rcParams['xtick.minor.size']=3
    plt.rcParams['ytick.major.width']=1.5
    plt.rcParams['ytick.minor.width']=1
    plt.rcParams['ytick.major.size']=4
    plt.rcParams['ytick.minor.size']=3
