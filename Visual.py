import matplotlib.pyplot as plt
from NodalAnalysis import intersection_point

def nodal_graph(irp, vlp):
    operating_point = intersection_point(irp, vlp)
    flowrate = operating_point[0]
    pressure = operating_point[1]

    line_width = 3.5
    fig, axs = plt.subplots(figsize=(12,8))

    text_size = 10
    plt.rc('font', size=text_size)          
    plt.rc('axes', titlesize=text_size)   
    plt.rc('axes', labelsize=text_size)    
    plt.rc('xtick', labelsize=text_size)    
    plt.rc('ytick', labelsize=text_size)   
    plt.rc('legend', fontsize=text_size)  
    plt.rc('figure', titlesize=text_size) 
    plt.xlabel('xlabel', fontsize=text_size)
    plt.ylabel('ylabel', fontsize=text_size)

    axs.plot(irp[0], irp[1], linewidth = line_width, label='IRP curve')
    axs.plot(vlp[0], vlp[1], color = 'r', linewidth =line_width, label='VLP curve')
    axs.scatter(flowrate, pressure, s=80, c='black', 
                label = 'Operating point:\nQ = {:0.0f} bbl/day\nP = {:0.0f} psi'.format(flowrate, pressure))
    
    axs.set_title(label = 'Nodal Analysis', fontsize = 18, fontweight="bold")    
    axs.legend(loc='upper right')
    axs.set(ylabel='Bottomhole pressure, psi', xlabel = 'Flowrate, bbl/day')
    axs.set_ylim([0, 1.1 * max(irp[1])])
    axs.set_xlim([0, 1.1 * max(irp[0])])
    axs.grid(True)

    plt.show()

