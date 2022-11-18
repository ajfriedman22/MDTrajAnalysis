def error_bar(x1, x2, y1, y2, p, h, col):
    import matplotlib.pyplot as plt
    y = max([y1, y2])
    if p < 0.05 and p > 0.01:
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "*" , ha='center', va='bottom', color=col)
    if p < 0.01 and p > 0.001:
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "**" , ha='center', va='bottom', color=col)
    if p < 0.001:
        plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
        plt.text((x1+x2)*0.5, y+h, "***" , ha='center', va='bottom', color=col)

def plot_torsion(dihe_dist, name, i, maxima):
    import matplotlib.pyplot as plt

    #Seperate dihedral angles
    dihe_name = name[i]
    
    #Histogram of the data
    n, bins, patches = plt.hist(dihe_dist, 30, density=True, facecolor='g', alpha=0.75)
    #Inidcate Maxima
    plt.axvline(x = maxima[0], color = 'k')
    if len(maxima) == 2:
        plt.axvline(x = maxima[1], color = 'b')

    plt.xlabel('Torsional Angle(rad)')
    plt.ylabel('Probability')
    plt.xlim(-180, 180)
    plt.title('Histogram of Torsion Angle ' + dihe_name)
    plt.grid(True)
    plt.savefig('dihedrals/dihe_angle_' + dihe_name + '.png')
    plt.close()

