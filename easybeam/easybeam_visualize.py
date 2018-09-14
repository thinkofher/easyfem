import matplotlib.pyplot as plt


def graph_moments(x, M, x_unit='m', y_unit='Nm'):

    fig = plt.figure(figsize=(10, 4), dpi=80, facecolor='w', edgecolor='k')

    plt.plot(x, M, color='red')
    plt.fill_between(x, M, color='red', alpha=0.3)
    plt.xlabel('x [{}]'.format(x_unit), fontsize=25)
    plt.ylabel('M [{}]'.format(y_unit), fontsize=25)
    plt.grid(True)
    plt.show()


def graph_shears(x, T, x_unit='m', y_unit='N'):

    fig = plt.figure(figsize=(10, 4), dpi=80, facecolor='w', edgecolor='k')

    plt.plot(x, T, color='green')
    plt.fill_between(x,T, color='green', alpha=0.3)
    plt.xlabel('x [{}]'.format(x_unit), fontsize=25)
    plt.ylabel('T [{}]'.format(y_unit), fontsize=25)
    plt.grid(True)
    plt.show()


def graph_disps(x, d, x_unit='m', y_unit='m'):

    fig = plt.figure(figsize=(10, 4), dpi=80, facecolor='w', edgecolor='k')

    plt.plot(x, d, color='purple')
    plt.fill_between(x, d, color='purple', alpha=0.3)
    plt.xlabel('x [{}]'.format(x_unit), fontsize=25)
    plt.ylabel('d [{}]'.format(y_unit), fontsize=25)
    plt.grid(True)
    plt.show()


def graph_rots(x, r, x_unit='m', y_unit='rad'):

    fig = plt.figure(figsize=(10, 4), dpi=80, facecolor='w', edgecolor='k')

    plt.plot(x, r, color='orange')
    plt.fill_between(x, r, color='orange', alpha=0.3)
    plt.xlabel('x [{}]'.format(x_unit), fontsize=25)
    plt.ylabel('r [{}]'.format(y_unit), fontsize=25)
    plt.grid(True)
    plt.show()
