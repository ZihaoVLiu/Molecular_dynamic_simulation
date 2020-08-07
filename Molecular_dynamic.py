# The libraries you needed to implement the Python code
import numpy as np
import pandas as pd
import h5py
import matplotlib.pyplot as plt

# def get_count_fast(x, y, x_interval, y_interval):
#     result = np.zeros((x_interval.shape[0], y_interval.shape[0]))
#     x_count = np.zeros((x.shape[0], x_interval.shape[0], x_interval.shape[0]))
#     y_count = np.zeros((y.shape[0], y_interval.shape[0], y_interval.shape[0]))
#     # padding the first unit of x_count and y_count
#     x_count[:,0,0] = np.sum(x <= x_interval[0])
#     y_count[:,0,0] = np.sum(y <= y_interval[0])
#     for i in range(1, x_count.shape[1]):
#         count = (x_interval[i - 1] < x) & (x <= x_interval[i])
#         for ii in range(1, x_count.shape[2]):
#             x_count[:, i, ii] = count
#     for j in range(1, y_count.shape[1]):
#         count = (y_interval[j - 1] < y) & (y <= y_interval[j])
#         for jj in range(1, y_count.shape[2]):
#             y_count[:, jj, j] = count
#
#     for row in range(result.shape[0]):
#         for col in range(result.shape[1]):
#             result[row, col] = np.sum(x_count[row] * y_count[col])
#         print('{} row count done.'.format(str(row + 1)))
#     return result


def get_count(x, y, x_interval, y_interval):
    '''
    Get the number of atoms of each pixel.
    Note: The higher the resolution, the more computation time consuming
    :param x: The X-coord location of your view
    :param y: The Y-coord location of your view
    :param x_interval: The intervel of X axis (the resolution of your visualization)
    :param y_interval: The intervel of Y axis (the resolution of your visualization)
    :return: a n * n matrix (n is the resolution you want to visualize)
    '''
    result = np.zeros((x_interval.shape[0], y_interval.shape[0]))
    x_count = np.zeros((x_interval.shape[0], x.shape[0]))
    y_count = np.zeros((y_interval.shape[0], y.shape[0]))
    # padding the first unit of x_count and y_count
    x_count[0] = (x <= x_interval[0])
    y_count[0] = (y <= y_interval[0])
    for i in range(1, x_interval.size):
        x_count[i] = (x_interval[i - 1] < x) & (x <= x_interval[i])
    for j in range(1, y_interval.size):
        y_count[j] = (y_interval[j - 1] < y) & (y <= y_interval[j])
    for row in range(result.shape[0]):
        for col in range(result.shape[1]):
            result[row, col] = np.sum(x_count[row] * y_count[col])
        print('{} row count done.'.format(str(row + 1)))
    return result


def draw_contour(x, y, interval, cmap, title, save_name='_.png'):
    '''
    Visualize the molecular dynamic and output a .png format image.
    :param x: The X-coord location of your view.
    :param y: The y-coord location of your view.
    :param interval: The resolution you want (The higher the resolution, the more computation time consuming).
    :param cmap: Color map, the color you want to visualize on your heat map.
    :param title: Title showed on .png image.
    :param save_name: File name of your output image.
    :return: None, but save visualization image on your work path.
    '''
    x_max = np.max(x)
    x_min = np.min(x)
    y_max = np.max(y)
    y_min = np.min(y)
    x_interval = np.linspace(x_min, x_max, interval)
    y_interval = np.linspace(y_min, y_max, interval)
    X, Y = np.meshgrid(x_interval, y_interval)
    Z = get_count(x, y, x_interval, y_interval)
    plt.contourf(X, Y, Z, cmap=cmap)
    plt.colorbar()
    plt.xlabel('* 10^-1 nm')
    plt.ylabel('* 10^-1 nm')
    plt.title(title)
    plt.savefig(save_name, dpi=300)
    plt.show()


def get_atom_xyz(data, atom_type='all'):
    '''
    Get the particular atom you want from the given time
    :param data: data at different time
    :param atom_type: just the name or index of atom.
    :return: x, y, z location information
    '''
    if atom_type == 'all':
        return data[:, 2], data[:, 3], data[:, 4]
    elif atom_type == 'H' or atom_type == 1:
        indexes = np.where(data[:, 1] == 1)
    elif atom_type == 'O' or atom_type == 2:
        indexes = np.where(data[:, 1] == 2)
    elif atom_type == 'N' or atom_type == 3:
        indexes = np.where(data[:, 1] == 3)
    elif atom_type == 'Na' or atom_type == 4:
        indexes = np.where(data[:, 1] == 4)
    elif atom_type == 'Cl' or atom_type == 5:
        indexes = np.where(data[:, 1] == 5)
    else:
        print('Input atom_type invalid.')
        return
    result = data[indexes[0]]
    return result[:, 2], result[:, 3], result[:, 4]



if __name__ == '__main__':
    # load the .excel file and save it to .h5 format to save loading time
    path = '/Users/zihaoliu/Desktop/Pony/2molLNaCl.xlsx'

    # comment the following code after you run the following code block first time
    # this code block is used to transfer the .xlsx file to .h5 format,
    # which could save a lot of I/O time each time you run the code.
    data_frame01 = pd.read_excel(path, sheet_name='0.1ns')
    data_frame5 = pd.read_excel(path, sheet_name='5ns')
    data_frame10 = pd.read_excel(path, sheet_name='10ns')
    with h5py.File('data.h5', 'w') as hf:
        hf.create_dataset("01ns", data=data_frame01.iloc[:, 0:5].values)
        hf.create_dataset("5ns", data=data_frame5.iloc[:, 0:5].values)
        hf.create_dataset("10ns", data=data_frame10.iloc[:, 0:5].values)

    # load the h5 file you get from last step
    with h5py.File('data.h5', 'r') as hf:
        data01 = hf['01ns'][:]
        data5 = hf['5ns'][:]
        data10 = hf['10ns'][:]

    # Molecular dynamic simulation visualization
    # you can modify parameters which you want to try
    # load the xyz of the atom on particular time
    x, y, z = get_atom_xyz(data=data10, atom_type='N')
    # draw the heat map (also known as contour(等高线) map)
    draw_contour(x=x, y=z, interval=50, cmap='Oranges',
                 title='10 ns, N atom molecular dynamic (Front view)', save_name='10_N_Front.png')

