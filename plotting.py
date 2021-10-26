import matplotlib.pyplot as plt
import os
import csv
import numpy as np

base_dir = './simulation_data'
plot_base_dir = './simulation_plots/'

def get_parameters(filenme):
    # Get simulation parameters from filename
    simulation_parameters = filename.split('_')
    simulation_parameters_dict = {}
    for simulation_parameter in simulation_parameters:
        if '=' in simulation_parameter and '.csv' not in simulation_parameter:
            simulation_parameters_dict[simulation_parameter.split('=')[0]] = simulation_parameter.split('=')[1]
        elif '.csv' in simulation_parameter:
            simulation_parameters_dict[simulation_parameter.split('=')[0]] = simulation_parameter.split('=')[1].split('.')[0]
    return simulation_parameters_dict

# Plot simulations with erausre_prob = 0.3 for varying n
plt.figure(1)
for filename in os.listdir(base_dir):

    simulation_parameters_dict = get_parameters(filename)

    if simulation_parameters_dict['BEC'] == '0.42':
        errors = []

        with open(base_dir + '/' + filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                if len(row) == 1:
                    errors += row
        errors = [float(error) for error in errors]

        plt.plot(errors, label=simulation_parameters_dict['n'])
        plt.yscale('log')
        plt.suptitle('Regular LDPC code simulation', fontsize='18')
        subtitle_text = ''
        for parameter in simulation_parameters_dict.keys():
            if parameter != 'n' and parameter != 'k':
                subtitle_text += parameter
                subtitle_text += '='
                subtitle_text += simulation_parameters_dict[parameter]
                subtitle_text += '  '
        plt.title(subtitle_text)
plt.xlabel('Iterations')
plt.ylabel('Log BER')

# Sort legend labels by ascending order
handles, labels = plt.gca().get_legend_handles_labels()
labels = [int(label) for label in labels]
zipped = list(zip(handles, labels))
zipped.sort(key=lambda tup: tup[1])
print(zipped)
handles, labels = [[i for i,j in zipped], [j for i,j in zipped]]
plt.legend(handles, labels)

plt.savefig(plot_base_dir + 'BER=0.42_varying_n.jpg', bbox_inches='tight')

# Plot simulations with erausre_prob = 0.4 for varying n
plt.figure(2)
for filename in os.listdir(base_dir):

    simulation_parameters_dict = get_parameters(filename)

    if simulation_parameters_dict['BEC'] == '0.43':
        errors = []

        with open(base_dir + '/' + filename) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            for row in csv_reader:
                if len(row) == 1:
                    errors += row
        errors = [float(error) for error in errors]

        plt.plot(errors, label=simulation_parameters_dict['n'])
        plt.yscale('log')
        plt.suptitle('Regular LDPC code simulation', fontsize='18')
        subtitle_text = ''
        for parameter in simulation_parameters_dict.keys():
            if parameter != 'n' and parameter != 'k':
                subtitle_text += parameter
                subtitle_text += '='
                subtitle_text += simulation_parameters_dict[parameter]
                subtitle_text += '  '
        plt.title(subtitle_text)
plt.xlabel('Iterations')
plt.ylabel('Log BER')

# Sort legend labels by ascending order
handles, labels = plt.gca().get_legend_handles_labels()
labels = [int(label) for label in labels]
zipped = list(zip(handles, labels))
zipped.sort(key=lambda tup: tup[1])
print(zipped)
handles, labels = [[i for i,j in zipped], [j for i,j in zipped]]
plt.legend(handles, labels)

plt.savefig(plot_base_dir + 'BER=0.43_varying_n.jpg', bbox_inches='tight')

# Plot final error probability against code length, n for constant erasure prob=0.4
# plt.figure(3)
# n_list = []
# error_list = []
# for filename in os.listdir(base_dir):

#     simulation_parameters_dict = get_parameters(filename)

#     if simulation_parameters_dict['BEC'] == '0.4':

#         with open(base_dir + '/' + filename) as csv_file:
#             csv_reader = csv.reader(csv_file, delimiter=',')
#             for row in csv_reader:
#                 if csv_reader.line_num == 51:
#                     error_list.append(float(row[0]))
#                     n_list.append(int(simulation_parameters_dict['n']))

# zipped = list(zip(n_list, error_list))
# zipped.sort(key=lambda tup: tup[0])
# n_list, error_list = [[i for i,j in zipped], [j for i,j in zipped]]

# plt.plot(n_list[:20], error_list[:20], label='Message passing')
# plt.yscale('log')
# plt.suptitle('Regular LDPC code simulation', fontsize='18')
# subtitle_text = ''
# for parameter in simulation_parameters_dict.keys():
#     if parameter != 'n' and parameter != 'k':
#         subtitle_text += parameter
#         subtitle_text += '='
#         subtitle_text += simulation_parameters_dict[parameter]
#         subtitle_text += '  '
# plt.title(subtitle_text)
# plt.xlabel('n')
# plt.ylabel('Log BER')

# plt.savefig(plot_base_dir + 'BER=0.4_varying_n.jpg', bbox_inches='tight')
