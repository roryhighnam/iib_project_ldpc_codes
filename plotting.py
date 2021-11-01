import matplotlib.pyplot as plt
import os
import csv
import numpy as np
from density_evolution import density_evolution

base_dir = './simulation_data'
plot_base_dir = './simulation_plots/'

def get_parameters(filename):
    # Get simulation parameters from filename
    simulation_parameters = filename.split('_')
    simulation_parameters_dict = {}
    for simulation_parameter in simulation_parameters:
        if '=' in simulation_parameter and '.csv' not in simulation_parameter:
            simulation_parameters_dict[simulation_parameter.split('=')[0]] = simulation_parameter.split('=')[1]
        elif '.csv' in simulation_parameter:
            simulation_parameters_dict[simulation_parameter.split('=')[0]] = simulation_parameter.split('=')[1].split('.')[0]
    return simulation_parameters_dict

def plot_error_vs_iteration_number(erasure_prob, ns, figure, save_as_filename):
    min_error = erasure_prob
    # Plot simulations with certain erausre_prob for various values of n
    plt.figure(figure)
    for filename in os.listdir(base_dir):

        simulation_parameters_dict = get_parameters(filename)

        if simulation_parameters_dict['BEC'] == str(erasure_prob) and int(simulation_parameters_dict['n']) in ns and 'number' not in simulation_parameters_dict.keys():
            errors = [erasure_prob]

            with open(base_dir + '/' + filename) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    if len(row) == 1:
                        errors += row
                        if float(row[0]) < min_error:
                            min_error = float(row[0])
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
    plt.plot(density_evolution(erasure_prob, int(simulation_parameters_dict['it']), int(simulation_parameters_dict['dv']), int(simulation_parameters_dict['dc']), min_error), '--', label='Density evolution')


    # Sort legend labels by ascending order
    handles, labels = plt.gca().get_legend_handles_labels()
    i = labels.index('Density evolution')
    density_handle = handles[i]
    density_label = labels[i]
    handles.pop(i)
    labels.pop(i)
    labels = [int(label) for label in labels]
    zipped = list(zip(handles, labels))
    zipped.sort(key=lambda tup: tup[1])
    handles, labels = [[i for i,j in zipped], [j for i,j in zipped]]
    handles.append(density_handle)
    labels.append(density_label)
    plt.legend(handles, labels)

    plt.savefig(plot_base_dir + save_as_filename, bbox_inches='tight')

def plot_error_vs_n(erasure_prob, figure, save_as_filename):

    # Plot simulations with certain erausre_prob against n for fixed code parameters, also plot ML decoder
    final_errors = []
    ns = []
    optimal_final_errors = []
    optimal_ns = []

    subtitle_text = ''
    plt.figure(figure)
    for filename in os.listdir(base_dir):

        simulation_parameters_dict = get_parameters(filename)

        if simulation_parameters_dict['BEC'] == str(erasure_prob) and 'number' not in simulation_parameters_dict.keys():
            if subtitle_text == '':
                for parameter in simulation_parameters_dict.keys():
                    if parameter != 'n' and parameter != 'k' and parameter != 'number':
                        subtitle_text += parameter
                        subtitle_text += '='
                        subtitle_text += simulation_parameters_dict[parameter]
                        subtitle_text += '  '

            with open(base_dir + '/' + filename) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    if csv_reader.line_num == 51:
                        final_errors.append(float(row[0]))
                        ns.append(int(simulation_parameters_dict['n']))
                    elif csv_reader.line_num == 55:
                        # Stored as a percentage so divide by 100
                        optimal_final_errors.append(float(row[1])/100)
                        optimal_ns.append(int(simulation_parameters_dict['n']))

            # plt.plot(errors, label=simulation_parameters_dict['n'])
    zipped = list(zip(ns, final_errors))
    zipped.sort(key=lambda tup: tup[0])
    ns, final_errors = [[i for i,j in zipped], [j for i,j in zipped]]
    
    zipped = list(zip(optimal_ns, optimal_final_errors))
    zipped.sort(key=lambda tup: tup[0])
    optimal_ns, optimal_final_errors = [[i for i,j in zipped], [j for i,j in zipped]]

    plt.plot(ns, final_errors, label='Iterative decoder')
    plt.plot(optimal_ns, optimal_final_errors, label='ML decoder')
    plt.yscale('log')
    plt.legend()
    plt.suptitle('Regular LDPC code simulation', fontsize='18')
    plt.title(subtitle_text)
    plt.xlabel('n')
    plt.ylabel('Log BER')

    plt.savefig(plot_base_dir + save_as_filename, bbox_inches='tight')

def plot_concentration(erasure_prob, n, figure, save_as_filename):
    count = 0

    # Plot simulations for many different codes all with the same parameters 
    plt.figure(figure)
    for filename in os.listdir(base_dir):

        simulation_parameters_dict = get_parameters(filename)

        if simulation_parameters_dict['BEC'] == str(erasure_prob) and int(simulation_parameters_dict['n']) == n and 'number' in simulation_parameters_dict.keys():
            errors = []
            count += 1

            with open(base_dir + '/' + filename) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    if len(row) == 1:
                        errors += row
            errors = [float(error) for error in errors]

            plt.plot(errors)
            plt.yscale('log')
            plt.suptitle('Regular LDPC code simulation', fontsize='18')
            subtitle_text = ''
            for parameter in simulation_parameters_dict.keys():
                if parameter != 'k' and parameter != 'number':
                    subtitle_text += parameter
                    subtitle_text += '='
                    subtitle_text += simulation_parameters_dict[parameter]
                    subtitle_text += '  '
    subtitle_text += 'count='+str(count)
    plt.title(subtitle_text)
    plt.xlabel('Iterations')
    plt.ylabel('Log BER')

    plt.savefig(plot_base_dir + save_as_filename, bbox_inches='tight')

def plot_error_vs_erasure_prob_fixed_code(n, figure, save_as_filename):
    final_errors = {}
    # Plot simulations of error prob vs erasure prob
    plt.figure(figure)
    for filename in os.listdir(base_dir):

        simulation_parameters_dict = get_parameters(filename)

        if simulation_parameters_dict['n'] == str(n) and 'number' in simulation_parameters_dict.keys():
            # final_errors = []

            with open(base_dir + '/' + filename) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    if csv_reader.line_num == 51:
                        if simulation_parameters_dict['number'] in final_errors.keys():
                            final_errors[simulation_parameters_dict['number']].append((float(simulation_parameters_dict['BEC']), float(row[0])))
                        else:
                            final_errors[simulation_parameters_dict['number']] = [(float(simulation_parameters_dict['BEC']), float(row[0]))]

            # plt.plot(errors)
            plt.yscale('log')
            plt.suptitle('Regular LDPC code simulation', fontsize='18')
            subtitle_text = ''
            for parameter in simulation_parameters_dict.keys():
                if parameter != 'k' and parameter != 'number' and parameter != 'BEC':
                    subtitle_text += parameter
                    subtitle_text += '='
                    subtitle_text += simulation_parameters_dict[parameter]
                    subtitle_text += '  '
    subtitle_text += 'count='+str(len(final_errors.keys()))
    for key,value in final_errors.items():
        value.sort(key=lambda tup: tup[0])
        ns, final_errors = [[i for i,j in value], [j for i,j in value]]
        plt.plot(ns, final_errors)
    plt.title(subtitle_text)
    plt.xlabel('Erasure Probability')
    plt.ylabel('Log BER')

    plt.savefig(plot_base_dir + save_as_filename, bbox_inches='tight')

def plot_error_vs_erasure_prob(n, figure, save_as_filename):
    # Plot simulations with certain erausre_prob against n for fixed code parameters
    final_errors = []
    erasure_probs = []
    optimal_final_errors = []
    optimal_erasure_probs = []
    subtitle_text = ''
    plt.figure(figure)
    for filename in os.listdir(base_dir):

        simulation_parameters_dict = get_parameters(filename)

        if simulation_parameters_dict['n'] == str(n) and 'number' not in simulation_parameters_dict.keys():
            if subtitle_text == '':
                for parameter in simulation_parameters_dict.keys():
                    if parameter != 'k' and parameter != 'number' and parameter != 'BEC':
                        subtitle_text += parameter
                        subtitle_text += '='
                        subtitle_text += simulation_parameters_dict[parameter]
                        subtitle_text += '  '

            with open(base_dir + '/' + filename) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                for row in csv_reader:
                    if csv_reader.line_num == 51:
                        final_errors.append(float(row[0]))
                        erasure_probs.append(float(simulation_parameters_dict['BEC']))
                    elif csv_reader.line_num == 55:
                        optimal_final_errors.append(float(row[1])/100)
                        optimal_erasure_probs.append(float(simulation_parameters_dict['BEC']))

    zipped = list(zip(erasure_probs, final_errors))
    zipped.sort(key=lambda tup: tup[0])
    erasure_probs, final_errors = [[i for i,j in zipped], [j for i,j in zipped]]
    zipped = list(zip(optimal_erasure_probs, optimal_final_errors))
    zipped.sort(key=lambda tup: tup[0])
    optimal_erasure_probs, optimal_final_errors = [[i for i,j in zipped], [j for i,j in zipped]]
    plt.plot(erasure_probs, final_errors, label='Iterative decoder')
    plt.plot(optimal_erasure_probs, optimal_final_errors, label='ML decoder')
    plt.yscale('log')
    plt.legend()
    plt.suptitle('Regular LDPC code simulation', fontsize='18')
    plt.title(subtitle_text)
    plt.xlabel('Erasure Prob')
    plt.ylabel('Log BER')

    plt.savefig(plot_base_dir + save_as_filename, bbox_inches='tight')

# erasure_probs = [0.3, 0.4, 0.42, 0.43]
# ns = [50,100,200,500,1000]
figure = 1
# for erasure_prob in erasure_probs:
#     plot_error_vs_iteration_number(erasure_prob, ns, figure, 'BEC='+str(erasure_prob)+'_varying_n.jpg')
#     figure += 1

# plot_error_vs_n(0.4, figure, 'BEC=0.4_error_against_n.jpg')
# figure += 1

# erasure_probs = [0.3, 0.4, 0.42]
# for erasure_prob in erasure_probs:
#     plot_concentration(erasure_prob, 200, figure, 'BEC='+str(erasure_prob)+'_concentration.jpg')
#     figure += 1

plot_error_vs_erasure_prob_fixed_code(512, figure, 'n=512_varying_erasure_prob.jpg')
figure += 1

# plot_error_vs_erasure_prob(200, figure, 'n=200_varying_erasure_prob.jpg')
# figure += 1