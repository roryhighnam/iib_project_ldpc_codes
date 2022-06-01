import matplotlib.pyplot as plt
import os
import csv
import numpy as np
from sympy import sec
from density_evolution import modified_density_evolution

base_dir = './random_ensemble_report_combined'
plot_base_dir = './random_ensemble_report_plots/'

colours = ['red', 'green', 'blue', 'orange']

def get_parameters(filename):
    # Get simulation parameters from filename
    simulation_parameters = filename.split('_')
    simulation_parameters_dict = {}
    for simulation_parameter in simulation_parameters:
        if '=' in simulation_parameter and '.csv' not in simulation_parameter:
            simulation_parameters_dict[simulation_parameter.split('=')[0]] = simulation_parameter.split('=')[1]
        elif 'combined.csv' in simulation_parameter:
            simulation_parameters_dict['combined'] = 'true'
        elif '.csv' in simulation_parameter:
            simulation_parameters_dict[simulation_parameter.split('=')[0]] = simulation_parameter.split('=')[1].split('.')[0]
    return simulation_parameters_dict

def plot_error_vs_iteration_number(erasure_prob, ns, figure, save_as_filename):
    colour_index = 0
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

            plt.plot(errors, label=simulation_parameters_dict['n'], color= colours[colour_index%len(colours)])
            # if simulation_parameters_dict['n'] == '100':
            #     if erasure_prob == 0.3:
            #         plt.hlines(0.00927103257501793, 0, 201, linestyles='dashed', colors=colours[colour_index%len(colours)], label='Finite length analysis n=100')
            #         plt.hlines(0.0007533780499659765, 0, 201, linestyles='dotted', colors=colours[colour_index%len(colours)], label='ML decoder n=100')
            #     if erasure_prob == 0.35:
            #         plt.hlines(0.0059553, 0, 201, linestyles='dotted', colors=colours[colour_index%len(colours)], label='ML decoder n=100')
            #         plt.hlines(0.0469878894379349, 0, 201, linestyles='dashed', colors=colours[colour_index%len(colours)], label='Finite length analysis n=100')
            #     if erasure_prob == 0.4:
            #         plt.hlines(0.145904820551055, 0, 201, linestyles='dashed', colors=colours[colour_index%len(colours)], label='Finite length analysis n=100')
            #         plt.hlines(0.045414847161572056, 0, 201, linestyles='dotted', colors=colours[colour_index%len(colours)], label='ML decoder n=100')

            # if simulation_parameters_dict['n'] == '1000':
            #     if erasure_prob == 0.3:
            #         plt.hlines(7.36669556389754e-6, 0, 201, linestyles='dashed', colors=colours[colour_index%len(colours)], label='Finite length analysis n=1000')
            #         # plt.hlines(1.28779716565522e-6, 0, 201, linestyles='dashed', colors=colours[colour_index%len(colours)], label='Finite length analysis n=1000')
            #         # ML decoder value needs changing
            #         plt.hlines(2.7454103745054328e-06, 0, 201, linestyles='dotted', colors=colours[colour_index%len(colours)], label='ML decoder n=1000')
            #     if erasure_prob == 0.35:
            #         plt.hlines(1.30239574884079e-5, 0, 201, linestyles='dashed', colors=colours[colour_index%len(colours)], label='Finite length analysis n=1000')
            #         plt.hlines(3.5772410602232966e-06, 0, 201, linestyles='dotted', colors=colours[colour_index%len(colours)], label='ML decoder n=1000')
            #     if erasure_prob == 0.4:
            #         plt.hlines(0.0214084026498356, 0, 201, linestyles='dashed', colors=colours[colour_index%len(colours)], label='Finite length analysis n=1000')
            #         # ML decoder updated but seems incorrect?
            #         plt.hlines(5.020577124995474e-06, 0, 201, linestyles='dotted', colors=colours[colour_index%len(colours)], label='ML decoder n=1000')

            colour_index += 1
            plt.yscale('log')
            plt.suptitle('Regular LDPC code simulation', fontsize='18')
            subtitle_text = ''
            for parameter in simulation_parameters_dict.keys():
                if parameter != 'n' and parameter != 'k' and parameter != 'num' and parameter != 'time':
                    subtitle_text += parameter
                    subtitle_text += '='
                    subtitle_text += simulation_parameters_dict[parameter]
                    subtitle_text += '  '
            plt.title(subtitle_text)
    plt.xlabel('Iterations')
    plt.ylabel('Log BER')
    plt.plot(modified_density_evolution(erasure_prob, int(simulation_parameters_dict['it']), int(simulation_parameters_dict['dv']), int(simulation_parameters_dict['dc']), min_error), '--', label='Density evolution')
    
    # Sort legend labels by ascending order
    handles, labels = plt.gca().get_legend_handles_labels()
    # text_titles = ['Density evolution', 'Finite length analysis n=100', 'Finite length analysis n=1000', 'ML decoder n=100',  'ML decoder n=1000']
    text_titles = ['Density evolution']
    text_labels = []
    text_handles = []
    text_indexes = []
    for text_title in text_titles:
        index = labels.index(text_title)
        text_indexes.append(index)
        text_labels.append(labels[index])
        text_handles.append(handles[index])
        handles.pop(index)
        labels.pop(index)
    labels = [int(label) for label in labels]
    zipped = list(zip(handles, labels))
    zipped.sort(key=lambda tup: tup[1])
    handles, labels = [[i for i,j in zipped], [j for i,j in zipped]]
    for i in range(len(labels)):
        labels[i] = '$n$ = ' + str(labels[i])
    for i in range(len(text_handles)):
        handles.append(text_handles[i])
        labels.append(text_labels[i])
    plt.legend(handles, labels, loc=(1.04,0))
    plt.grid('both')
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

def plot_error_vs_erasure_prob_fixed_code_combined(n, figure, save_as_filename):
    final_errors = {}
    # Plot simulations of error prob vs erasure prob
    plt.figure(figure)
    for filename in os.listdir(base_dir):

        simulation_parameters_dict = get_parameters(filename)

        if True:
            print(filename)
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
                    if parameter != 'k' and parameter != 'number' and parameter != 'BEC' and parameter != 'combined':
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

# plot_error_vs_erasure_prob_fixed_code_combined(512, figure, 'n=512_varying_erasure_prob.jpg')
# figure += 1
# plot_error_vs_iteration_number(0.3, [100,1000,10000], figure, 'BEC=0.3_n=10000_vs_iteration_number.jpg')
# figure += 1
plot_error_vs_iteration_number(0.42, [100000], figure, 'BEC=0.42_n=100000_vs_iteration_number_latest.jpg')
figure += 1

# plot_error_vs_iteration_number(0.42, [100,1000,10000], figure, 'BEC=0.42_n=10000_vs_iteration_number.jpg')
# figure += 1
# plot_error_vs_iteration_number(0.43, [100,1000,10000], figure, 'BEC=0.43_n=10000_vs_iteration_number.jpg')
# figure += 1
# plot_error_vs_iteration_number(0.5, [100,1000,10000], figure, 'BEC=0.5_n=10000_vs_iteration_number.jpg')
# figure += 1
# plot_error_vs_erasure_prob(200, figure, 'n=200_varying_erasure_prob.jpg')
# figure += 1