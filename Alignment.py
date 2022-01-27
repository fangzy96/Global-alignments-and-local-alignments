import numpy as np
import argparse
import copy

def solution_construction(align_type, decision_array, current_alignment, sequence_x, sequence_y, len_y, len_x):

    if align_type == 'g':
        if len_x == 0 and len_y == 0:
            solution_array.append(current_alignment)
        else:
            # renew the alignment
            for decision in decision_array[len_y][len_x]:
                if decision == 0:
                    renew_array = copy.deepcopy(current_alignment)
                    renew_array[0] = renew_array[0] + sequence_x[len_x - 1]
                    renew_array[1] = renew_array[1] + '-'
                    solution_construction(align_type, decision_array, renew_array, sequence_x, sequence_y, len_y, len_x - 1)
                if decision == 1:
                    renew_array = copy.deepcopy(current_alignment)
                    renew_array[0] = renew_array[0] + '-'
                    renew_array[1] = renew_array[1] + sequence_y[len_y-1]
                    solution_construction(align_type, decision_array, renew_array, sequence_x, sequence_y, len_y - 1, len_x)
                if decision == 2:
                    renew_array = copy.deepcopy(current_alignment)
                    renew_array[0] = renew_array[0] + sequence_x[len_x-1]
                    renew_array[1] = renew_array[1] + sequence_y[len_y-1]
                    solution_construction(align_type, decision_array, renew_array, sequence_x, sequence_y, len_y - 1, len_x - 1)

    elif align_type == 'l':
        if 3 in decision_array[len_y][len_x]:
            solution_array.append(current_alignment)
        else:
            # renew the alignment
            for decision in decision_array[len_y][len_x]:
                if decision == 0:
                    renew_array = copy.deepcopy(current_alignment)
                    renew_array[0] = renew_array[0] + sequence_x[len_x - 1]
                    renew_array[1] = renew_array[1] + '-'
                    solution_construction(align_type, decision_array, renew_array, sequence_x, sequence_y, len_y, len_x - 1)
                if decision == 1:
                    renew_array = copy.deepcopy(current_alignment)
                    renew_array[0] = renew_array[0] + '-'
                    renew_array[1] = renew_array[1] + sequence_y[len_y-1]
                    solution_construction(align_type, decision_array, renew_array, sequence_x, sequence_y, len_y - 1, len_x)
                if decision == 2:
                    renew_array = copy.deepcopy(current_alignment)
                    renew_array[0] = renew_array[0] + sequence_x[len_x-1]
                    renew_array[1] = renew_array[1] + sequence_y[len_y-1]
                    solution_construction(align_type, decision_array, renew_array, sequence_x, sequence_y, len_y - 1, len_x - 1)

def global_alignment(score_table, sequence_x, sequence_y, len_x, len_y):
    match = score_table[0]
    mismatch = score_table[1]
    gap = score_table[2]

    # Initialization
    memory_array = np.zeros((len_y + 1, len_x + 1))
    decision_array = [[[]for i in range(len_x + 1)] for j in range(len_y + 1)]
    # decision_array = np.array([[[]]])
    # decision_array.reshape(len_sequence_y + 1, len_sequence_x + 1, 1)
    for i in range(len_x + 1):
        memory_array[0][i] = -i * score_table[2]
        decision_array[0][i] = [0]
        # np.append(decision_array[0][i], [0])
    for j in range(len_y + 1):
        memory_array[j][0] = -j * score_table[2]
        # np.append(decision_array[j][0], [1])
        decision_array[j][0] = [1]

    # Calculate the score
    for i in range(1, len_y + 1):
        for j in range(1, len_x + 1):
            left = memory_array[i][j - 1] - gap
            up = memory_array[i - 1][j] - gap
            if sequence_x[j - 1] == sequence_y[i - 1]:
                diag = memory_array[i - 1][j - 1] + match
            else:
                diag = memory_array[i - 1][j - 1] - mismatch
            current_score = [left, up, diag]
            current_max = max(current_score)
            # Left: 0, Up: 1, Diagonal: 2
            for k in range(3):
                if current_max == current_score[k]:
                    decision_array[i][j].append(k)
            memory_array[i][j] = current_max
    # print(decision_array)
    return memory_array[len_y][len_x], decision_array



def local_alignment(score_table, sequence_x, sequence_y, len_x, len_y):
    match = score_table[0]
    mismatch = score_table[1]
    gap = score_table[2]
    temp_score = 0
    current_state = [[0, 0]]

    # Initialization
    memory_array = np.zeros((len_y + 1, len_x + 1))
    decision_array = [[[]for i in range(len_x + 1)] for j in range(len_y + 1)]
    for i in range(len_x + 1):
        memory_array[0][i] = 0
        decision_array[0][i] = [3]
    for j in range(len_y + 1):
        memory_array[j][0] = 0
        decision_array[j][0] = [3]

    for i in range(1, len_y + 1):
        for j in range(1, len_x + 1):
            left = memory_array[i][j - 1] - gap
            up = memory_array[i - 1][j] - gap
            if sequence_x[j - 1] == sequence_y[i - 1]:
                diag = memory_array[i - 1][j - 1] + match
            else:
                diag = memory_array[i - 1][j - 1] - mismatch
            current_score = [left, up, diag, 0]
            current_max = max(current_score)
            # Left: 0, Up: 1, Diagonal: 2, Self: 3
            for k in range(4):
                if current_max == current_score[k]:
                    decision_array[i][j].append(k)
            memory_array[i][j] = current_max
            if memory_array[i][j] > temp_score:
                temp_score = memory_array[i][j]
                current_state = [[i, j]]
            elif memory_array[i][j] == temp_score:
                current_state.append([i, j])

    return temp_score, decision_array, current_state

if __name__ == '__main__':
    # Set up argparse arguments
    parser = argparse.ArgumentParser(description='Run the alignment algorithm.')
    parser.add_argument('path', metavar='PATH', type=str, help='The path to the data.')
    args = parser.parse_args()
    data_path = args.path

    data = open(data_path)
    data_line = data.readlines()  # Get each data line
    align_type = data_line[0].rstrip('\n')  # Get alignment type
    score_table = data_line[1].rstrip('\n')  # Get score table
    sequence_x = data_line[2].rstrip('\n')  # Get sequence x
    sequence_y = data_line[3].rstrip('\n')  # Get sequence y

    score_table = score_table.split(' ')
    score_table = list(map(int, score_table))  # Convert str to int
    len_x = len(sequence_x)
    len_y = len(sequence_y)
    solution_array = []

    # Global alignment
    if align_type == 'g':
        score, decision_array = global_alignment(score_table, sequence_x, sequence_y, len_x, len_y)
        current_alignment = ['', '']
        solution_construction(align_type, decision_array, current_alignment, sequence_x, sequence_y, len_y, len_x)
        for k in range(len(solution_array)):
            solution_array[k][0] = solution_array[k][0][::-1]
            solution_array[k][1] = solution_array[k][1][::-1]
        for i, j in solution_array:
            print('S1:', i)
            print('S2:', j)
        print('Score', score)
        print('Number of solutions:', len(solution_array))
        # print(solution_array)

    # Local alignment
    elif align_type == 'l':
        score, decision_array, current_state = local_alignment(score_table, sequence_x, sequence_y, len_x, len_y)
        for state in current_state:
            current_alignment = ['', '']
            solution_construction(align_type, decision_array, current_alignment, sequence_x, sequence_y, state[0], state[1])
        for k in range(len(solution_array)):
            solution_array[k][0] = solution_array[k][0][::-1]
            solution_array[k][1] = solution_array[k][1][::-1]
        for i, j in solution_array:
            print('S1:', i)
            print('S2:', j)
        print('Score', score)
        print('Number of solutions:', len(solution_array))
