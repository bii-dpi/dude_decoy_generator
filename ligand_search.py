def search_in_range(value, delta, sorted_arr):
    selected_indices = []
    min_value = value - delta
    max_value = value + delta
    for i, other_value in enumerate(sorted_arr):
        if other_value > max_value:
            break
        elif other_value < min_value:
            continue
        selected_indices.append(i)

    return set(selected_indices)


search_in_range(10, 5, list(range(1, 10000)))



