DICT_TRESHOLD = {33: 0, 34: 1, 35: 2, 36: 3, 37: 4, 38: 5, 39: 6, 40: 7, 41: 8, 42: 9, 43: 10,
                 44: 11, 45: 12, 46: 13, 47: 14, 48: 15, 49: 16, 50: 17, 51: 18, 52: 19, 53: 20,
                 54: 21, 55: 22, 56: 23, 57: 24, 58: 25, 59: 26, 60: 27, 61: 28, 62: 29, 63: 30,
                 64: 31, 65: 32, 66: 33, 67: 34, 68: 35, 69: 36, 70: 37, 71: 38, 72: 39, 73: 40}


def filter_gc(filtered_dict, key, min_gc_bound, max_gc_bound):
    """
        The function filters sequences by gc composition
            Parameters:
                filtered_dict - a dictionary with sequences that we filter
                key - current dictionary key
                min_gc_bound - lower bound gc_bounds
                max_gc_bound - upper limit gc_bounds
        Return:
            dictionary with only those sequences that have been filtered
    """
    gc_count = filtered_dict[key][0].count('G') + filtered_dict[key][0].count('C')+filtered_dict[key][0].count('g') + filtered_dict[key][0].count('c')
    gc_persent = (gc_count / len(filtered_dict[key][0]))*100
    if min_gc_bound > gc_persent  or gc_persent > max_gc_bound:
        del filtered_dict[key]
    return filtered_dict


def filter_length(filtered_dict, key, min_len_bound, max_len_bound):
    """
        The function filters sequences by length
            Parameters:
                filtered_dict - a dictionary with sequences that we filter
                key - current dictionary key
                min_len_bound - lower limit len_bounds
                max_len_bound - upper limit len_bounds
        Return:
           dictionary with only those sequences that have been filtered
       """
    len_seq = len(filtered_dict[key][0])
    if min_len_bound > len_seq or len_seq > max_len_bound:
        del filtered_dict[key]
    return filtered_dict


def filter_quality_threshold(filtered_dict, key, quality_threshold):
    """
        The function filters sequences by average read quality
            Parameters:
                filtered_dict - a dictionary with sequences that we filter
                key - current dictionary key
                quality_threshold - average quality threshold
        Return:
            dictionary with only those sequences that have been filtered
    """
    threshold = 0
    for nuc in filtered_dict[key][1]:
        threshold += int(DICT_TRESHOLD[ord(nuc)])
    mean_thresholds = threshold / len(filtered_dict[key][0])
    if mean_thresholds < quality_threshold:
        del filtered_dict[key]
    return filtered_dict
